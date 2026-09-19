use crate::bed::BedRecord;
use std::collections::HashSet;

/// Stage 1: bedmap --max-element --fraction-either 0.1
///
/// For each reference interval, finds all overlapping intervals (>=10% overlap
/// of either element) and maps the reference interval to the highest-scoring element
/// among those overlapping intervals. Unlike transitive clustering, this performs
/// local non-maximum suppression, ensuring that non-overlapping intervals linked
/// only through an intermediate lower-scoring element are not dropped.
/// Input MUST be sorted by chrom + start + end.
pub fn bedmap_max_element(records: &[BedRecord]) -> Vec<BedRecord> {
    if records.is_empty() {
        return vec![];
    }

    let mut result: Vec<BedRecord> = Vec::with_capacity(records.len());
    let mut chrom_start = 0;

    while chrom_start < records.len() {
        let current_chrom = &records[chrom_start].chrom;
        let mut chrom_end = chrom_start + 1;
        while chrom_end < records.len() && &records[chrom_end].chrom == current_chrom {
            chrom_end += 1;
        }

        let chrom_records = &records[chrom_start..chrom_end];
        let max_len = chrom_records.iter().map(|r| r.length()).max().unwrap_or(0);
        let n = chrom_records.len();

        for i in 0..n {
            let rec_i = &chrom_records[i];
            let lower_bound = rec_i.start.saturating_sub(max_len);
            let left = chrom_records[..i].partition_point(|r| r.start < lower_bound);
            let right = i + chrom_records[i..].partition_point(|r| r.start < rec_i.end);

            let mut best = rec_i;
            for cand in &chrom_records[left..right] {
                if has_reciprocal_overlap(rec_i, cand, 0.1) {
                    let is_better = cand.score > best.score
                        || ((cand.score - best.score).abs() < 1e-6
                            && (cand.start > best.start
                                || (cand.start == best.start && cand.end > best.end)));
                    if is_better {
                        best = cand;
                    }
                }
            }
            result.push(best.clone());
        }

        chrom_start = chrom_end;
    }

    result
}

/// Check if two intervals overlap by at least `frac` of *either* element,
/// matching `bedmap --fraction-either` (an OR, not a reciprocal AND).
fn has_reciprocal_overlap(a: &BedRecord, b: &BedRecord, frac: f64) -> bool {
    let overlap_start = a.start.max(b.start);
    let overlap_end = a.end.min(b.end);
    if overlap_start >= overlap_end {
        return false;
    }
    let overlap = (overlap_end - overlap_start) as f64;
    let len_a = a.length() as f64;
    let len_b = b.length() as f64;
    if len_a == 0.0 || len_b == 0.0 {
        return false;
    }
    overlap / len_a >= frac || overlap / len_b >= frac
}

/// Stage 2: exact deduplication.
/// Remove records with identical BED9 string representation.
pub fn exact_dedup(records: &[BedRecord]) -> Vec<BedRecord> {
    let mut seen = HashSet::new();
    let mut result = Vec::new();
    for rec in records {
        let key = rec.to_bed9_string();
        if seen.insert(key) {
            result.push(rec.clone());
        }
    }
    result
}

/// Stage 3: near-dedup (from overlap_filter.py).
///
/// Sequential pass: if start or end is within +/-10 bp of the previous kept
/// record, keep the one with the higher score (tie-break by length).
///
/// The Python original carried its `prev_*` state across the whole file, which
/// was safe only because every input file held one chromosome. This port takes
/// multi-sequence FASTA, so the comparison is reset at each sequence boundary;
/// otherwise the first record of a sequence could be dropped because of where
/// the *previous* sequence happened to end.
pub fn near_dedup(records: &[BedRecord]) -> Vec<BedRecord> {
    let mut result: Vec<BedRecord> = Vec::new();
    // chrom, start, end, score, length
    let mut prev: Option<(&str, i64, i64, f64, i64)> = None;

    for rec in records {
        let start = rec.start as i64;
        let end = rec.end as i64;
        let length = end - start;

        if let Some((prev_chrom, prev_start, prev_end, prev_score, prev_length)) = prev {
            if prev_chrom == rec.chrom
                && ((start - prev_start).abs() < 10 || (end - prev_end).abs() < 10)
            {
                if (prev_score - rec.score).abs() < 0.01 {
                    // Same score: keep the longer interval.
                    if length > prev_length {
                        result.pop();
                    } else {
                        continue;
                    }
                } else if rec.score > prev_score {
                    result.pop();
                } else {
                    continue;
                }
            }
        }

        result.push(rec.clone());
        prev = Some((&rec.chrom, start, end, rec.score, length));
    }

    result
}

/// Run all 3 filtering stages in sequence.
pub fn filter_overlaps(records: &[BedRecord]) -> Vec<BedRecord> {
    let stage1 = bedmap_max_element(records);
    let stage2 = exact_dedup(&stage1);
    near_dedup(&stage2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bed::test_record as make_rec;

    #[test]
    fn test_bedmap_max_non_overlapping() {
        let recs = vec![
            make_rec("chr1", 0, 100, "a", 50.0),
            make_rec("chr1", 200, 300, "b", 60.0),
        ];
        let result = bedmap_max_element(&recs);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_bedmap_max_overlapping() {
        let recs = vec![
            make_rec("chr1", 0, 100, "a", 50.0),
            make_rec("chr1", 10, 110, "b", 80.0),
        ];
        let stage1 = bedmap_max_element(&recs);
        assert_eq!(stage1.len(), 2);
        assert_eq!(stage1[0].name, "b");
        assert_eq!(stage1[1].name, "b");

        let result = filter_overlaps(&recs);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].name, "b");
    }

    #[test]
    fn test_bedmap_reviewer_counterexample() {
        // Reviewer counterexample:
        // A: [0, 171), score 200
        // B: [150, 321), score 150
        // C: [300, 471), score 190
        // Transitive clustering would incorrectly merge A, B, C into one cluster and drop C.
        // True BEDOPS bedmap --max-element + dedup outputs both A and C, dropping only B.
        let recs = vec![
            make_rec("chr1", 0, 171, "A", 200.0),
            make_rec("chr1", 150, 321, "B", 150.0),
            make_rec("chr1", 300, 471, "C", 190.0),
        ];
        let stage1 = bedmap_max_element(&recs);
        assert_eq!(stage1.len(), 3);
        assert_eq!(stage1[0].name, "A");
        assert_eq!(stage1[1].name, "A");
        assert_eq!(stage1[2].name, "C");

        let result = filter_overlaps(&recs);
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].name, "A");
        assert_eq!(result[1].name, "C");
    }

    #[test]
    fn test_exact_dedup() {
        let recs = vec![
            make_rec("chr1", 0, 100, "a", 50.0),
            make_rec("chr1", 0, 100, "a", 50.0),
            make_rec("chr1", 200, 300, "b", 60.0),
        ];
        let result = exact_dedup(&recs);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_near_dedup_keeps_higher_score() {
        let recs = vec![
            make_rec("chr1", 100, 270, "a", 50.0),
            make_rec("chr1", 105, 272, "b", 80.0), // within +/-10bp
        ];
        let result = near_dedup(&recs);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].name, "b");
    }

    #[test]
    fn test_near_dedup_same_score_keeps_longer() {
        let recs = vec![
            make_rec("chr1", 100, 270, "a", 50.0),
            make_rec("chr1", 105, 280, "b", 50.0), // same score, longer
        ];
        let result = near_dedup(&recs);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].name, "b");
    }

    /// Records on different sequences must never suppress each other, however
    /// close their coordinates happen to be.
    #[test]
    fn test_near_dedup_resets_across_sequences() {
        let recs = vec![
            make_rec("chr1", 100, 270, "a", 80.0),
            make_rec("chr2", 100, 270, "b", 50.0), // identical coords, lower score
        ];
        let result = near_dedup(&recs);
        assert_eq!(result.len(), 2);
        assert_eq!(result[1].chrom, "chr2");
    }
}
