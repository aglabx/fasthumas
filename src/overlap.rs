use crate::bed::BedRecord;
use std::collections::HashSet;

/// Stage 1: bedmap --max-element --fraction-either 0.1
///
/// Groups overlapping intervals (>=10% overlap of either element) and keeps the
/// highest-scoring element from each group.
/// Input MUST be sorted by chrom + start.
pub fn bedmap_max_element(records: &[BedRecord]) -> Vec<BedRecord> {
    if records.is_empty() {
        return vec![];
    }

    let mut result: Vec<BedRecord> = Vec::new();

    // Sweep-line: for each record, check whether it overlaps any record in the
    // current cluster. If it does, join the cluster; otherwise emit the
    // cluster's max-scoring element and start a new one.
    let mut cluster: Vec<BedRecord> = vec![records[0].clone()];

    for rec in records.iter().skip(1) {
        let overlaps = cluster
            .iter()
            .any(|c| c.chrom == rec.chrom && has_reciprocal_overlap(c, rec, 0.1));

        if overlaps {
            cluster.push(rec.clone());
        } else {
            if let Some(best) = max_by_score(&cluster) {
                result.push(best.clone());
            }
            cluster.clear();
            cluster.push(rec.clone());
        }
    }

    if let Some(best) = max_by_score(&cluster) {
        result.push(best.clone());
    }

    result
}

fn max_by_score(records: &[BedRecord]) -> Option<&BedRecord> {
    records.iter().max_by(|a, b| {
        a.score
            .partial_cmp(&b.score)
            .unwrap_or(std::cmp::Ordering::Equal)
    })
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
        let result = bedmap_max_element(&recs);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].name, "b");
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
