use crate::bed::BedRecord;
use std::collections::HashSet;

/// Stage 1: bedmap --max-element --fraction-either 0.1
///
/// Groups overlapping intervals (≥10% reciprocal overlap) and keeps the
/// highest-scoring element from each group.
/// Input MUST be sorted by chrom + start.
pub fn bedmap_max_element(records: &[BedRecord]) -> Vec<BedRecord> {
    if records.is_empty() {
        return vec![];
    }

    let mut result: Vec<BedRecord> = Vec::new();

    // Sweep-line: for each record, find all overlapping records in the current
    // cluster and keep the one with the highest score.
    // A cluster is a set of records that have ≥10% reciprocal overlap with
    // at least one other record in the cluster.

    // Simple approach: for each record, check if it overlaps with the current
    // best in the cluster. If yes, keep the higher-scoring one. If no overlap,
    // emit the current best and start a new cluster.

    let mut cluster: Vec<BedRecord> = vec![records[0].clone()];

    for rec in records.iter().skip(1) {
        // Check if this record overlaps with ANY record in the current cluster
        let overlaps = cluster
            .iter()
            .any(|c| c.chrom == rec.chrom && has_reciprocal_overlap(c, rec, 0.1));

        if overlaps {
            cluster.push(rec.clone());
        } else {
            // Emit the max-score element from the cluster
            if let Some(best) = cluster.iter().max_by(|a, b| {
                a.score
                    .partial_cmp(&b.score)
                    .unwrap_or(std::cmp::Ordering::Equal)
            }) {
                result.push(best.clone());
            }
            cluster.clear();
            cluster.push(rec.clone());
        }
    }

    // Emit last cluster
    if let Some(best) = cluster.iter().max_by(|a, b| {
        a.score
            .partial_cmp(&b.score)
            .unwrap_or(std::cmp::Ordering::Equal)
    }) {
        result.push(best.clone());
    }

    result
}

/// Check if two intervals have ≥ frac reciprocal overlap.
/// "fraction-either" means: overlap / min(len_a, len_b) >= frac
/// Actually bedmap --fraction-either means: overlap >= frac * len_a AND overlap >= frac * len_b
/// i.e., overlap must be at least frac of BOTH intervals.
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
    // fraction-either: overlap >= frac * len_a OR overlap >= frac * len_b
    // bedops bedmap --fraction-either: the overlap must be >= frac of EITHER element
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
/// Sequential pass: if start or end is within ±10bp of previous record,
/// keep the one with higher score (tie-break by length).
pub fn near_dedup(records: &[BedRecord]) -> Vec<BedRecord> {
    if records.is_empty() {
        return vec![];
    }

    let mut result: Vec<BedRecord> = Vec::new();
    let mut prev_start: i64 = -100;
    let mut prev_end: i64 = -100;
    let mut prev_score: f64 = 0.0;
    let mut prev_length: i64 = 0;

    for rec in records {
        let start = rec.start as i64;
        let end = rec.end as i64;
        let length = end - start;
        let score = rec.score;

        let start_near = (start - prev_start).abs() < 10;
        let end_near = (end - prev_end).abs() < 10;

        if start_near || end_near {
            // Scores nearly equal (within 0.01)?
            if (prev_score - score).abs() < 0.01 {
                // Same score → keep longer
                if length > prev_length {
                    result.pop();
                    // fall through to add current
                } else {
                    continue;
                }
            } else if score > prev_score {
                // Current is better → replace previous
                if !result.is_empty() {
                    result.pop();
                }
                // fall through to add current
            } else {
                // Previous is better → skip current
                continue;
            }
        }

        result.push(rec.clone());
        prev_start = start;
        prev_end = end;
        prev_score = score;
        prev_length = length;
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
    use crate::bed::Rgb;

    fn make_rec(chrom: &str, start: u64, end: u64, name: &str, score: f64) -> BedRecord {
        BedRecord {
            chrom: chrom.into(),
            start,
            end,
            name: name.into(),
            score,
            strand: '+',
            thick_start: start,
            thick_end: end,
            color: Rgb(0, 0, 0),
        }
    }

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
            make_rec("chr1", 105, 272, "b", 80.0), // within ±10bp
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
}
