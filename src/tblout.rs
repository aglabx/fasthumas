use crate::bed::BedRecord;
use crate::color::get_color;

/// A parsed hit from nhmmer --tblout output.
#[derive(Debug)]
pub struct HmmHit {
    pub target_name: String,
    pub query_name: String,
    pub hmm_from: u64,
    pub hmm_to: u64,
    pub ali_from: u64,
    pub ali_to: u64,
    pub strand: char,
    pub score: f64,
}

/// Parse a single tblout line into an HmmHit.
/// Returns None for comment lines or unparseable lines.
#[allow(dead_code)]
pub fn parse_tblout_line(line: &str) -> Option<HmmHit> {
    if line.starts_with('#') || line.trim().is_empty() {
        return None;
    }

    let fields: Vec<&str> = line.split_whitespace().collect();
    if fields.len() < 16 {
        return None;
    }

    let target_name = fields[0].to_string();
    let query_name = fields[2].to_string();
    // nhmmer tblout: fields are 1-based
    let hmm_from: u64 = fields[4].parse().ok()?;
    let hmm_to: u64 = fields[5].parse().ok()?;
    let ali_from: u64 = fields[6].parse().ok()?;
    let ali_to: u64 = fields[7].parse().ok()?;
    // fields[8]/fields[9] are envelope coordinates; this port always uses ali
    // coordinates, matching how hmmer-run.sh invokes hmmertblout2bed.awk.
    let strand = fields[11].chars().next()?;
    let score: f64 = fields[13].parse().ok()?;

    Some(HmmHit {
        target_name,
        query_name,
        hmm_from,
        hmm_to,
        ali_from,
        ali_to,
        strand,
        score,
    })
}

/// Convert an HmmHit to a BedRecord, applying the score/length threshold.
///
/// nhmmer reports 1-based inclusive alignment coordinates in *query* order, so
/// on the minus strand `ali_from > ali_to` — column 7 is the end of the hit on
/// the target, not its start. The hit spans `[min, max]` 1-based inclusive on
/// the target either way, which in 0-based half-open BED is `[min - 1, max)`.
///
/// The original `hmmertblout2bed.awk` computed `aliFrom = $7 - 1; aliTo = $8`
/// and then printed `aliTo, aliFrom` for the minus strand — putting the -1 on
/// the end instead of the start. That shifts every minus-strand interval one
/// base to the right and makes it two bases too short, and the same wrong
/// length then feeds the score/length threshold. This port does not reproduce
/// that; minus-strand output therefore differs from the original pipeline.
#[allow(dead_code)]
pub fn hit_to_bed(hit: &HmmHit, threshold: f64) -> Option<BedRecord> {
    let (lo, hi) = if hit.ali_from <= hit.ali_to {
        (hit.ali_from, hit.ali_to)
    } else {
        (hit.ali_to, hit.ali_from)
    };

    // 1-based inclusive -> 0-based half-open
    let start = lo.checked_sub(1)?;
    let end = hi;

    let ali_length = end.checked_sub(start)?;
    if ali_length == 0 {
        return None;
    }

    let score_per_len = hit.score / ali_length as f64;
    // Match AWK: sprintf("%.1f", score/aliLength) < th
    let formatted: f64 = format!("{:.1}", score_per_len)
        .parse()
        .unwrap_or(score_per_len);
    if formatted < threshold {
        return None;
    }

    let color = get_color(&hit.query_name);

    Some(BedRecord {
        chrom: hit.target_name.clone(),
        start,
        end,
        name: hit.query_name.clone(),
        score: hit.score,
        strand: hit.strand,
        thick_start: start,
        thick_end: end,
        color,
        hmm_from: hit.hmm_from,
        hmm_to: hit.hmm_to,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hit(ali_from: u64, ali_to: u64, strand: char, score: f64) -> HmmHit {
        HmmHit {
            target_name: "chr1".into(),
            query_name: "S1C3H1L.1".into(),
            hmm_from: 1,
            hmm_to: 171,
            ali_from,
            ali_to,
            strand,
            score,
        }
    }

    #[test]
    fn test_parse_comment() {
        assert!(parse_tblout_line("# this is a comment").is_none());
        assert!(parse_tblout_line("").is_none());
    }

    #[test]
    fn test_parse_valid_line() {
        // Fake tblout line with 16+ fields
        let line = "chr1 - S1C3H1L.1 - 3 168 1000 1171 998 1173 - + - 120.5 0.0001 -";
        let h = parse_tblout_line(line).unwrap();
        assert_eq!(h.target_name, "chr1");
        assert_eq!(h.query_name, "S1C3H1L.1");
        assert_eq!((h.hmm_from, h.hmm_to), (3, 168));
        assert_eq!((h.ali_from, h.ali_to), (1000, 1171));
        assert_eq!(h.strand, '+');
        assert!((h.score - 120.5).abs() < 0.01);
    }

    #[test]
    fn test_hit_to_bed_plus_strand() {
        let bed = hit_to_bed(&hit(1000, 1171, '+', 120.5), 0.7).unwrap();
        assert_eq!(bed.start, 999); // 1000-1 = 999
        assert_eq!(bed.end, 1171);
        assert_eq!(bed.length(), 172);
        assert_eq!(bed.strand, '+');
    }

    /// On the minus strand nhmmer swaps the columns: $7 is the end, $8 the
    /// start. The interval must come out the same length as the equivalent
    /// plus-strand hit, not two bases shorter (which is what the original AWK
    /// produced).
    #[test]
    fn test_hit_to_bed_minus_strand() {
        let bed = hit_to_bed(&hit(2125777, 2125613, '-', 123.1), 0.7).unwrap();
        assert_eq!(bed.start, 2125612);
        assert_eq!(bed.end, 2125777);
        assert_eq!(bed.length(), 165);
        assert_eq!(bed.strand, '-');
    }

    /// A plus hit and the same span reported on the minus strand must yield
    /// identical coordinates.
    #[test]
    fn test_strands_agree_on_span() {
        let p = hit_to_bed(&hit(500, 670, '+', 130.0), 0.7).unwrap();
        let m = hit_to_bed(&hit(670, 500, '-', 130.0), 0.7).unwrap();
        assert_eq!((p.start, p.end), (m.start, m.end));
    }

    #[test]
    fn test_threshold_filter() {
        // 10/172 = 0.058, way below 0.7
        assert!(hit_to_bed(&hit(100, 271, '+', 10.0), 0.7).is_none());
    }

    #[test]
    fn test_model_coordinates_are_kept() {
        let mut h = hit(1000, 1171, '+', 120.5);
        h.hmm_from = 3;
        h.hmm_to = 168;
        let bed = hit_to_bed(&h, 0.7).unwrap();
        assert_eq!((bed.hmm_from, bed.hmm_to), (3, 168));
    }
}
