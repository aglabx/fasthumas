use crate::bed::BedRecord;
use crate::color::get_color;

/// A parsed hit from nhmmer --tblout output.
#[derive(Debug)]
pub struct HmmHit {
    pub target_name: String,
    pub query_name: String,
    pub ali_from: u64,
    pub ali_to: u64,
    pub strand: char,
    pub score: f64,
}

/// Parse a single tblout line into an HmmHit.
/// Returns None for comment lines or unparseable lines.
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
    let ali_from: u64 = fields[6].parse().ok()?;
    let ali_to: u64 = fields[7].parse().ok()?;
    // fields[8]/fields[9] are envelope coordinates; this port always uses ali
    // coordinates, matching how hmmer-run.sh invokes hmmertblout2bed.awk.
    let strand = fields[11].chars().next()?;
    let score: f64 = fields[13].parse().ok()?;

    Some(HmmHit {
        target_name,
        query_name,
        ali_from,
        ali_to,
        strand,
        score,
    })
}

/// Convert an HmmHit to a BedRecord, applying score/length threshold filter.
/// Uses ali coordinates (not env), matching the default AWK behavior.
/// Returns None if the hit fails the threshold.
pub fn hit_to_bed(hit: &HmmHit, threshold: f64) -> Option<BedRecord> {
    // Convert 1-based ali_from to 0-based BED start
    // For + strand: start = ali_from-1, end = ali_to
    // For - strand: swap so start < end (ali_to-1, ali_from) — matching AWK behavior
    let (start, end) = if hit.strand == '+' {
        (hit.ali_from.saturating_sub(1), hit.ali_to)
    } else {
        // AWK prints: aliTo, aliFrom for minus strand
        // aliTo is the smaller coordinate for minus strand in nhmmer output
        // But AWK does aliTo, aliFrom which means it swaps
        // Actually in nhmmer, for - strand: ali_from > ali_to
        // AWK: print targetName, aliTo, aliFrom — so end coord first, then start
        // Wait, re-reading AWK: for minus: print targetName, aliTo, aliFrom
        // aliTo is the smaller number for - strand, aliFrom is bigger
        // So BED start = aliTo (smaller) -1 for 0-based, end = aliFrom
        // Actually the AWK doesn't do -1 on aliTo for minus strand
        // Let me re-read: aliFrom = $7-1, aliTo = $8 (no adjustment)
        // For -: print targetName, aliTo, aliFrom
        // So start=aliTo (which is $8, no -1), end=aliFrom (which is $7-1)
        // But $7 > $8 for minus strand, so aliFrom-1 > aliTo
        // print: aliTo, aliFrom → smaller, larger → correct BED
        (hit.ali_to, hit.ali_from.saturating_sub(1))
    };

    // Threshold: score / length >= th
    let ali_length = if end > start {
        end - start
    } else {
        return None;
    };
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
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_comment() {
        assert!(parse_tblout_line("# this is a comment").is_none());
        assert!(parse_tblout_line("").is_none());
    }

    #[test]
    fn test_parse_valid_line() {
        // Fake tblout line with 16+ fields
        let line = "chr1 - S1C3H1L.1 - - - 1000 1171 998 1173 - + - 120.5 0.0001 -";
        let hit = parse_tblout_line(line).unwrap();
        assert_eq!(hit.target_name, "chr1");
        assert_eq!(hit.query_name, "S1C3H1L.1");
        assert_eq!(hit.ali_from, 1000);
        assert_eq!(hit.ali_to, 1171);
        assert_eq!(hit.strand, '+');
        assert!((hit.score - 120.5).abs() < 0.01);
    }

    #[test]
    fn test_hit_to_bed_plus_strand() {
        let hit = HmmHit {
            target_name: "chr1".into(),
            query_name: "S1C3H1L.1".into(),
            ali_from: 1000,
            ali_to: 1171,
            strand: '+',
            score: 120.5,
        };
        let bed = hit_to_bed(&hit, 0.7).unwrap();
        assert_eq!(bed.start, 999); // 1000-1 = 999
        assert_eq!(bed.end, 1171);
        assert_eq!(bed.strand, '+');
    }

    #[test]
    fn test_threshold_filter() {
        let hit = HmmHit {
            target_name: "chr1".into(),
            query_name: "Aa".into(),
            ali_from: 100,
            ali_to: 271,
            strand: '+',
            score: 10.0, // 10/171 = 0.058, way below 0.7
        };
        assert!(hit_to_bed(&hit, 0.7).is_none());
    }
}
