use std::fmt;

#[derive(Debug, Clone, PartialEq)]
pub struct Rgb(pub u8, pub u8, pub u8);

impl fmt::Display for Rgb {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{},{},{}", self.0, self.1, self.2)
    }
}

#[derive(Debug, Clone)]
pub struct BedRecord {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub name: String,
    pub score: f64,
    pub strand: char,
    pub thick_start: u64,
    pub thick_end: u64,
    pub color: Rgb,
    /// 1-based model coordinates of the alignment, straight from nhmmer's
    /// tblout. Not written to the BED; the junction pass uses them to work out
    /// how many model positions were trimmed off each end.
    pub hmm_from: u64,
    pub hmm_to: u64,
}

impl BedRecord {
    pub fn to_bed9_string(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{:.1}\t{}\t{}\t{}\t{}",
            self.chrom,
            self.start,
            self.end,
            self.name,
            self.score,
            self.strand,
            self.thick_start,
            self.thick_end,
            self.color,
        )
    }

    pub fn length(&self) -> u64 {
        self.end.saturating_sub(self.start)
    }
}

/// Sort BED records by sequence name, then by start coordinate.
///
/// The original pipeline used `sort -k 1.4,1 -k 2,2n`, i.e. it compared names
/// from the 4th character on — a shortcut that assumes a `chr` prefix and that
/// only ever ran on single-chromosome files, where the name key is a no-op.
/// This port takes multi-sequence FASTA, so it compares the whole name; for
/// uniformly prefixed names the resulting order is the same, and unlike the
/// original it cannot interleave records from different sequences (which the
/// order-dependent filtering stages downstream rely on).
pub fn sort_bed_records(records: &mut [BedRecord]) {
    records.sort_by(|a, b| {
        a.chrom
            .cmp(&b.chrom)
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
    });
}

#[cfg(test)]
pub(crate) fn test_record(chrom: &str, start: u64, end: u64, name: &str, score: f64) -> BedRecord {
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
        hmm_from: 1,
        hmm_to: end - start,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bed9_format() {
        let r = BedRecord {
            chrom: "chr1".to_string(),
            start: 100,
            end: 200,
            name: "S1C3H1L.1".to_string(),
            score: 123.4,
            strand: '+',
            thick_start: 100,
            thick_end: 200,
            color: Rgb(171, 171, 7),
            hmm_from: 1,
            hmm_to: 100,
        };
        assert_eq!(
            r.to_bed9_string(),
            "chr1\t100\t200\tS1C3H1L.1\t123.4\t+\t100\t200\t171,171,7"
        );
    }

    #[test]
    fn test_sort() {
        let mut recs = vec![
            test_record("chr2", 50, 100, "a", 1.0),
            test_record("chr1", 200, 300, "b", 1.0),
            test_record("chr1", 100, 200, "c", 1.0),
        ];
        sort_bed_records(&mut recs);
        assert_eq!(recs[0].chrom, "chr1");
        assert_eq!(recs[0].start, 100);
        assert_eq!(recs[1].chrom, "chr1");
        assert_eq!(recs[1].start, 200);
        assert_eq!(recs[2].chrom, "chr2");
    }

    /// Sequence names of different shapes must still group cleanly, which the
    /// old `chrom[3..]` key could not guarantee.
    #[test]
    fn test_sort_groups_mixed_names() {
        let mut recs = vec![
            test_record("scaffold_2", 10, 20, "a", 1.0),
            test_record("chr1", 30, 40, "b", 1.0),
            test_record("scaffold_2", 5, 8, "c", 1.0),
            test_record("chr1", 10, 20, "d", 1.0),
        ];
        sort_bed_records(&mut recs);
        let names: Vec<&str> = recs.iter().map(|r| r.chrom.as_str()).collect();
        assert_eq!(names, vec!["chr1", "chr1", "scaffold_2", "scaffold_2"]);
        assert_eq!(recs[0].start, 10);
        assert_eq!(recs[2].start, 5);
    }
}
