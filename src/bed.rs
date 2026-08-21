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

/// Sort BED records by chrom (from 4th char onward to match `sort -k 1.4,1`),
/// then by start coordinate.
pub fn sort_bed_records(records: &mut [BedRecord]) {
    records.sort_by(|a, b| {
        let a_key = if a.chrom.len() >= 3 {
            &a.chrom[3..]
        } else {
            &a.chrom[..]
        };
        let b_key = if b.chrom.len() >= 3 {
            &b.chrom[3..]
        } else {
            &b.chrom[..]
        };
        a_key.cmp(b_key).then(a.start.cmp(&b.start))
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bed9_format() {
        let rec = BedRecord {
            chrom: "chr1".to_string(),
            start: 100,
            end: 200,
            name: "S1C3H1L.1".to_string(),
            score: 123.4,
            strand: '+',
            thick_start: 100,
            thick_end: 200,
            color: Rgb(171, 171, 7),
        };
        assert_eq!(
            rec.to_bed9_string(),
            "chr1\t100\t200\tS1C3H1L.1\t123.4\t+\t100\t200\t171,171,7"
        );
    }

    #[test]
    fn test_sort() {
        let mut recs = vec![
            BedRecord {
                chrom: "chr2".into(),
                start: 50,
                end: 100,
                name: "a".into(),
                score: 1.0,
                strand: '+',
                thick_start: 50,
                thick_end: 100,
                color: Rgb(0, 0, 0),
            },
            BedRecord {
                chrom: "chr1".into(),
                start: 200,
                end: 300,
                name: "b".into(),
                score: 2.0,
                strand: '+',
                thick_start: 200,
                thick_end: 300,
                color: Rgb(0, 0, 0),
            },
            BedRecord {
                chrom: "chr1".into(),
                start: 100,
                end: 200,
                name: "c".into(),
                score: 3.0,
                strand: '+',
                thick_start: 100,
                thick_end: 200,
                color: Rgb(0, 0, 0),
            },
        ];
        sort_bed_records(&mut recs);
        assert_eq!(recs[0].chrom, "chr1");
        assert_eq!(recs[0].start, 100);
        assert_eq!(recs[1].chrom, "chr1");
        assert_eq!(recs[1].start, 200);
        assert_eq!(recs[2].chrom, "chr2");
    }
}
