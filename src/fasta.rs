use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::error::PipelineError;

/// One sequence of the input, written out as its own FASTA file.
#[derive(Debug)]
pub struct SplitSeq {
    /// Name as nhmmer will report it: the header up to the first whitespace.
    pub name: String,
    /// Single-sequence FASTA in the temporary directory.
    pub path: PathBuf,
}

#[derive(Debug)]
pub struct Split {
    pub seqs: Vec<SplitSeq>,
    /// Total residues across the whole input, used to pin nhmmer's `-Z`.
    pub residues: u64,
}

/// Stream `input` and write every sequence to its own FASTA under `dir`.
///
/// The input is never held in memory: it is copied line by line into the
/// current output file. Headers are written through unchanged, so nhmmer
/// reports the same target names it would for the undivided file.
pub fn split_fasta(input: &Path, dir: &Path) -> Result<Split, PipelineError> {
    let reader = BufReader::new(File::open(input)?);

    let mut seqs: Vec<SplitSeq> = Vec::new();
    let mut residues: u64 = 0;
    let mut writer: Option<BufWriter<File>> = None;

    for line in reader.lines() {
        let line = line?;

        if let Some(header) = line.strip_prefix('>') {
            if let Some(mut w) = writer.take() {
                w.flush()?;
            }
            let name = header.split_whitespace().next().unwrap_or("").to_string();
            let path = dir.join(format!("seq_{:06}.fa", seqs.len()));
            let mut w = BufWriter::new(File::create(&path)?);
            writeln!(w, "{}", line)?;
            writer = Some(w);
            seqs.push(SplitSeq { name, path });
        } else if let Some(w) = writer.as_mut() {
            residues += line.trim_end().len() as u64;
            writeln!(w, "{}", line)?;
        }
        // Anything before the first header is not part of a sequence.
    }

    if let Some(mut w) = writer.take() {
        w.flush()?;
    }

    if seqs.is_empty() {
        return Err(PipelineError::NoSequences(input.to_path_buf()));
    }

    Ok(Split { seqs, residues })
}

/// Pull the bases over each of `intervals` (0-based, half-open, sorted by
/// start) out of a single-sequence FASTA.
///
/// The file is streamed once and only the requested bases are retained, so a
/// handful of junction gaps costs a pass over the file rather than a copy of
/// the chromosome in memory.
pub fn read_intervals(
    path: &Path,
    intervals: &[(u64, u64)],
) -> Result<Vec<Vec<u8>>, PipelineError> {
    let mut out: Vec<Vec<u8>> = intervals
        .iter()
        .map(|(s, e)| Vec::with_capacity(e.saturating_sub(*s) as usize))
        .collect();
    if intervals.is_empty() {
        return Ok(out);
    }

    let reader = BufReader::new(File::open(path)?);
    let mut offset: u64 = 0;
    let mut first = 0usize;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            continue;
        }
        let bases = line.trim_end().as_bytes();
        let line_end = offset + bases.len() as u64;

        // Intervals that ended before this line will never be seen again.
        while first < intervals.len() && intervals[first].1 <= offset {
            first += 1;
        }
        let mut i = first;
        while i < intervals.len() && intervals[i].0 < line_end {
            let from = intervals[i].0.max(offset);
            let to = intervals[i].1.min(line_end);
            if from < to {
                out[i].extend_from_slice(&bases[(from - offset) as usize..(to - offset) as usize]);
            }
            i += 1;
        }

        offset = line_end;
    }

    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    /// Write `content` to a scratch dir under the crate's target directory and
    /// return both paths. Kept out of the system temp dir on purpose.
    fn scratch(name: &str, content: &str) -> (PathBuf, PathBuf) {
        let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("target")
            .join("test-scratch")
            .join(name);
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        let input = dir.join("in.fa");
        fs::write(&input, content).unwrap();
        (input, dir)
    }

    #[test]
    fn test_split_writes_one_file_per_sequence() {
        let (input, dir) = scratch("split_basic", ">chr1 first\nACGT\nACGT\n>chr2\nTTTT\n");
        let split = split_fasta(&input, &dir).unwrap();

        assert_eq!(split.seqs.len(), 2);
        assert_eq!(split.seqs[0].name, "chr1");
        assert_eq!(split.seqs[1].name, "chr2");
        assert_eq!(split.residues, 12);

        // The header must survive verbatim, description and all.
        let first = fs::read_to_string(&split.seqs[0].path).unwrap();
        assert_eq!(first, ">chr1 first\nACGT\nACGT\n");
        let second = fs::read_to_string(&split.seqs[1].path).unwrap();
        assert_eq!(second, ">chr2\nTTTT\n");
    }

    #[test]
    fn test_read_intervals_spans_line_breaks() {
        let (input, _) = scratch("read_intervals", ">chr1\nACGTA\nCGTAC\nGTACG\n");
        // 0-based half-open over the concatenated residues "ACGTACGTACGTACG"
        let got = read_intervals(&input, &[(0, 3), (4, 7), (13, 15)]).unwrap();
        assert_eq!(got[0], b"ACG");
        assert_eq!(got[1], b"ACG");
        assert_eq!(got[2], b"CG");
    }

    #[test]
    fn test_read_intervals_empty_request() {
        let (input, _) = scratch("read_intervals_none", ">chr1\nACGT\n");
        assert!(read_intervals(&input, &[]).unwrap().is_empty());
    }

    #[test]
    fn test_split_rejects_empty_input() {
        let (input, dir) = scratch("split_empty", "\n\n");
        assert!(split_fasta(&input, &dir).is_err());
    }
}
