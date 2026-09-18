use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::error::PipelineError;

/// One sequence of the input, written out as its own FASTA file.
#[allow(dead_code)]
#[derive(Debug)]
pub struct SplitSeq {
    /// Name as nhmmer will report it: the header up to the first whitespace.
    pub name: String,
    /// Single-sequence FASTA in the temporary directory.
    pub path: PathBuf,
    /// Residues in this sequence, which is what the work of scanning it costs.
    pub len: u64,
}

#[allow(dead_code)]
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
#[allow(dead_code)]
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
            seqs.push(SplitSeq { name, path, len: 0 });
        } else if let Some(w) = writer.as_mut() {
            let n = line.trim_end().len() as u64;
            residues += n;
            if let Some(last) = seqs.last_mut() {
                last.len += n;
            }
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
#[allow(dead_code)]
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

/// One FASTA sequence stored completely in memory.
#[derive(Debug, Clone)]
pub struct MemoryFastaSeq {
    pub name: String,
    pub seq: Vec<u8>,
}

/// Read all sequences from a FASTA file directly into memory.
pub fn read_fasta_in_memory(path: &Path) -> Result<Vec<MemoryFastaSeq>, PipelineError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut seqs = Vec::new();
    let mut name = String::new();
    let mut seq = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            if !name.is_empty() {
                seqs.push(MemoryFastaSeq {
                    name: std::mem::take(&mut name),
                    seq: std::mem::take(&mut seq),
                });
            }
            let first = header.split_whitespace().next().unwrap_or("").to_string();
            name = first;
        } else if !name.is_empty() {
            seq.extend(
                line.trim_end()
                    .as_bytes()
                    .iter()
                    .filter(|b| b.is_ascii_alphabetic())
                    .map(|b| b.to_ascii_uppercase()),
            );
        }
    }

    if !name.is_empty() {
        seqs.push(MemoryFastaSeq { name, seq });
    }

    if seqs.is_empty() {
        return Err(PipelineError::NoSequences(path.to_path_buf()));
    }

    Ok(seqs)
}

/// Extract intervals from an in-memory sequence without any disk reading.
pub fn extract_intervals_in_memory(seq: &[u8], intervals: &[(u64, u64)]) -> Vec<Vec<u8>> {
    let mut out = Vec::with_capacity(intervals.len());
    for &(s, e) in intervals {
        let start = (s as usize).min(seq.len());
        let end = (e as usize).min(seq.len());
        if start < end {
            out.push(seq[start..end].to_vec());
        } else {
            out.push(Vec::new());
        }
    }
    out
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
        assert_eq!((split.seqs[0].len, split.seqs[1].len), (8, 4));

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

    #[test]
    fn test_read_fasta_in_memory() {
        let (input, _) = scratch("mem_fasta", ">seq1 header info\nACGT\nACGT\n>seq2\nTTTT\n");
        let seqs = read_fasta_in_memory(&input).unwrap();
        assert_eq!(seqs.len(), 2);
        assert_eq!(seqs[0].name, "seq1");
        assert_eq!(seqs[0].seq, b"ACGTACGT");
        assert_eq!(seqs[1].name, "seq2");
        assert_eq!(seqs[1].seq, b"TTTT");

        let intervals = extract_intervals_in_memory(&seqs[0].seq, &[(0, 3), (4, 8)]);
        assert_eq!(intervals[0], b"ACG");
        assert_eq!(intervals[1], b"ACGT");
    }
}
