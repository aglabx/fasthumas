use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use std::process::{Command, Stdio};

use crate::bed::{sort_bed_records, BedRecord};
use crate::error::PipelineError;
use crate::overlap::filter_overlaps;
use crate::tblout::{hit_to_bed, parse_tblout_line};

/// Run nhmmer with --tblout piped to stdout, parse hits on the fly,
/// apply the score threshold, sort, and filter overlaps.
///
/// `dbsize_mb` is passed through as `-Z`. The pipeline searches one sequence
/// per nhmmer call but pins the database size to the whole input, so E-values
/// — and therefore which hits clear nhmmer's own reporting threshold — do not
/// depend on how the input happened to be divided.
///
/// `--dna` settles the alphabet instead of letting nhmmer guess it. The guess
/// is made from the first few thousand residues of the target and a
/// low-complexity start can defeat it: the human telomere, `TAACCCTAACCCT...`,
/// uses only T, A and C, each equally valid as an amino acid, and nhmmer then
/// refuses the whole file. In CHM13 that hits chr22 and chrY, whose leading
/// telomere runs past 4 kb. The flag only reaches the target from HMMER 3.4
/// on; under 3.3.2 it is accepted and ignored, and those two chromosomes still
/// fail.
///
/// Returns filtered, sorted BED records.
#[allow(dead_code)]
pub fn run_nhmmer(
    hmm_path: &Path,
    fasta_path: &Path,
    cpu_threads: usize,
    dbsize_mb: f64,
    score_threshold: f64,
) -> Result<Vec<BedRecord>, PipelineError> {
    let failed = |message: String| PipelineError::NhmmerFailed {
        path: fasta_path.to_path_buf(),
        message,
    };

    let mut child = Command::new("nhmmer")
        .arg("--dna")
        .arg("--cpu")
        .arg(cpu_threads.to_string())
        .arg("-Z")
        .arg(format!("{}", dbsize_mb))
        .arg("--notextw")
        .arg("--noali")
        .arg("--tblout")
        .arg("/dev/stdout")
        .arg("-o")
        .arg("/dev/null")
        .arg(hmm_path)
        .arg(fasta_path)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| {
            if e.kind() == std::io::ErrorKind::NotFound {
                PipelineError::NhmmerNotFound
            } else {
                failed(e.to_string())
            }
        })?;

    // Drain stderr on its own thread. nhmmer says why it gave up there, and a
    // pipe nobody reads will eventually block the process that writes to it.
    let stderr = child.stderr.take().expect("stderr was piped");
    let stderr_reader = std::thread::spawn(move || {
        let mut text = String::new();
        let _ = BufReader::new(stderr).read_to_string(&mut text);
        text
    });

    let stdout = child.stdout.take().expect("stdout was piped");
    let reader = BufReader::new(stdout);

    let mut records: Vec<BedRecord> = Vec::new();
    let mut read_error = None;

    for line in reader.lines() {
        let line = match line {
            Ok(line) => line,
            Err(e) => {
                read_error = Some(e.to_string());
                break;
            }
        };

        if let Some(hit) = parse_tblout_line(&line) {
            if let Some(bed) = hit_to_bed(&hit, score_threshold) {
                records.push(bed);
            }
        }
    }

    let status = child.wait().map_err(|e| failed(e.to_string()))?;
    let complaint = stderr_reader.join().unwrap_or_default();

    if !status.success() {
        return Err(failed(match last_line(&complaint) {
            Some(line) => format!("nhmmer exited with status {}: {}", status, line),
            None => format!("nhmmer exited with status {}", status),
        }));
    }
    if let Some(e) = read_error {
        return Err(failed(e));
    }

    // Sort by sequence name + start, matching `sort -k 1.4,1 -k 2,2n`
    sort_bed_records(&mut records);

    // 3-stage overlap filtering
    let filtered = filter_overlaps(&records);

    Ok(filtered)
}

/// The last thing nhmmer said before giving up, which is the part that explains
/// why.
#[allow(dead_code)]
fn last_line(text: &str) -> Option<String> {
    text.lines()
        .map(str::trim)
        .rfind(|l| !l.is_empty())
        .map(str::to_string)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_last_line_picks_the_complaint() {
        let text = "\n# nhmmer :: search a DNA model\n\nError: Invalid alphabet type in target for nhmmer. Expect DNA or RNA.\n\n";
        assert_eq!(
            last_line(text).unwrap(),
            "Error: Invalid alphabet type in target for nhmmer. Expect DNA or RNA."
        );
        assert!(last_line("   \n\n").is_none());
    }
}
