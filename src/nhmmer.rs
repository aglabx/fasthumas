use std::io::{BufRead, BufReader};
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
/// Returns filtered, sorted BED records.
pub fn run_nhmmer(
    hmm_path: &Path,
    fasta_path: &Path,
    cpu_threads: usize,
    dbsize_mb: f64,
    score_threshold: f64,
) -> Result<Vec<BedRecord>, PipelineError> {
    let mut child = Command::new("nhmmer")
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
                PipelineError::NhmmerFailed {
                    path: fasta_path.to_path_buf(),
                    message: e.to_string(),
                }
            }
        })?;

    let stdout = child.stdout.take().unwrap();
    let reader = BufReader::new(stdout);

    let mut records: Vec<BedRecord> = Vec::new();

    for line in reader.lines() {
        let line = line.map_err(|e| PipelineError::NhmmerFailed {
            path: fasta_path.to_path_buf(),
            message: e.to_string(),
        })?;

        if let Some(hit) = parse_tblout_line(&line) {
            if let Some(bed) = hit_to_bed(&hit, score_threshold) {
                records.push(bed);
            }
        }
    }

    let status = child.wait().map_err(|e| PipelineError::NhmmerFailed {
        path: fasta_path.to_path_buf(),
        message: e.to_string(),
    })?;

    if !status.success() {
        return Err(PipelineError::NhmmerFailed {
            path: fasta_path.to_path_buf(),
            message: format!("nhmmer exited with status {}", status),
        });
    }

    // Sort by sequence name + start, matching `sort -k 1.4,1 -k 2,2n`
    sort_bed_records(&mut records);

    // 3-stage overlap filtering
    let filtered = filter_overlaps(&records);

    Ok(filtered)
}
