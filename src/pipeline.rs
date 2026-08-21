use std::fs;
use std::path::Path;

use log::{error, info};
use rayon::prelude::*;

use crate::bed::{BedRecord, Rgb};
use crate::config::Config;
use crate::error::PipelineError;
use crate::nhmmer::run_nhmmer;

/// Run the full pipeline on all FASTA files.
pub fn run(config: &Config) -> Result<(), PipelineError> {
    fs::create_dir_all(&config.output_dir)?;

    let fasta_files = config.discover_fasta_files()?;
    info!(
        "Found {} FASTA files; {} CPU budget split into {} parallel jobs \u{00d7} {} nhmmer threads",
        fasta_files.len(),
        config.total_threads,
        config.parallel_jobs,
        config.nhmmer_threads,
    );

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(config.parallel_jobs)
        .build()
        .map_err(|e| PipelineError::Io(std::io::Error::other(e.to_string())))?;

    let errors: Vec<String> = pool.install(|| {
        fasta_files
            .par_iter()
            .filter_map(|fasta| match process_one_file(config, fasta) {
                Ok(()) => None,
                Err(e) => {
                    let msg = format!("{}: {}", fasta.display(), e);
                    error!("{}", msg);
                    Some(msg)
                }
            })
            .collect()
    });

    if errors.is_empty() {
        info!("Pipeline complete: {} files processed", fasta_files.len());
        Ok(())
    } else {
        Err(PipelineError::NhmmerFailed {
            path: config.input_dir.clone(),
            message: format!("{} files failed:\n{}", errors.len(), errors.join("\n")),
        })
    }
}

/// Process a single FASTA file: HOR track + SF track.
fn process_one_file(config: &Config, fasta: &Path) -> Result<(), PipelineError> {
    let stem = fasta.file_stem().unwrap_or_default().to_string_lossy();
    info!("Processing {}", stem);

    // --- HOR track ---
    let hor_records = run_nhmmer(
        &config.hor_hmm,
        fasta,
        config.nhmmer_threads,
        config.score_threshold,
    )?;

    // AS-HOR+SF: all filtered records
    let hor_sf_path = config.output_path(fasta, "AS-HOR+SF");
    write_bed(&hor_sf_path, &hor_records)?;
    info!(
        "  HOR+SF: {} records → {}",
        hor_records.len(),
        hor_sf_path.display()
    );

    // AS-HOR: skip SF monomers (2-char model names)
    let hor_only: Vec<BedRecord> = hor_records
        .into_iter()
        .filter(|r| r.name.split('.').next().map_or(true, |base| base.len() > 2))
        .collect();

    let hor_path = config.output_path(fasta, "AS-HOR");
    write_bed(&hor_path, &hor_only)?;
    info!(
        "  HOR:    {} records → {}",
        hor_only.len(),
        hor_path.display()
    );

    // --- SF track ---
    let sf_records = run_nhmmer(
        &config.sf_hmm,
        fasta,
        config.nhmmer_threads,
        config.score_threshold,
    )?;

    // AS-SF: all SF records
    let sf_path = config.output_path(fasta, "AS-SF");
    write_bed(&sf_path, &sf_records)?;
    info!(
        "  SF:     {} records → {}",
        sf_records.len(),
        sf_path.display()
    );

    // AS-strand: recolor by strand (+→blue, -→red)
    let strand_records: Vec<BedRecord> = sf_records
        .into_iter()
        .map(|mut r| {
            r.color = match r.strand {
                '+' => Rgb(0, 0, 255),
                '-' => Rgb(255, 0, 0),
                _ => r.color,
            };
            r
        })
        .collect();

    let strand_path = config.output_path(fasta, "AS-strand");
    write_bed(&strand_path, &strand_records)?;
    info!(
        "  strand: {} records → {}",
        strand_records.len(),
        strand_path.display()
    );

    Ok(())
}

/// Write BED records to a file.
fn write_bed(path: &Path, records: &[BedRecord]) -> Result<(), PipelineError> {
    let mut content = String::with_capacity(records.len() * 80);
    for rec in records {
        content.push_str(&rec.to_bed9_string());
        content.push('\n');
    }
    fs::write(path, content)?;
    Ok(())
}
