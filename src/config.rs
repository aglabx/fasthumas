use std::path::{Path, PathBuf};

use crate::error::PipelineError;

#[derive(Debug, Clone)]
pub struct Config {
    pub input_dir: PathBuf,
    pub hor_hmm: PathBuf,
    pub sf_hmm: PathBuf,
    pub output_dir: PathBuf,
    pub total_threads: usize,
    pub parallel_jobs: usize,
    pub nhmmer_threads: usize,
    pub score_threshold: f64,
}

impl Config {
    pub fn new(
        input_dir: PathBuf,
        hor_hmm: PathBuf,
        sf_hmm: PathBuf,
        output_dir: PathBuf,
        total_threads: usize,
        score_threshold: f64,
    ) -> Result<Self, PipelineError> {
        if !input_dir.is_dir() {
            return Err(PipelineError::NoInputFiles(input_dir));
        }
        if !hor_hmm.exists() {
            return Err(PipelineError::HmmNotFound(hor_hmm));
        }
        if !sf_hmm.exists() {
            return Err(PipelineError::HmmNotFound(sf_hmm));
        }

        // Auto-compute: target ~4 nhmmer threads per job
        let nhmmer_threads = 4.min(total_threads);
        let parallel_jobs = (total_threads / nhmmer_threads).max(1);

        Ok(Config {
            input_dir,
            hor_hmm,
            sf_hmm,
            output_dir,
            total_threads,
            parallel_jobs,
            nhmmer_threads,
            score_threshold,
        })
    }

    /// Discover all .fa files in input_dir.
    pub fn discover_fasta_files(&self) -> Result<Vec<PathBuf>, PipelineError> {
        let pattern = self.input_dir.join("*.fa");
        let pattern_str = pattern.to_string_lossy();

        let mut files: Vec<PathBuf> = glob::glob(&pattern_str)
            .map_err(|e| PipelineError::Io(std::io::Error::new(
                std::io::ErrorKind::Other,
                e.to_string(),
            )))?
            .filter_map(|entry| entry.ok())
            .filter(|p| p.is_file())
            .collect();

        if files.is_empty() {
            return Err(PipelineError::NoInputFiles(self.input_dir.clone()));
        }

        files.sort();
        Ok(files)
    }

    /// Get output path for a given FASTA file and track name.
    pub fn output_path(&self, fasta: &Path, track: &str) -> PathBuf {
        let stem = fasta.file_stem().unwrap_or_default().to_string_lossy();
        self.output_dir.join(format!("{}-vs-{}.bed", track, stem))
    }
}
