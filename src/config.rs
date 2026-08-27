use std::path::{Path, PathBuf};

use crate::error::PipelineError;

/// Everything the CLI collects, before validation.
#[derive(Debug, Clone)]
pub struct Options {
    pub input: PathBuf,
    pub hmm: PathBuf,
    pub output: Option<PathBuf>,
    pub threads: usize,
    pub threads_per_job: usize,
    pub temp_dir: Option<PathBuf>,
    pub keep_temp: bool,
    pub score_threshold: f64,
}

#[derive(Debug, Clone)]
pub struct Config {
    /// The FASTA file to annotate. May hold any number of sequences.
    pub input: PathBuf,
    /// The HMM profile to scan with — HOR or SF, one per run.
    pub hmm: PathBuf,
    /// The BED file to write.
    pub output: PathBuf,
    /// Total CPU budget across all concurrent nhmmer jobs.
    pub threads: usize,
    /// Threads handed to a single nhmmer job.
    pub threads_per_job: usize,
    /// Directory the private scratch directory is created in.
    pub temp_base: PathBuf,
    pub keep_temp: bool,
    pub score_threshold: f64,
}

impl Config {
    pub fn new(o: Options) -> Result<Self, PipelineError> {
        if !o.input.is_file() {
            return Err(PipelineError::InputNotFile(o.input));
        }
        if !o.hmm.exists() {
            return Err(PipelineError::HmmNotFound(o.hmm));
        }

        let stem = o.input.file_stem().unwrap_or_default().to_os_string();
        let output = match o.output {
            // No -o: write `<input stem>.bed` into the current directory.
            None => with_extension(&PathBuf::from(&stem)),
            // -o names a directory: keep the input stem as the file name.
            Some(p) if p.is_dir() || ends_with_separator(&p) => with_extension(&p.join(&stem)),
            // -o is the path prefix itself.
            Some(p) => with_extension(&p),
        };

        let temp_base = o.temp_dir.unwrap_or_else(|| match output.parent() {
            Some(p) if !p.as_os_str().is_empty() => p.to_path_buf(),
            _ => PathBuf::from("."),
        });

        Ok(Config {
            input: o.input,
            hmm: o.hmm,
            output,
            threads: o.threads.max(1),
            threads_per_job: o.threads_per_job.max(1),
            temp_base,
            keep_temp: o.keep_temp,
            score_threshold: o.score_threshold,
        })
    }

    /// A scratch directory this process owns outright, so that deleting it at
    /// the end cannot touch anything the user put there.
    pub fn temp_dir_path(&self) -> PathBuf {
        self.temp_base
            .join(format!("humas-hmmer-tmp-{}", std::process::id()))
    }

    /// Split the CPU budget over the sequences at hand: never more jobs than
    /// sequences, and whatever threads that leaves go back into each job.
    pub fn job_layout(&self, sequences: usize) -> (usize, usize) {
        let sequences = sequences.max(1);
        let jobs = (self.threads / self.threads_per_job).clamp(1, sequences);
        let threads_per_job = (self.threads / jobs).max(1);
        (jobs, threads_per_job)
    }
}

/// `res/chm13` becomes `res/chm13.bed`; an explicit `.bed` is left alone.
fn with_extension(prefix: &Path) -> PathBuf {
    if prefix
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("bed"))
    {
        return prefix.to_path_buf();
    }
    let mut name = prefix.to_path_buf().into_os_string();
    name.push(".bed");
    PathBuf::from(name)
}

fn ends_with_separator(p: &Path) -> bool {
    p.as_os_str()
        .to_string_lossy()
        .ends_with(std::path::MAIN_SEPARATOR)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn config(threads: usize, per_job: usize) -> Config {
        Config {
            input: PathBuf::from("chm13.fa"),
            hmm: PathBuf::from("hor.hmm"),
            output: PathBuf::from("chm13.bed"),
            threads,
            threads_per_job: per_job,
            temp_base: PathBuf::from("."),
            keep_temp: false,
            score_threshold: 0.7,
        }
    }

    #[test]
    fn test_output_gets_a_bed_extension() {
        assert_eq!(
            with_extension(&PathBuf::from("res/chm13")),
            PathBuf::from("res/chm13.bed")
        );
        assert_eq!(
            with_extension(&PathBuf::from("res/chm13.bed")),
            PathBuf::from("res/chm13.bed")
        );
    }

    #[test]
    fn test_job_layout_full_budget() {
        // 96 threads, 4 per job, plenty of sequences: 24 jobs of 4.
        assert_eq!(config(96, 4).job_layout(34), (24, 4));
    }

    #[test]
    fn test_job_layout_gives_spare_threads_back() {
        // Only 10 sequences, so 10 jobs share all 96 threads.
        assert_eq!(config(96, 4).job_layout(10), (10, 9));
        // A single sequence gets the whole budget.
        assert_eq!(config(96, 4).job_layout(1), (1, 96));
    }

    #[test]
    fn test_job_layout_never_zero() {
        assert_eq!(config(1, 4).job_layout(34), (1, 1));
    }
}
