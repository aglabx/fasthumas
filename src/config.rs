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

/// nhmmer's own `--block_length` default: the target is read in blocks of
/// this many residues and one block goes to one worker thread.
const BLOCK_LENGTH: u64 = 1024 * 256;

/// How many jobs run at once, and how many threads each sequence gets.
#[derive(Debug, Clone, PartialEq)]
pub struct Layout {
    pub jobs: usize,
    /// One entry per sequence, in the order they were given.
    pub threads: Vec<usize>,
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

    /// Split the CPU budget over the sequences at hand.
    ///
    /// Sequences are wildly uneven — a centromeric HOR array can hold 86 % of
    /// an input's residues while a dozen small fields hold the rest — so an
    /// equal share per job leaves the whole run waiting on one starved job.
    /// Threads are handed out in proportion to length instead.
    ///
    /// Two bounds keep that honest. The share is normalised over the longest
    /// `jobs` sequences, which is the heaviest set that can run at once, so the
    /// budget holds when it matters. And no sequence is given more threads than
    /// it has work for: nhmmer reads a target in blocks of `BLOCK_LENGTH` and
    /// hands one block to a worker, so threads past that have nothing to do.
    pub fn plan(&self, lengths: &[u64]) -> Layout {
        let n = lengths.len().max(1);
        let jobs = (self.threads / self.threads_per_job).clamp(1, n);

        let mut longest: Vec<u64> = lengths.to_vec();
        longest.sort_unstable_by(|a, b| b.cmp(a));
        let window: u64 = longest.iter().take(jobs).sum::<u64>().max(1);

        let threads = lengths
            .iter()
            .map(|&len| {
                let share = (self.threads as u64).saturating_mul(len) / window;
                // Twice the block count, leaving room for parallelism finer
                // than the block loop rather than starving a job on a guess.
                let useful = len.div_ceil(BLOCK_LENGTH).saturating_mul(2).max(1);
                (share.clamp(1, useful) as usize).min(self.threads)
            })
            .collect();

        Layout { jobs, threads }
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
    fn test_plan_gives_a_lone_sequence_everything() {
        let plan = config(96, 4).plan(&[248_387_328]);
        assert_eq!(plan.jobs, 1);
        assert_eq!(plan.threads, vec![96]);
    }

    /// The array that holds most of the residues must not be left on the same
    /// thin share as the small fields beside it.
    #[test]
    fn test_plan_follows_length() {
        // chr1's alpha fields: one 4.5 Mb array and twelve small ones.
        let mut lengths = vec![4_504_439u64, 290_315, 120_272, 67_989, 54_629];
        lengths.extend([
            42_728u64, 31_100, 21_756, 21_198, 20_603, 16_051, 6_044, 5_686,
        ]);
        let plan = config(96, 4).plan(&lengths);

        assert_eq!(plan.jobs, 13);
        assert!(
            plan.threads[0] >= 8 * plan.threads[1],
            "the big array got {} against {} for the next one",
            plan.threads[0],
            plan.threads[1]
        );
        // Nobody is starved, and nobody exceeds the budget.
        assert!(plan.threads.iter().all(|&t| (1..=96).contains(&t)));
    }

    /// A sequence smaller than one nhmmer block has nothing to spread over
    /// more than a couple of threads, however much budget is going spare.
    #[test]
    fn test_plan_caps_at_the_work_available() {
        let plan = config(96, 4).plan(&[60_000, 60_000, 60_000]);
        assert!(
            plan.threads.iter().all(|&t| t <= 2),
            "handed out {:?} for single-block sequences",
            plan.threads
        );
    }

    #[test]
    fn test_plan_never_hands_out_zero() {
        let plan = config(1, 4).plan(&[10_000_000, 1]);
        assert!(plan.threads.iter().all(|&t| t >= 1));
        assert_eq!(plan.jobs, 1);
    }
}
