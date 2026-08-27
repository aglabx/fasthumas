mod bed;
mod color;
mod config;
mod consensus;
mod error;
mod fasta;
mod junction;
mod nhmmer;
mod overlap;
mod pipeline;
mod tblout;

use std::path::PathBuf;
use std::process;

use clap::Parser;

use config::{Config, Options};

/// HumAS-HMMER: alpha-satellite annotation.
///
/// Scans one FASTA file with one HMM profile — HOR or SF — and writes a
/// colour-coded BED9 track: hits filtered by score, overlaps resolved, and the
/// small gaps HMMER's alignment trimming leaves between adjacent monomers
/// closed where the sequence bears it out.
#[derive(Parser, Debug)]
#[command(name = "humas-hmmer", version, about)]
struct Cli {
    /// Input FASTA file (may contain any number of sequences)
    #[arg(short = 'i', long, value_name = "FASTA")]
    input: PathBuf,

    /// HMM profile to scan with, HOR or SF (e.g. AS-HORs-hmmer3.3.2-120124.hmm)
    #[arg(long, value_name = "HMM")]
    hmm: PathBuf,

    /// Output path; `.bed` is appended unless it is already there. If it names
    /// an existing directory the input file stem is used inside it. Defaults to
    /// the input file stem in the current directory
    #[arg(short = 'o', long, value_name = "PATH")]
    output: Option<PathBuf>,

    /// Total CPU budget across all concurrent nhmmer jobs
    #[arg(short = 't', long, default_value = "48")]
    threads: usize,

    /// Threads per nhmmer job; jobs run concurrently, THREADS / this at a time.
    /// Lower it for more parallelism, raise it to cap peak memory
    #[arg(long, default_value = "4", value_name = "N")]
    threads_per_job: usize,

    /// Where to put the scratch directory (a private subdirectory is created
    /// inside it and removed when the run finishes). Defaults to the output
    /// directory
    #[arg(long, value_name = "DIR")]
    temp_dir: Option<PathBuf>,

    /// Keep the scratch directory instead of deleting it
    #[arg(long)]
    keep_temp: bool,

    /// Score-to-length threshold for filtering hits (0.0-1.0)
    #[arg(long, default_value = "0.7")]
    score_threshold: f64,
}

fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info"))
        .format_timestamp_secs()
        .init();

    let cli = Cli::parse();

    let config = match Config::new(Options {
        input: cli.input,
        hmm: cli.hmm,
        output: cli.output,
        threads: cli.threads,
        threads_per_job: cli.threads_per_job,
        temp_dir: cli.temp_dir,
        keep_temp: cli.keep_temp,
        score_threshold: cli.score_threshold,
    }) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("Error: {}", e);
            process::exit(1);
        }
    };

    if let Err(e) = pipeline::run(&config) {
        eprintln!("Pipeline failed: {}", e);
        process::exit(1);
    }
}
