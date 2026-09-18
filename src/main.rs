mod bed;
mod color;
mod config;
mod consensus;
mod error;
mod fasta;
mod hmm;
mod junction;
mod nhmmer;
mod overlap;
mod pipeline;
mod scanner;
mod tblout;

use std::path::PathBuf;
use std::process;

use clap::Parser;

use config::{Config, Options};

/// FastHumAS: Ultra-fast in-memory centromeric alpha-satellite HOR/SF annotation engine.
///
/// Scans FASTA sequences using HMM profiles — HOR and/or SF — and writes
/// colour-coded BED9 tracks: hits filtered by score, overlaps resolved, and gaps
/// healed with 0 bp seamless boundaries.
#[derive(Parser, Debug)]
#[command(name = "fasthumas", version, about)]
struct Cli {
    /// Input FASTA file (may contain any number of sequences)
    #[arg(short = 'i', long, value_name = "FASTA")]
    input: PathBuf,

    /// Single HMM profile to scan with, HOR or SF (e.g. AS-HORs-hmmer3.3.2-120124.hmm)
    #[arg(long, value_name = "HMM")]
    hmm: Option<PathBuf>,

    /// Higher-Order Repeat (HOR) HMM profile (e.g. AS-HORs-hmmer3.3.2-120124.hmm)
    #[arg(long, value_name = "HOR_HMM")]
    hor: Option<PathBuf>,

    /// Superfamily (SF) HMM profile (e.g. AS-SFs-hmmer3.0.290621.hmm)
    #[arg(long, value_name = "SF_HMM")]
    sf: Option<PathBuf>,

    /// Output path or prefix. In single-HMM mode: output BED path. In combined HOR+SF mode: output prefix for standard tracks (.AS-HOR+SF.bed, .AS-HOR.bed, .AS-SF.bed, .AS-strand.bed)
    #[arg(short = 'o', long, value_name = "PATH")]
    output: Option<PathBuf>,

    /// Number of worker threads for parallel scanning
    #[arg(short = 't', long, default_value = "48")]
    threads: usize,

    /// Score-to-length threshold for filtering hits (0.0-1.0)
    #[arg(long, default_value = "0.7")]
    score_threshold: f64,

    /// Legacy option preserved for compatibility
    #[arg(long, default_value = "4", value_name = "N", hide = true)]
    threads_per_job: usize,

    /// Legacy option preserved for compatibility
    #[arg(long, value_name = "DIR", hide = true)]
    temp_dir: Option<PathBuf>,

    /// Legacy option preserved for compatibility
    #[arg(long, hide = true)]
    keep_temp: bool,
}

fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info"))
        .format_timestamp_secs()
        .init();

    let cli = Cli::parse();

    if cli.threads > 0 {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global();
    }

    let config = match Config::new(Options {
        input: cli.input,
        hmm: cli.hmm,
        hor: cli.hor,
        sf: cli.sf,
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
