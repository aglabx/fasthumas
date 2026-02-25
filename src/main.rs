mod bed;
mod color;
mod config;
mod error;
mod nhmmer;
mod overlap;
mod pipeline;
mod tblout;

use std::path::PathBuf;
use std::process;

use clap::Parser;

use config::Config;

/// HumAS-HMMER: parallel alpha-satellite annotation pipeline.
///
/// Runs nhmmer with HOR and SF HMM profiles against FASTA files,
/// converts to BED9 with color coding, and filters overlaps.
#[derive(Parser, Debug)]
#[command(name = "humas-hmmer", version, about)]
struct Cli {
    /// Input directory containing .fa files
    #[arg(short = 'i', long)]
    input_dir: PathBuf,

    /// Path to HOR HMM profile (e.g., AS-HORs-hmms-v3.0.hmm)
    #[arg(long)]
    hor_hmm: PathBuf,

    /// Path to SF HMM profile (e.g., AS-SFs-hmms-v3.0.hmm)
    #[arg(long)]
    sf_hmm: PathBuf,

    /// Output directory for BED files
    #[arg(short = 'o', long, default_value = ".")]
    output_dir: PathBuf,

    /// Total number of CPU threads to use
    #[arg(short = 't', long, default_value = "48")]
    threads: usize,

    /// Score-to-length threshold for filtering hits (0.0-1.0)
    #[arg(long, default_value = "0.7")]
    score_threshold: f64,
}

fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info"))
        .format_timestamp_secs()
        .init();

    let cli = Cli::parse();

    let config = match Config::new(
        cli.input_dir,
        cli.hor_hmm,
        cli.sf_hmm,
        cli.output_dir,
        cli.threads,
        cli.score_threshold,
    ) {
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
