use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use log::{debug, error, info, warn};
use rayon::prelude::*;

use crate::bed::BedRecord;
use crate::config::Config;
use crate::consensus::Consensus;
use crate::error::PipelineError;
use crate::fasta::{read_intervals, split_fasta, SplitSeq};
use crate::junction::{self, GapStats, Side, Tally, Verdict, POLICY};
use crate::nhmmer::run_nhmmer;

/// Bases owed to a record, accumulated on each side.
type Edits = HashMap<usize, (u64, u64)>;

/// What one sequence produced, before any junction is acted on.
#[derive(Debug, Default)]
struct Pending {
    sides: Vec<Side>,
    stats: GapStats,
}

/// Annotate one FASTA file with one profile and write the BED9 track.
///
/// The input is split into one file per sequence and those are annotated
/// concurrently. Junctions are only *described* while that runs: deciding
/// whether a trimmed base is believable needs to know how often the same base
/// turns up at the same model position across the whole input, so the decision
/// waits until every sequence has been seen.
pub fn run(config: &Config) -> Result<(), PipelineError> {
    if let Some(parent) = config.output.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)?;
        }
    }

    fs::create_dir_all(&config.temp_base)?;
    let temp_dir = config.temp_dir_path();
    // create_dir, not create_dir_all: if it already exists we do not own it,
    // and we must not delete a directory we did not create.
    fs::create_dir(&temp_dir)?;

    let result = run_in(config, &temp_dir);

    match (&result, config.keep_temp) {
        (Ok(()), false) => fs::remove_dir_all(&temp_dir)?,
        (Ok(()), true) => info!("temporary files kept in {}", temp_dir.display()),
        (Err(_), _) => warn!(
            "temporary files left in {} for inspection",
            temp_dir.display()
        ),
    }

    result
}

fn run_in(config: &Config, temp_dir: &Path) -> Result<(), PipelineError> {
    let split = split_fasta(&config.input, temp_dir)?;
    let n = split.seqs.len();

    // Pin nhmmer's database size to the whole input so that splitting it does
    // not shift E-values, and with them which hits get reported at all.
    let dbsize_mb = (split.residues as f64 / 1.0e6).max(1.0e-6);
    let (jobs, threads_per_job) = config.job_layout(n);

    info!(
        "{}: {} sequences, {:.1} Mb; {} concurrent jobs x {} nhmmer threads",
        config.input.display(),
        n,
        dbsize_mb,
        jobs,
        threads_per_job
    );

    let consensus = Consensus::from_profile(&config.hmm)?;
    info!(
        "{}: {} model names ({} covering models of more than one length)",
        config.hmm.display(),
        consensus.names(),
        consensus.length_ambiguous()
    );

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(jobs)
        .build()
        .map_err(|e| PipelineError::Io(std::io::Error::other(e.to_string())))?;

    // --- annotate every sequence, describing junctions but not acting yet ---
    let outcomes: Vec<Result<Pending, String>> = pool.install(|| {
        split
            .seqs
            .par_iter()
            .enumerate()
            .map(|(i, seq)| {
                annotate(
                    config,
                    temp_dir,
                    i,
                    seq,
                    threads_per_job,
                    dbsize_mb,
                    &consensus,
                )
                .map_err(|e| format!("{}: {}", seq.name, e))
            })
            .collect()
    });

    let mut pendings: Vec<Pending> = Vec::with_capacity(n);
    let mut failures: Vec<String> = Vec::new();
    for outcome in outcomes {
        match outcome {
            Ok(p) => pendings.push(p),
            Err(message) => {
                error!("{}", message);
                failures.push(message);
            }
        }
    }
    if !failures.is_empty() {
        return Err(PipelineError::NhmmerFailed {
            path: config.input.clone(),
            message: format!(
                "{} of {} sequences failed:\n{}",
                failures.len(),
                n,
                failures.join("\n")
            ),
        });
    }

    // --- count what the whole run saw at each model position ---
    let mut tally = Tally::default();
    for pending in &pendings {
        for side in &pending.sides {
            tally.add(side);
        }
    }

    // --- decide, apply, and write ---
    let mut out = BufWriter::new(File::create(&config.output)?);
    let mut stats = GapStats::default();
    let mut records = 0usize;

    for (index, pending) in pendings.iter().enumerate() {
        stats.merge(&pending.stats);
        let mut edits = Edits::new();

        for side in &pending.sides {
            let verdict = tally.verdict(side, &POLICY);
            stats.sides += 1;
            match verdict {
                Verdict::Consensus => stats.by_consensus += 1,
                Verdict::Recurrent { .. } => stats.by_recurrence += 1,
                Verdict::Rejected => stats.rejected += 1,
            }
            if verdict.accepted() {
                stats.bases += side.bases;
                let slot = edits.entry(side.record).or_insert((0, 0));
                if side.extend_right {
                    slot.1 += side.bases;
                } else {
                    slot.0 += side.bases;
                }
            }
        }

        records += emit(&part_path(temp_dir, index), &edits, &mut out)?;
    }
    out.flush()?;

    info!(
        "junctions: {} gaps <= {} bp, {} accounted for by trimming, {} unresolved; sides: {} by consensus, {} by recurring variant, {} left alone; {} bases returned",
        stats.candidates, POLICY.max_gap, stats.accounted, stats.unresolved,
        stats.by_consensus, stats.by_recurrence, stats.rejected, stats.bases
    );
    info!("{} records -> {}", records, config.output.display());

    Ok(())
}

/// Scan one sequence, write the raw records, and describe the junctions without
/// touching them.
#[allow(clippy::too_many_arguments)]
fn annotate(
    config: &Config,
    temp_dir: &Path,
    index: usize,
    seq: &SplitSeq,
    threads: usize,
    dbsize_mb: f64,
    consensus: &Consensus,
) -> Result<Pending, PipelineError> {
    debug!("{}: starting", seq.name);

    let records = run_nhmmer(
        &config.hmm,
        &seq.path,
        threads,
        dbsize_mb,
        config.score_threshold,
    )?;

    let (candidates, stats) = junction::candidates(&records, consensus, &POLICY);
    let mut sides = Vec::new();
    if !candidates.is_empty() {
        let intervals: Vec<(u64, u64)> = candidates.iter().map(|c| (c.start, c.end)).collect();
        let bases = read_intervals(&seq.path, &intervals)?;
        sides = junction::sides(&records, &candidates, &bases, consensus);
    }

    write_part(&part_path(temp_dir, index), &records)?;

    // The split FASTA has served its purpose — the junction pass was the last
    // thing that needed the sequence itself.
    let _ = fs::remove_file(&seq.path);
    Ok(Pending { sides, stats })
}

/// Stream one sequence's records, widening the ones that earned it.
fn emit(part: &Path, edits: &Edits, out: &mut BufWriter<File>) -> Result<usize, PipelineError> {
    let reader = BufReader::new(File::open(part)?);
    let mut written = 0usize;

    for (index, line) in reader.lines().enumerate() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        match edits.get(&index) {
            None => writeln!(out, "{}", line)?,
            Some((left, right)) => {
                let mut fields: Vec<String> = line.split('\t').map(str::to_string).collect();
                if fields.len() < 9 {
                    continue;
                }
                widen(&mut fields, *left, *right);
                writeln!(out, "{}", fields.join("\t"))?;
            }
        }
        written += 1;
    }

    let _ = fs::remove_file(part);
    Ok(written)
}

/// Move a BED9 record's edges out by the bases it was owed.
fn widen(fields: &mut [String], left: u64, right: u64) {
    let start: u64 = fields[1].parse().unwrap_or(0);
    let end: u64 = fields[2].parse().unwrap_or(0);
    fields[1] = start.saturating_sub(left).to_string();
    fields[2] = (end + right).to_string();
    fields[6] = fields[1].clone();
    fields[7] = fields[2].clone();
}

fn part_path(temp_dir: &Path, index: usize) -> PathBuf {
    temp_dir.join(format!("seq_{:06}.bed", index))
}

fn write_part(path: &Path, records: &[BedRecord]) -> Result<(), PipelineError> {
    let mut content = String::with_capacity(records.len() * 80);
    for rec in records {
        content.push_str(&rec.to_bed9_string());
        content.push('\n');
    }
    fs::write(path, content)?;
    Ok(())
}
