use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use log::info;
use rayon::prelude::*;

use crate::bed::{BedRecord, Rgb};
use crate::config::{Config, RunMode};
use crate::consensus::Consensus;
use crate::error::PipelineError;
use crate::fasta::{extract_intervals_in_memory, read_fasta_in_memory, MemoryFastaSeq};
use crate::hmm::HmmCollection;
use crate::junction::{self, GapStats, Side, Tally, Verdict, POLICY};
use crate::scanner::{scan_sequence_in_memory, MasterProfile};

/// Bases owed to a record, accumulated on each side.
type Edits = HashMap<usize, (u64, u64)>;

/// What one sequence produced, before any junction is acted on.
#[derive(Debug, Default)]
struct Pending {
    records: Vec<BedRecord>,
    sides: Vec<Side>,
    stats: GapStats,
}

fn track_path(prefix: &Path, suffix: &str) -> PathBuf {
    let stem = prefix.file_name().unwrap_or_default().to_string_lossy();
    prefix.with_file_name(format!("{}.{}.bed", stem, suffix))
}

fn write_bed_file<P: AsRef<Path>>(path: P, records: &[BedRecord]) -> Result<(), PipelineError> {
    let path = path.as_ref();
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)?;
        }
    }
    let mut out = BufWriter::new(File::create(path)?);
    for rec in records {
        writeln!(out, "{}", rec.to_bed9_string())?;
    }
    out.flush()?;
    Ok(())
}

/// Run in-memory HMM profile scan and junction healing for a slice of FASTA records.
pub fn annotate_sequences(
    seqs: &[MemoryFastaSeq],
    hmm_path: &Path,
    score_threshold: f64,
) -> Result<(Vec<BedRecord>, GapStats), PipelineError> {
    let start_hmm = std::time::Instant::now();
    let coll = HmmCollection::load_from_file(hmm_path)?;
    let master = MasterProfile::from_collection(&coll);
    let consensus = Consensus::from_collection(&coll);
    info!(
        "{}: loaded {} models ({} unique names) in {:.2}s",
        hmm_path.display(),
        coll.models.len(),
        consensus.names(),
        start_hmm.elapsed().as_secs_f64()
    );

    let start_scan = std::time::Instant::now();
    let pendings: Vec<Pending> = seqs
        .par_iter()
        .map(|seq| {
            let records =
                scan_sequence_in_memory(&seq.name, &seq.seq, &coll, &master, score_threshold);

            let (candidates, stats) = junction::candidates(&records, &consensus, &POLICY);
            let mut sides = Vec::new();
            if !candidates.is_empty() {
                let intervals: Vec<(u64, u64)> =
                    candidates.iter().map(|c| (c.start, c.end)).collect();
                let bases = extract_intervals_in_memory(&seq.seq, &intervals);
                sides = junction::sides(&records, &candidates, &bases, &consensus);
            }

            Pending {
                records,
                sides,
                stats,
            }
        })
        .collect();

    info!(
        "{}: in-memory scan completed in {:.2}s",
        hmm_path.display(),
        start_scan.elapsed().as_secs_f64()
    );

    let mut tally = Tally::default();
    for pending in &pendings {
        for side in &pending.sides {
            tally.add(side);
        }
    }

    let mut total_stats = GapStats::default();
    let mut all_records = Vec::new();

    for pending in &pendings {
        total_stats.merge(&pending.stats);
        let mut edits = Edits::new();

        for side in &pending.sides {
            let verdict = tally.verdict(side, &POLICY);
            total_stats.sides += 1;
            match verdict {
                Verdict::Consensus => total_stats.by_consensus += 1,
                Verdict::Recurrent { .. } => total_stats.by_recurrence += 1,
                Verdict::Rejected => total_stats.rejected += 1,
            }
            if verdict.accepted() {
                total_stats.bases += side.bases;
                let slot = edits.entry(side.record).or_insert((0, 0));
                if side.extend_right {
                    slot.1 += side.bases;
                } else {
                    slot.0 += side.bases;
                }
            }
        }

        for (rec_idx, rec) in pending.records.iter().enumerate() {
            let mut modified = rec.clone();
            if let Some(&(left, right)) = edits.get(&rec_idx) {
                modified.start = modified.start.saturating_sub(left);
                modified.end += right;
                modified.thick_start = modified.start;
                modified.thick_end = modified.end;
            }
            all_records.push(modified);
        }
    }

    info!(
        "{}: junctions: {} gaps <= {} bp, {} accounted for by trimming, {} unresolved; sides: {} by consensus, {} by recurring variant, {} left alone; {} bases returned",
        hmm_path.display(),
        total_stats.candidates, POLICY.max_gap, total_stats.accounted, total_stats.unresolved,
        total_stats.by_consensus, total_stats.by_recurrence, total_stats.rejected, total_stats.bases
    );

    Ok((all_records, total_stats))
}

/// Run legacy nhmmer backend on the input FASTA.
/// Splits FASTA per sequence into temp_dir, runs nhmmer concurrently with proportional threading,
/// parses hits in memory, heals junctions, and returns final BedRecords.
pub fn annotate_sequences_legacy(
    config: &Config,
    hmm_path: &Path,
    temp_dir: &Path,
) -> Result<(Vec<BedRecord>, GapStats), PipelineError> {
    let start_legacy = std::time::Instant::now();
    let split = crate::fasta::split_fasta(&config.input, temp_dir)?;
    let n = split.seqs.len();
    let dbsize_mb = (split.residues as f64 / 1.0e6).max(1.0e-6);
    let lengths: Vec<u64> = split.seqs.iter().map(|s| s.len).collect();
    let plan = config.plan(&lengths);

    info!(
        "{}: {} sequences, {:.1} Mb; {} concurrent legacy nhmmer jobs, {}-{} threads each",
        config.input.display(),
        n,
        dbsize_mb,
        plan.jobs,
        plan.threads.iter().min().copied().unwrap_or(0),
        plan.threads.iter().max().copied().unwrap_or(0),
    );

    let consensus = Consensus::from_profile(hmm_path)?;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(plan.jobs)
        .build()
        .map_err(|e| PipelineError::Io(std::io::Error::other(e.to_string())))?;

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by_key(|&i| std::cmp::Reverse(lengths[i]));

    let outcomes: Vec<(usize, Result<Pending, String>)> = pool.install(|| {
        order
            .par_iter()
            .map(|&i| {
                let seq = &split.seqs[i];
                let threads = plan.threads[i];
                let records = match crate::nhmmer::run_nhmmer(
                    hmm_path,
                    &seq.path,
                    threads,
                    dbsize_mb,
                    config.score_threshold,
                ) {
                    Ok(r) => r,
                    Err(e) => return (i, Err(format!("{}: {}", seq.name, e))),
                };

                let (candidates, stats) = junction::candidates(&records, &consensus, &POLICY);
                let mut sides = Vec::new();
                if !candidates.is_empty() {
                    let intervals: Vec<(u64, u64)> =
                        candidates.iter().map(|c| (c.start, c.end)).collect();
                    let bases = match crate::fasta::read_intervals(&seq.path, &intervals) {
                        Ok(b) => b,
                        Err(e) => return (i, Err(format!("{}: {}", seq.name, e))),
                    };
                    sides = junction::sides(&records, &candidates, &bases, &consensus);
                }

                if !config.keep_temp {
                    let _ = fs::remove_file(&seq.path);
                }

                (
                    i,
                    Ok(Pending {
                        records,
                        sides,
                        stats,
                    }),
                )
            })
            .collect()
    });

    let mut pendings: Vec<Pending> = Vec::with_capacity(n);
    let mut slots: Vec<Option<Pending>> = (0..n).map(|_| None).collect();
    let mut failures: Vec<String> = Vec::new();

    for (i, outcome) in outcomes {
        match outcome {
            Ok(pending) => slots[i] = Some(pending),
            Err(msg) => failures.push(msg),
        }
    }

    if !failures.is_empty() {
        return Err(PipelineError::NhmmerFailed {
            path: config.input.clone(),
            message: format!(
                "{} sequence(s) failed in nhmmer:\n{}",
                failures.len(),
                failures.join("\n")
            ),
        });
    }

    for p in slots.into_iter().flatten() {
        pendings.push(p);
    }

    let mut tally = Tally::default();
    for pending in &pendings {
        for side in &pending.sides {
            tally.add(side);
        }
    }

    let mut total_stats = GapStats::default();
    let mut all_records = Vec::new();

    for pending in &pendings {
        total_stats.merge(&pending.stats);
        let mut edits = Edits::new();

        for side in &pending.sides {
            let verdict = tally.verdict(side, &POLICY);
            total_stats.sides += 1;
            match verdict {
                Verdict::Consensus => total_stats.by_consensus += 1,
                Verdict::Recurrent { .. } => total_stats.by_recurrence += 1,
                Verdict::Rejected => total_stats.rejected += 1,
            }
            if verdict.accepted() {
                total_stats.bases += side.bases;
                let slot = edits.entry(side.record).or_insert((0, 0));
                if side.extend_right {
                    slot.1 += side.bases;
                } else {
                    slot.0 += side.bases;
                }
            }
        }

        for (rec_idx, rec) in pending.records.iter().enumerate() {
            let mut modified = rec.clone();
            if let Some(&(left, right)) = edits.get(&rec_idx) {
                modified.start = modified.start.saturating_sub(left);
                modified.end += right;
                modified.thick_start = modified.start;
                modified.thick_end = modified.end;
            }
            all_records.push(modified);
        }
    }

    info!(
        "{}: legacy nhmmer run finished: {} records in {:.2}s",
        hmm_path.display(),
        all_records.len(),
        start_legacy.elapsed().as_secs_f64()
    );

    Ok((all_records, total_stats))
}

/// FastHumAS execution entry point.
///
/// Supports two backends:
/// - Fast in-memory pure-Rust scanner (default)
/// - Legacy nhmmer backend via `--legacy`
pub fn run(config: &Config) -> Result<(), PipelineError> {
    let start_total = std::time::Instant::now();

    let (seqs, temp_dir) = if config.legacy {
        fs::create_dir_all(&config.temp_base)?;
        let temp_dir = config.temp_dir_path();
        fs::create_dir_all(&temp_dir)?;
        info!(
            "Running in --legacy mode with nhmmer backend (temp dir: {})",
            temp_dir.display()
        );
        (None, Some(temp_dir))
    } else {
        let start_io = std::time::Instant::now();
        let s = read_fasta_in_memory(&config.input)?;
        let total_residues: usize = s.iter().map(|item| item.seq.len()).sum();
        info!(
            "{}: loaded {} sequences, {:.2} Mb in {:.2}s",
            config.input.display(),
            s.len(),
            total_residues as f64 / 1.0e6,
            start_io.elapsed().as_secs_f64()
        );
        (Some(s), None)
    };

    let annotate = |hmm: &Path| -> Result<(Vec<BedRecord>, GapStats), PipelineError> {
        if config.legacy {
            annotate_sequences_legacy(config, hmm, temp_dir.as_ref().unwrap())
        } else {
            annotate_sequences(seqs.as_ref().unwrap(), hmm, config.score_threshold)
        }
    };

    let res = match config.mode {
        RunMode::Single {
            ref hmm,
            ref output,
        } => {
            let (records, _) = annotate(hmm)?;
            write_bed_file(output, &records)?;
            info!(
                "{} records -> {} (total time: {:.2}s)",
                records.len(),
                output.display(),
                start_total.elapsed().as_secs_f64()
            );
            Ok(())
        }
        RunMode::Combined {
            ref hor,
            ref sf,
            ref prefix,
        } => {
            if let Some(ref hor_path) = hor {
                let (hor_records, _) = annotate(hor_path)?;

                let hor_sf_path = track_path(prefix, "AS-HOR+SF");
                write_bed_file(&hor_sf_path, &hor_records)?;
                info!("{} records -> {}", hor_records.len(), hor_sf_path.display());

                let hor_path_out = track_path(prefix, "AS-HOR");
                let hor_only: Vec<BedRecord> = hor_records
                    .into_iter()
                    .filter(|r| r.name.len() != 2)
                    .collect();
                write_bed_file(&hor_path_out, &hor_only)?;
                info!("{} records -> {}", hor_only.len(), hor_path_out.display());
            }

            if let Some(ref sf_path) = sf {
                let (sf_records, _) = annotate(sf_path)?;

                let sf_path_out = track_path(prefix, "AS-SF");
                write_bed_file(&sf_path_out, &sf_records)?;
                info!("{} records -> {}", sf_records.len(), sf_path_out.display());

                let strand_path_out = track_path(prefix, "AS-strand");
                let strand_records: Vec<BedRecord> = sf_records
                    .into_iter()
                    .map(|mut r| {
                        if r.strand == '+' {
                            r.color = Rgb(0, 0, 255);
                        } else if r.strand == '-' {
                            r.color = Rgb(255, 0, 0);
                        }
                        r
                    })
                    .collect();
                write_bed_file(&strand_path_out, &strand_records)?;
                info!(
                    "{} records -> {}",
                    strand_records.len(),
                    strand_path_out.display()
                );
            }

            info!(
                "Combined HOR+SF pipeline completed in {:.2}s",
                start_total.elapsed().as_secs_f64()
            );
            Ok(())
        }
    };

    if let Some(ref td) = temp_dir {
        if !config.keep_temp {
            let _ = fs::remove_dir_all(td);
        } else {
            info!("Temporary files preserved in: {}", td.display());
        }
    }

    res
}
