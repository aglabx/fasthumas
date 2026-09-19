# FastHumAS Benchmark & Reproducibility Package

This directory contains the machine logs, run scripts, manifests, and validation tools supporting the empirical performance and accuracy claims of the FastHumAS manuscript.

## Environment & Hardware
- **Server:** `aglab0`
- **CPU Architecture:** Dual-socket AMD EPYC 9654 (total 192 physical cores, 384 threads @ 2.40–3.70 GHz)
- **Allocated Worker Threads:** 96 threads per FastHumAS run (`-t 96`)
- **System Memory:** 2.9 TiB DDR5-4800 ECC (3,019,898,880 KB)
- **Operating System:** Ubuntu 22.04 LTS (Linux kernel 6.8.0-40-generic x86_64)
- **HMM Profile Libraries:**
  - `AS-HORs-hmmer3.3.2-120124.hmm`: 1,963 monomer-level HOR profile models (lengths 88–176 bp; MD5: `4c642b44df914a47bed2cd3d0998f8e6`)
  - `AS-SFs-hmmer3.0.290621.hmm`: 39 SF profile models representing suprachromosomal family classes SF1–SF5 (lengths 167–174 bp, mean 170.8 bp, 21 models exactly 171 bp; MD5: `b7c211e2e72edfbe9f9068a5432f4025`)

## Directory Contents

### Manifests and Automated Runners
- `manifest.tsv`: Complete registry of the 7 evaluated primate T2T reference assemblies, including source accessions, sequence counts, exact nucleotide lengths from `.fai`, and FASTA MD5 checksums.
- `run_fasthumas_suite.sh`: Automated multi-genome benchmark harness with isolated per-assembly result directories, strict exit status checking, and GNU `/usr/bin/time -v` resource logging.
- `fasthumas_summary.tsv`: Complete machine-generated summary of wall-clock times, CPU user/system times, peak resident memory (RSS in KB), and output record counts across all 7 evaluated assemblies.

### Empirical Validation & Accuracy Scripts
- `compare_fasthumas_vs_nhmmer.py`: Comprehensive concordance evaluator comparing FastHumAS with legacy HumAS-HMMER (nhmmer) across whole chromosomes.
  - **Results on 9 finished CHM13 chromosomes (174,925 legacy loci):**
    - Sensitivity (Recall): **99.37%** (173,827 / 174,925)
    - Precision: **99.96%** (173,827 / 173,899)
    - Strand agreement: **100.00%**
    - Subfamily / HOR label agreement: **99.28%**
    - Boundary concordance within $\le 3$ bp: **97.33%**; within $\le 10$ bp: **99.72%**
- `calibrate_scoring.py`: Calibration script evaluating the empirical 0.85 scaling factor between Gotoh PSSM log-odds sums and HMMER bit scores across 58,183 exact-boundary monomers with a strict chromosome split:
  - **Training set (chr1, chr3, chr8; $n = 23,527$):** Mean $0.8495 \pm 0.0074$ (SD), SEM $0.000048$, Median $0.8504$
  - **Independent validation set (chr10, 11, 12, 14, 22, Y; $n = 34,656$):** Mean $0.8483 \pm 0.0085$ (SD), SEM $0.000046$, Median $0.8500$
- `analyze_junction_healing.py`: Analysis of inter-monomer junction gaps:
  - Demonstrates that legacy HumAS-HMMER contained 99,561 micro-gaps (1–3 bp) across 9 chromosomes (56.9% of monomers), of which 68,211 (68.5%) were artificial 2-bp gaps on the minus strand resulting from the 1-bp coordinate offset in `hmmertblout2bed.awk` (~190,000–200,000 genome-wide).
  - FastHumAS eliminates coordinate shifts and heals trimming artifacts, increasing flush abutting junctions (gap = 0 bp) from 32.7% to 93.0% and reducing micro-gaps by 92.1% (from 99,561 to 7,911), while preserving authentic biological insertions.
- `compare_chm13_dp_fix.py`: Validates consistency before and after the Gotoh DP recurrence fix on CHM13 (488,636 concordant records, 99.976%).

### Machine Logs (`logs/`)
- Individual `.time.log` and `.pipeline.log` for all 7 assemblies:
  - `T2T-CHM13_v2.0` (Homo sapiens, 3,117.28 Mb, 38m 55s, 31.05 GB RSS)
  - `T2T-HG002_mat` (Homo sapiens, 3,051.51 Mb, 40m 44s, 30.44 GB RSS)
  - `T2T-RPE1_hap1` (Homo sapiens, 3,064.05 Mb, 47m 21s, 30.55 GB RSS)
  - `Chimpanzee_mPanTro3` (Pan troglodytes, 3,034.92 Mb, 39m 51s, 30.23 GB RSS)
  - `Bonobo_mPanPan1` (Pan paniscus, 3,210.03 Mb, 38m 56s, 30.74 GB RSS)
  - `Gorilla_mGorGor1` (Gorilla gorilla, 3,554.71 Mb, 1h 05m 14s, 35.28 GB RSS)
  - `Orangutan_mPonPyg2` (Pongo pygmaeus, 3,171.52 Mb, 1h 04m 52s, 31.76 GB RSS)
- `nhmmer_concordance.summary.txt`: Detailed per-chromosome and overall concordance statistics between FastHumAS and legacy nhmmer.
- `scoring_calibration.summary.txt`: Detailed statistics of Gotoh vs HMMER bit scores across training and validation chromosome splits.
- `paired_calibration_scores_sample.tsv`: 5,000 sampled paired monomer scores for direct independent inspection.
- `junction_healing_analysis.summary.txt`: Per-chromosome counts of flush junctions, 1–3 bp micro-gaps, and 2-bp minus-strand gaps.
- `chm13_dp_fix_diff.summary.txt`: Summary of the Gotoh DP recurrence fix validation.
- `chm13_dp_fix_discordant_loci.tsv`: Full classification of all 243 discordant loci between pre-fix and post-fix DP runs on CHM13.
