# FastHumAS Benchmark & Reproducibility Package

This directory contains the machine logs, run scripts, manifests, and validation tools supporting the empirical performance and accuracy claims of the FastHumAS manuscript.

## Environment
- **Server:** `aglab0`
- **CPU:** Dual AMD EPYC 9654 (96 physical cores, 192 threads @ 2.40–3.70 GHz)
- **RAM:** 2.9 TiB DDR5-4800 ECC (3,019,898,880 KB)
- **OS:** Ubuntu 22.04 LTS (Linux kernel 6.8.0-40-generic x86_64)
- **Profiles:**
  - `AS-HORs-hmmer3.3.2-120124.hmm` (1,963 models, MD5: `4c642b44df914a47bed2cd3d0998f8e6`)
  - `AS-SFs-hmmer3.0.290621.hmm` (39 models, MD5: `b7c211e2e72edfbe9f9068a5432f4025`)

## Files in this Directory
- `manifest.tsv`: Complete list of the 7 evaluated T2T reference assemblies, with accession, sequence counts, exact nucleotide lengths from `.fai`, and MD5 checksums.
- `run_fasthumas_suite.sh`: The automated benchmark harness executing end-to-end multi-genome runs with GNU `/usr/bin/time -v`.
- `fasthumas_summary.tsv`: Machine-generated summary of wall-clock times, CPU user/system times, peak resident memory (RSS in KB), and track record counts.
- `compare_chm13_dp_fix.py`: Validation script demonstrating that 99.976% of CHM13 records are identical before and after the Gotoh DP recurrence fix, isolating changes to loci with genuine insertions.
- `logs/`: Directory containing raw machine `.time.log` and `.pipeline.log` records.
