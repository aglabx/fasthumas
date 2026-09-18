# FastHumAS

A fast, single-binary Rust engine for centromeric alpha-satellite Higher-Order Repeat (HOR)
and Suprachromosomal Family (SF) annotation.

It replaces the legacy [HumAS-HMMER](https://github.com/enigene/HumAS-HMMER) / [HumAS-HMMER_for_AnVIL](https://github.com/fedorrik/HumAS-HMMER_for_AnVIL)
shell pipeline (`nhmmer` → `awk` → `sort` → `bedmap` → `awk` → `overlap_filter.py`). Everything
runs in memory in pure Rust: no external dependencies, no intermediate temporary files, and no
subprocess round-trips.

> **Key differences from the original bash pipeline:**
> - **Zero external dependencies:** Pure Rust, no `nhmmer`, `bedtools`, `awk`, or `python` needed on `PATH`.
> - **Minus-strand coordinate fix:** The original pipeline shifted minus-strand intervals by 1 bp to the right and made them 2 bp too short. FastHumAS fixes this, restoring exact forward/reverse-complement symmetry.
> - **MEA junction gap healing:** Closes the 1–3 bp artificial gaps left by HMMER alignment trimming between adjacent tandem monomers, meeting flush (>98% closed).
> - **Combined HOR + SF mode:** Can scan both HOR and SF profiles in a single pass, directly producing all four standard T2T CenSat BED9 tracks.
> - See [Differences from the original bash pipeline](#differences-from-the-original-bash-pipeline) for details.

## Install

Requires a Rust toolchain (Rust ≥ 1.74):

```bash
git clone https://github.com/aglabx/fasthumas.git
cd fasthumas
cargo build --release
# compiled binary is at ./target/release/fasthumas
```

HMM profiles can be obtained from [HumAS-HMMER_for_AnVIL](https://github.com/fedorrik/HumAS-HMMER_for_AnVIL):
- `AS-HORs-hmmer3.3.2-120124.hmm` (HOR models)
- `AS-SFs-hmmer3.0.290621.hmm` (SF models, 39 classes)

## Usage

### 1. Combined HOR + SF mode (Recommended)

Pass both `--hor` and `--sf` to annotate Higher-Order Repeats and Suprachromosomal Families simultaneously:

```bash
fasthumas \
  -i genome.fa \
  -t 48 \
  --hor AS-HORs-hmmer3.3.2-120124.hmm \
  --sf AS-SFs-hmmer3.0.290621.hmm \
  -o results/genome
```

This writes all four standard T2T CenSat tracks directly:

- `results/genome.AS-HOR+SF.bed` — composite annotation (HOR hits with SF fallback)
- `results/genome.AS-HOR.bed` — HOR-only annotation (SF-only monomer calls filtered out)
- `results/genome.AS-SF.bed` — SF monomer classification across the genome
- `results/genome.AS-strand.bed` — strand-resolved track (Blue `+`, Red `-`)

### 2. Single-profile mode

Scan against a single profile:

```bash
fasthumas -i genome.fa --hmm AS-HORs-hmmer3.3.2-120124.hmm -o results/genome_hor.bed -t 48
fasthumas -i genome.fa --hmm AS-SFs-hmmer3.0.290621.hmm    -o results/genome_sf.bed  -t 48
```

### Command-line options

| Flag | Default | Description |
|---|---|---|
| `-i`, `--input` | *required* | Input FASTA file (whole genome or extracted regions) |
| `--hor` | — | Path to HOR HMM profile (used with `--sf`) |
| `--sf` | — | Path to SF HMM profile (used with `--hor`) |
| `--hmm` | — | Path to single HMM profile (alternative to `--hor` + `--sf`) |
| `-o`, `--output` | input file stem | Output path (single mode) or prefix (combined mode) |
| `-t`, `--threads` | `48` | Total worker threads across all parallel jobs |
| `--score-threshold` | `0.7` | Minimum score-to-length ratio for a hit to be kept |

Set `RUST_LOG=debug` or `RUST_LOG=info` for per-sequence runtime details.

### What to feed it

You can feed FastHumAS:
- **A whole-genome FASTA** directly (e.g. T2T-CHM13, HG002, RPE1).
- **Extracted alpha-satellite regions** (e.g. from CenSat annotation: `hor`, `dhor`, and `mon` intervals).

Sequences are processed concurrently using Rayon in-memory work-stealing, scaling linearly with available CPU cores.

### Output format

Standard BED9 format. Column 9 is an RGB triplet from the standardized 218-entry HOR color table, ready for immediate visualization in the UCSC Genome Browser or IGV.

```tsv
chr1    121535434    121535605    S1.13/21    1000    +    121535434    121535605    0,210,0
```

---

## Differences from the original bash pipeline

### Minus-strand coordinates are corrected

`nhmmer` reports alignment coordinates in query order, so on the minus strand `alifrom > alito` — column 7 is the **end** of the hit on the target, not its start. `hmmertblout2bed.awk` computes:
```awk
aliFrom = $7 - 1; aliTo = $8;
...
print $3 "\t" aliTo "\t" aliFrom ...
```
Because `$7` is the upper coordinate on the minus strand, putting the `- 1` there adjusts the wrong end, leaving the start unadjusted. Every minus-strand interval therefore comes out shifted 1 base to the right and **2 bases too short**.

In the original pipeline output over 34 annotated CDR regions (12,568 records, 4,255 on the minus strand), tandem monomers butt-joined on the plus strand met exactly, while the same junctions on the minus strand were all 2 bp apart. Median monomer length was 171 bp on the plus strand and 169 bp on the minus strand.

FastHumAS converts coordinates as `[min(alifrom, alito) - 1, max(alifrom, alito))` on both strands, so forward and reverse-complement alignments yield identical length and symmetric boundaries.

### Gaps left by MEA trimming are closed

Profile HMM alignments constructed by maximum expected accuracy (MEA) drop residues whose posterior probability is low, regardless of emission score. At a tandem repeat junction between near-identical monomers, register ambiguity causes terminal residues to be trimmed. This leaves an artificial gap of 1–3 bases between monomers that are in fact flush in the sequence.

FastHumAS inspects the model coordinates (`hmmfrom`, `hmmto`):
1. Calculates the number of trimmed positions at each end.
2. Hands bases back to the monomer if:
   - **The model expects them:** the bases match the model's consensus (`hmmemit -c`) at the trimmed positions (reverse-complemented for minus-strand hits) at ≥ 50% identity; or
   - **The array vouches for them:** the same base recurs at that model position at least 10 times across the run at ≥ 50% local frequency.

Because the two sides of a junction own only the bases in the gap between them, this gap closure never produces overlaps.

### Other improvements

- **In-memory execution:** Direct Plan7 HMM parser and in-memory bounded banded DP alignment; zero temporary disk files.
- **Deterministic colors:** The AWK `getColor()` iterated with `for (i in cnA)` whose order is undefined in POSIX awk. FastHumAS matches rules deterministically in table source order.
- **Sequence boundary reset:** `overlap_filter.py` carried state across chromosome boundaries, which allowed the end of one chromosome to suppress hits at the start of the next. FastHumAS strictly scopes near-dedup to individual sequences.
- **Symmetric overlap window:** `overlap_filter.py` used `start in range(prev_start-10, prev_start+10)` (window `[-10, +9]`). FastHumAS uses symmetric `abs(delta) < 10`.

---

## Validation

Tested against the legacy bash pipeline on 34 CDR regions of a human haplotype assembly (2.1 Mb, 12,568 records):

| Check | Result |
|---|---|
| Plus strand vs bash pipeline (SF profile, 12,518 records) | Identical record for record |
| Plus strand vs original pipeline (HOR profile, 8,313 records) | Identical except 1 record where two equally scoring models tie |
| Minus strand vs original pipeline (4,255 records) | Differs by design: 1 bp wider at each end, correcting the legacy 2 bp truncation |
| Flush inter-monomer junctions (HOR profile) | Original: 58.1% (7,286 / 12,534) → **FastHumAS: 98.5% (12,350 / 12,534)** |
| Flush inter-monomer junctions (SF profile) | Original: 35.3% (4,428 / 12,531) → **FastHumAS: 87.6% (10,975 / 12,531)** |

---

## Development

```bash
cargo test
cargo clippy --all-targets -- -D warnings
cargo fmt --check
```

---

## Attribution and citation

This tool is based on the algorithm and methodology of **HumAS-HMMER**:
> Uralsky LI, Shepelev VA, Alexandrov AA, Yurov YB, Rogaev EI, Alexandrov IA (2019).
> Classification and monomer-by-monomer annotation dataset of suprachromosomal family 1 alpha satellite
> higher-order repeats in hg38 human genome assembly. *Data in Brief* 24:103708.
> [doi:10.1016/j.dib.2019.103708](https://doi.org/10.1016/j.dib.2019.103708)

and HMMER:
> Eddy SR (2011). Accelerated profile HMM searches. *PLoS Comput Biol* 7(10):e1002195.
> [doi:10.1371/journal.pcbi.1002195](https://doi.org/10.1371/journal.pcbi.1002195)

See [NOTICE.md](NOTICE.md) for provenance details.

## License

MIT License — see [LICENSE](LICENSE).
