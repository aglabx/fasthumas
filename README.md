# humas_hmmer

A single-binary Rust port of the [HumAS-HMMER](https://github.com/enigene/HumAS-HMMER) alpha-satellite
annotation pipeline. It scans FASTA files with `nhmmer` against alpha-satellite HOR and SF HMM profiles
and emits the four standard BED9 tracks — `AS-HOR+SF`, `AS-HOR`, `AS-SF`, `AS-strand`.

Everything except `nhmmer` itself runs in memory: no temp files, no `awk`/`sort`/`bedmap`/`python`
round-trips, and many genomes are processed in parallel.

> **Status: pre-release (v0.1.0).** The port is feature-complete and unit-tested, but output has not yet
> been diffed against the original bash pipeline on T2T-CHM13. See [Validation status](#validation-status).

## Why

The original pipeline is a chain of shell steps run once per FASTA file:

```
nhmmer → hmmertblout2bed.awk → sort → bedmap → awk → overlap_filter.py
```

Each step writes a temp file, and the driver script loops over inputs sequentially. On a corpus the size
of the VGP collection (~605 genomes) that is roughly a week of wall-clock time on a 48-core machine.
This port keeps `nhmmer` (optimized C with SIMD — there is nothing to gain by rewriting it) and replaces
everything around it with in-process Rust, running N `nhmmer` jobs concurrently over a `rayon` pool.

## Install

Requires a Rust toolchain and `nhmmer` (HMMER ≥ 3.3) on `PATH`.

```bash
git clone https://github.com/aglabx/humas_hmmer
cd humas_hmmer
cargo build --release
# binary at ./target/release/humas-hmmer
```

HMM profiles are not bundled (they are ~60 GB combined). Get them from
[HumAS-HMMER_for_AnVIL](https://github.com/fedorrik/HumAS-HMMER_for_AnVIL):
`AS-HORs-hmmer3.3.2-120124.hmm` (HOR, 1963 models) and `AS-SFs-hmmer3.0.290621.hmm` (SF, 39 models).

## Usage

```bash
humas-hmmer \
  -i /path/to/fasta_dir \
  --hor-hmm AS-HORs-hmmer3.3.2-120124.hmm \
  --sf-hmm  AS-SFs-hmmer3.0.290621.hmm \
  -o /path/to/output \
  -t 48
```

| Flag | Default | Meaning |
|---|---|---|
| `-i`, `--input-dir` | — | Directory scanned for `*.fa` files (one FASTA per job) |
| `--hor-hmm` | — | HOR HMM profile |
| `--sf-hmm` | — | SF HMM profile |
| `-o`, `--output-dir` | `.` | Where BED files are written (created if missing) |
| `-t`, `--threads` | `48` | Total CPU budget; split into `t/4` jobs × 4 `nhmmer` threads |
| `--score-threshold` | `0.7` | Minimum score-to-length ratio for a hit to be kept |

Set `RUST_LOG=debug` for per-file detail; the default level is `info`.

### Output

Four BED9 files per input FASTA, named `{track}-vs-{fasta_stem}.bed`:

| Track | Contents |
|---|---|
| `AS-HOR+SF` | All hits from the HOR profile (HOR models + SF monomers) |
| `AS-HOR` | `AS-HOR+SF` minus SF monomers (model names whose base is ≤ 2 characters) |
| `AS-SF` | All hits from the SF profile |
| `AS-strand` | `AS-SF` recolored by strand (`+` → blue `0,0,255`, `−` → red `255,0,0`) |

Column 9 is an RGB triple from the 218-entry HOR color table, so the files load directly as UCSC
Genome Browser custom tracks.

## How it works

```
FASTA files ─────────────────────────┐
                                     ├──▶ rayon pool (t/4 jobs, 4 nhmmer threads each)
HMM profiles (HOR 1963, SF 39) ──────┘      job i: nhmmer ──pipe──▶ parse ─▶ filter ─▶ BED
```

Per FASTA, per profile:

1. **Scan** — `nhmmer --tblout /dev/stdout -o /dev/null`, stdout consumed by a `BufReader` and parsed line by line.
2. **Threshold** — keep hits where `score / alignment_length ≥ 0.7` (formatted to one decimal first, matching the AWK `sprintf("%.1f")`).
3. **Color** — first matching pattern from the 218-rule table wins; unmatched names get grey `85,85,85`.
4. **Sort** — by target name from the 4th character onward, then start (equivalent to `sort -k 1.4,1 -k 2,2n`).
5. **Overlap filter**, three stages matching the original:
   - `bedmap_max_element` — cluster intervals with ≥10 % overlap of either element, keep the max-scoring one (replaces `bedmap --max-element --fraction-either 0.1`)
   - `exact_dedup` — drop byte-identical BED9 lines (replaces `awk '!seen[$0]++'`)
   - `near_dedup` — if start *or* end is within ±10 bp of the previous record, keep the higher score; on a tie keep the longer interval (replaces `overlap_filter.py`)

### Source layout

| File | Responsibility |
|---|---|
| `main.rs` | CLI (clap derive) |
| `config.rs` | Validation, thread allocation, `*.fa` discovery, output naming |
| `pipeline.rs` | Orchestration: discover → parallel map → 4 BED files per genome |
| `nhmmer.rs` | Spawn `nhmmer`, stream and parse stdout |
| `tblout.rs` | `--tblout` line → `HmmHit` → `BedRecord`, threshold filter, coordinate conversion |
| `color.rs` | 218 case-insensitive regex color rules |
| `overlap.rs` | The three filtering stages |
| `bed.rs` | `BedRecord`, BED9 formatting, sorting |
| `error.rs` | `PipelineError` |

## Differences from the original bash pipeline

Deliberate behavioral changes, not accidents:

- **Deterministic color assignment.** The AWK `getColor()` walks the color table with `for (i in cnA)`,
  whose iteration order is unspecified in POSIX awk. When a model name matches more than one pattern the
  original can therefore pick different colors between runs or between awk implementations. This port
  iterates in source order, so the first pattern listed always wins.
- **No `bedops` / `bedtools` / Python dependency.** The `bedmap` step is a sweep-line in Rust.
- **`--fraction-either` semantics** are implemented as *overlap ≥ frac × length of either element* (OR,
  not AND), matching `bedmap --fraction-either`.
- **Alignment coordinates only.** The AWK script has a `coords=env` mode; this port always uses `ali`
  coordinates, which is what both `hmmer-run.sh` and `hmmer-run_SF.sh` invoke.
- **Input glob is `*.fa`**, not `*.fasta` — matching `HumAS-HMMER_for_AnVIL`.

## Validation status

| Check | State |
|---|---|
| Unit tests (15) | ✅ `cargo test` |
| Byte-identical output vs. bash pipeline on T2T-CHM13 | ⬜ not yet run |
| Full VGP corpus run | ⬜ not yet run |

Until the diff against the reference pipeline is done, treat the output as *expected to match* rather
than *known to match*.

## Development

```bash
cargo test        # 15 unit tests, no nhmmer required
cargo clippy --all-targets -- -D warnings
cargo fmt --check
```

Tests cover BED9 formatting and sorting, tblout parsing and ±strand coordinate conversion, the score
threshold, color lookup (known/case-insensitive/default/FEDOR additions), and all three overlap stages.
They do not require `nhmmer` or HMM profiles.

## Attribution and citation

This is a port, not original science. The pipeline logic, the score/length threshold, and the entire HOR
color table come from **HumAS-HMMER** by Lev I. Uralsky (Institute of Molecular Genetics, Moscow) and the
`HumAS-HMMER_for_AnVIL` modifications by Fedor Ryabov. See [NOTICE.md](NOTICE.md) for the full
provenance chain and licensing situation.

If you use this tool, cite the original work:

> Uralsky LI, Shepelev VA, Alexandrov AA, Yurov YB, Rogaev EI, Alexandrov IA (2019).
> Classification and monomer-by-monomer annotation dataset of suprachromosomal family 1 alpha satellite
> higher-order repeats in hg38 human genome assembly. *Data in Brief* 24:103708.
> [doi:10.1016/j.dib.2019.103708](https://doi.org/10.1016/j.dib.2019.103708)

and HMMER:

> Eddy SR (2011). Accelerated profile HMM searches. *PLoS Comput Biol* 7(10):e1002195.

## License

The Rust code in this repository is MIT licensed — see [LICENSE](LICENSE).

The color table transplanted into `src/color.rs` is derived from upstream work that carries **no license
file**; see [NOTICE.md](NOTICE.md) before redistributing or using it commercially.
