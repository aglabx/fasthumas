# humas_hmmer

A single-binary Rust port of the [HumAS-HMMER](https://github.com/enigene/HumAS-HMMER) alpha-satellite
annotation pipeline. It scans a FASTA file with `nhmmer` against one alpha-satellite HMM profile — HOR
or SF — and writes a colour-coded BED9 track.

Everything except `nhmmer` itself runs in memory: no temp files for the annotation, no
`awk`/`sort`/`bedmap`/`python` round-trips.

> **Status: pre-release (v0.1.0).** Output on the plus strand has been diffed against the original bash
> pipeline and matches record for record. On the minus strand it deliberately does **not** match — the
> original coordinates there are wrong, and this port fixes them. It also closes the small gaps HMMER's
> alignment trimming leaves between adjacent monomers. See
> [Differences from the original bash pipeline](#differences-from-the-original-bash-pipeline).

## Why

The original pipeline is a chain of shell steps:

```
nhmmer → hmmertblout2bed.awk → sort → bedmap → awk → overlap_filter.py
```

Each step writes a temp file. It also expects you to chop your assembly into one-sequence `*.fa` files
and hand it a directory of them. This port keeps `nhmmer` (optimized C with SIMD — there is nothing to
gain by rewriting it), replaces everything around it with in-process Rust that consumes nhmmer's
`--tblout` stream directly, and takes one FASTA of any number of sequences.

Inside, the input is divided into one file per sequence and those are annotated concurrently, which is
what keeps the cores busy: `nhmmer` parallelises over blocks of the target database, so handing one
process a single sequence and 48 threads leaves most of them idle. On 34 CDR regions (2.1 Mb) a single
whole-file `nhmmer` call took 62 minutes against 26 for the split run.

## Install

Requires a Rust toolchain and **HMMER ≥ 3.4** on `PATH` — both `nhmmer` and `hmmemit` are used.

3.4 is not optional. nhmmer guesses the alphabet of its target from the first few thousand residues and
a low-complexity start defeats that guess: the human telomere, `TAACCCTAACCCT...`, uses only T, A and C,
each equally valid as an amino acid, so nhmmer refuses the file outright. In T2T-CHM13 that hits chr22
and chrY, whose leading telomere runs past 4 kb. This port passes `--dna` to settle the alphabet, and
that flag only reaches the target from 3.4 on; under 3.3.2 it is accepted and ignored, and those two
chromosomes still fail with `Invalid alphabet type in target for nhmmer`.

```bash
git clone https://github.com/aglabx/humas_hmmer
cd humas_hmmer
cargo build --release
# binary at ./target/release/humas-hmmer
```

HMM profiles are not bundled. Get them from
[HumAS-HMMER_for_AnVIL](https://github.com/fedorrik/HumAS-HMMER_for_AnVIL):
`AS-HORs-hmmer3.3.2-120124.hmm` (HOR) and `AS-SFs-hmmer3.0.290621.hmm` (SF, 39 models).

## Usage

One profile per run, one BED out.

```bash
humas-hmmer -i chm13.fa --hmm AS-HORs-hmmer3.3.2-120124.hmm -o results/chm13_hor
humas-hmmer -i chm13.fa --hmm AS-SFs-hmmer3.0.290621.hmm    -o results/chm13_sf
```

```
results/chm13_hor.bed
results/chm13_sf.bed
```

| Flag | Default | Meaning |
|---|---|---|
| `-i`, `--input` | — | One FASTA file, with any number of sequences |
| `--hmm` | — | The HMM profile to scan with |
| `-o`, `--output` | input file stem | Output path; `.bed` is appended unless already there. An existing directory means "put it in here" |
| `-t`, `--threads` | `48` | Total CPU budget across all concurrent jobs |
| `--threads-per-job` | `4` | Sets how many jobs run at once: `threads / this` |
| `--temp-dir` | output directory | Where the scratch directory is created |
| `--keep-temp` | off | Keep the scratch directory instead of deleting it |
| `--score-threshold` | `0.7` | Minimum score-to-length ratio for a hit to be kept |

Set `RUST_LOG=debug` for per-sequence detail; the default level is `info`.

`-t` is a budget, not a per-process setting. Threads are handed to each sequence **in proportion to its
length**, because sequences are wildly uneven: among chr1's 13 alpha-satellite fields one 4.5 Mb array
holds 86 % of the residues, and an equal share leaves the whole run waiting on that one starved job. Two
bounds keep the split honest — the share is normalised over the longest `threads / threads-per-job`
sequences, the heaviest set that can run at once, and no sequence is given more threads than it has work
for, since nhmmer reads a target in ¼ Mb blocks and hands one block to a worker.

On those 13 fields at `-t 96` that is worth roughly a factor of two: 70 s against 149 s for an equal
share, and better than 79 s for the best hand-tuned `--threads-per-job` we could find. Output is
identical either way; only the scheduling changes.

`--threads-per-job` now only sets how many jobs run at once. Lower it for more concurrency, raise it to
cut peak memory — each running job holds one sequence's records, so memory scales with the number of
jobs, not with the assembly.

Every job is given `-Z` set to the *whole* input's size, so E-values, and with them which hits clear
nhmmer's own reporting threshold, do not depend on how the input was divided. Splitting is an
implementation detail, not a change of results.

Scratch files live in a private `humas-hmmer-tmp-<pid>` directory. Only that directory is ever deleted,
and only if the run succeeds; a failed run leaves it behind and says where.

One caveat for assemblies made of many small contigs: each job pays the cost of loading the HMM profile,
about 30 s for the HOR profile. Thousands of short sequences will spend a visible share of the run
loading the same profile over and over.

### What to feed it

Give it the alpha-satellite regions, **all of them in one FASTA**, one record per region. Not a
directory of files, which is what the original pipeline wants, and not the whole genome.

Cutting to alpha satellite saves the scan from grinding through sequence that cannot contain any: on
chr1, alpha satellite is 5.2 Mb of 248 Mb. Keeping the regions together in one file is what lets the
tool run them concurrently — hand it a single whole chromosome and there is one sequence, so one job,
and `nhmmer`'s own threading over a single target scales poorly (96 threads bought 5.5× over 4 on
chr1, not 24×).

The regions can be cut straight from a CenSat annotation — the `hor`, `dhor` and `mon` classes are the
alpha satellite; `ct` is chromosome arm, and `hsat`/`bsat`/`gsat` are other satellite families.

Whole chromosomes do work, and with `--dna` they no longer fail on telomeric starts, but two things are
worth knowing. Annotation inside the regions is the same either way — re-running chr1 both ways, 30 457
of ~30 550 records are identical down to the model name and score. What differs sits at the edges: a
whole chromosome also finds alpha satellite outside the annotated regions (86 records on chr1), while a
cut region loses the flanking context its first and last monomers would have had, which moves a handful
of records at each boundary. If those edges matter, cut the regions with a couple of hundred bases of
margin and trim the output back afterwards.

### Output

One BED9 file. Sequences appear in the order they occur in the input FASTA; within a sequence, records
are sorted by start coordinate. Column 9 is an RGB triple from the 218-entry HOR color table, so the
file loads directly as a UCSC Genome Browser custom track.

The two derived tracks the original pipeline shipped are one-liners over this file:

```bash
awk -F'\t' 'length($4)!=2' out.bed > AS-HOR.bed                                  # drop SF monomers
awk -F'\t' 'BEGIN{OFS=FS} {$9=($6=="+")?"0,0,255":"255,0,0"; print}' out.bed \
    > AS-strand.bed                                                              # colour by strand
```

## How it works

```
FASTA ─▶ split ─┬─▶ seq 1 ─▶ nhmmer ─▶ filter ─┐
                ├─▶ seq 2 ─▶ nhmmer ─▶ filter ─┼─▶ junction pass ─▶ BED
                └─▶ ...   (jobs run concurrently)
```

Per sequence:

1. **Scan** — `nhmmer --tblout /dev/stdout -o /dev/null`, stdout consumed by a `BufReader` and parsed line by line.
2. **Threshold** — keep hits where `score / alignment_length ≥ 0.7` (formatted to one decimal first, matching the AWK `sprintf("%.1f")`).
3. **Color** — first matching pattern from the 218-rule table wins; unmatched names get grey `85,85,85`.
4. **Sort** — by start.
5. **Overlap filter**, three stages matching the original:
   - `bedmap_max_element` — cluster intervals with ≥10 % overlap of either element, keep the max-scoring one (replaces `bedmap --max-element --fraction-either 0.1`)
   - `exact_dedup` — drop byte-identical BED9 lines (replaces `awk '!seen[$0]++'`)
   - `near_dedup` — if start *or* end is within ±10 bp of the previous record on the same sequence, keep the higher score; on a tie keep the longer interval (replaces `overlap_filter.py`)
6. **Junction pass** — see below. It runs once the whole input has been scanned, because deciding a
   junction needs counts gathered across every sequence.

### Source layout

| File | Responsibility |
|---|---|
| `main.rs` | CLI (clap derive) |
| `config.rs` | Validation, output path, CPU budget split |
| `fasta.rs` | Streaming split into one FASTA per sequence, interval extraction |
| `pipeline.rs` | Orchestration: split → concurrent jobs → junction pass → BED |
| `nhmmer.rs` | Spawn `nhmmer`, stream and parse stdout |
| `tblout.rs` | `--tblout` line → `HmmHit` → `BedRecord`, threshold filter, coordinate conversion |
| `consensus.rs` | Model consensus sequences via `hmmemit -c` |
| `junction.rs` | Giving back MEA-trimmed bases at tandem junctions |
| `color.rs` | 218 case-insensitive regex color rules |
| `overlap.rs` | The three filtering stages |
| `bed.rs` | `BedRecord`, BED9 formatting, sorting |
| `error.rs` | `PipelineError` |

## Differences from the original bash pipeline

Deliberate behavioral changes, not accidents.

### Minus-strand coordinates are corrected

`nhmmer` reports alignment coordinates in query order, so on the minus strand `alifrom > alito` —
column 7 is the **end** of the hit on the target, not its start. `hmmertblout2bed.awk` computes
`aliFrom = $7 - 1; aliTo = $8` and then prints `aliTo, aliFrom` for the minus strand, which puts the
`-1` on the end instead of the start. Every minus-strand interval therefore comes out shifted one base
to the right and **two bases too short**, and the same wrong length is what feeds the score/length
threshold.

The signature is plain in any original output with minus-strand calls. Over 34 annotated CDR regions
(12 568 records, 4 255 on the minus strand), tandem monomers butt-joined on the plus strand meet
exactly, while the *same* junctions on the minus strand are all 2 bp apart — the plus-strand
distribution shifted by two. Median monomer length came out 171 bp on the plus strand and 169 on the
minus.

This port converts `[min(alifrom, alito) - 1, max(alifrom, alito))` on both strands, so a hit and its
reverse-complement twin get identical coordinates. Re-running those regions, all 4 255 minus-strand
records moved one base left at the start and one base right at the end, and nothing else about them
changed.

### Gaps left by MEA trimming are closed

HMMER builds its alignment by maximum expected accuracy and drops residues whose posterior probability
is low, however good their emission score. At a tandem junction the neighbour is another copy of the
same monomer, so the register is genuinely ambiguous and a few terminal residues get trimmed. What is
left is a gap of one to three bases between monomers that in the sequence are flush against each other.

Those gaps need not be guessed at. The tblout carries the model coordinates of every alignment
(`hmmfrom`, `hmmto`), so a hit that stops short of the model's first or last position says exactly how
many positions were trimmed, and at which end. Over 308 junctions on a 52 kb CENP-B region, 99.4 % of
gaps equal `trimmed(left) + trimmed(right)` exactly.

The pass hands those bases back when either test passes:

- **the model expects them** — the bases match that model's consensus (`hmmemit -c`) at the trimmed
  positions, walked backwards and complemented for minus-strand hits, at ≥ 50 % identity;
- **the array vouches for them** — the same base has turned up at the same model position at least 10
  times across the whole run, holding at least half the observations there. A base seen 147 times at
  position 1 of `J5` is a variant of that monomer, not noise, even where the model's consensus says
  otherwise.

The two sides of a junction are judged separately, so a junction may close from one side only. When the
trimmed positions add up to *more* than the gap, the genome is short of bases the model expects — a
deletion — and the gap is still closed with the bases that are there, provided only one side is missing
anything. Since the two sides between them own exactly the gap and nothing more, closing either or both
can never produce an overlap.

None of this is tunable. Monomer boundaries are a property of the models and the sequence, not a knob
for the caller.

One profile-side wrinkle: the HOR profile carries 1923 models under only 1238 distinct names, 58 of
which cover models of *different lengths*, and a tblout hit reports its model by name alone. One end of
every hit is free of the model length (`hmm_from - 1` bases were dropped there whatever the model), so
the pass tries each candidate length for the other end and acts only when exactly one combination adds
up. On the regions above that left 1 junction of 1516 unresolved.

### The rest

- **One profile in, one BED out.** The original globs a directory for `*.fa` and writes four files per
  input; this takes a single file and writes one. `AS-HOR` and `AS-strand` are `awk` one-liners over it.
- **Deterministic color assignment.** The AWK `getColor()` walks the color table with `for (i in cnA)`,
  whose iteration order is unspecified in POSIX awk, so a name matching more than one pattern can get
  different colors between runs. This port iterates in source order: the first pattern listed wins.
- **Sorting is by whole sequence name.** The original used `sort -k 1.4,1 -k 2,2n`, comparing names from
  the 4th character on — safe only because each input file held one chromosome.
- **`near_dedup` resets at each sequence boundary.** `overlap_filter.py` carried its `prev_*` state
  across the whole file; on multi-sequence input that lets the last record of one sequence suppress the
  first record of the next.
- **The ±10 bp window in `near_dedup` is symmetric.** `overlap_filter.py` tests
  `start in range(prev_start-10, prev_start+10)`, a window of `[-10, +9]`. This port uses
  `abs(delta) < 10`. The two disagree only when a coordinate sits exactly 10 bp below the previous one.
- **Ties inside a `bedmap` cluster are broken differently.** When two models score identically over the
  same interval, `bedmap --max-element` and this port can keep different ones. Over 12 568 records that
  happened once, and only the model name differed.
- **No `bedops` / `bedtools` / Python dependency.** The `bedmap` step is a sweep-line in Rust.
- **`--fraction-either` semantics** are implemented as *overlap ≥ frac × length of either element* (OR,
  not AND), matching `bedmap --fraction-either`.
- **Alignment coordinates only.** The AWK script has a `coords=env` mode; this port always uses `ali`
  coordinates, which is what both `hmmer-run.sh` and `hmmer-run_SF.sh` invoke.
- **Score column keeps one decimal.** The AWK emits `sprintf("%.1f", score)`, but `bedmap` in the
  original chain rewrites the column to `%f` (`154.400000`). This port keeps `154.4`.

## Validation status

Measured on 34 CDR regions of one haplotype assembly (2.1 Mb, 12 568 records), against the original
pipeline run on the same input with the same profile.

| Check | State |
|---|---|
| Unit tests (41) | ✅ `cargo test` |
| Plus strand vs. bash pipeline — SF profile, 2.1 Mb CENP-B region, 12 518 records | ✅ identical record for record |
| Plus strand vs. original pipeline — HOR profile, 8 313 records | ✅ identical but for 1 record, where two equally-scoring models tie |
| Minus strand vs. original pipeline — 4 255 records | ⛔ differs by design: one base wider at each end, the original being 2 bp short |
| Split, concurrent run vs. one whole-file `nhmmer` call | ✅ identical record for record |
| Multi-sequence input — 34 sequences in one FASTA | ✅ 12 568 records, matching the original's 34 single-sequence runs |
| Full VGP corpus run | ⬜ not yet run |

Junctions meeting flush, same input, same profile:

| | original | this port |
|---|---|---|
| HOR profile | 7 286 / 12 534 (58.1 %) | **12 350 / 12 534 (98.5 %)** |
| SF profile | 4 428 / 12 531 (35.3 %) | **10 975 / 12 531 (87.6 %)** |

Monomers of exactly 171 bp went from 4 089 to 7 332 of 12 568, and the mean length gap between strands
(169.37 plus vs 167.75 minus) closed to 169.52 vs 169.96.

What is left open on the HOR profile is 218 junctions of 12 534: 91 gaps wider than 3 bp, which look
like unannotated sequence rather than trimming, 79 where the bases matched neither the consensus nor a
recurring variant, and 14 overlaps that come from the `bedmap` stage and are present in the original
output too. The SF profile leaves more (1 090 one-base gaps) because its 39 family-level models cover
monomers from many arrays at once, so their terminal positions are genuinely polymorphic — at `D1`
position 171, for instance, `A` appears 89 times and `C` 86 times out of 642.

## Development

```bash
cargo test        # 41 unit tests, no nhmmer required
cargo clippy --all-targets -- -D warnings
cargo fmt --check
```

Tests cover BED9 formatting and sorting, tblout parsing and ±strand coordinate conversion, the score
threshold, output path construction, the proportional thread plan, FASTA splitting and interval extraction,
`hmmemit` output parsing, the junction pass on both strands including recurrence and deletions, color
lookup, and all three overlap stages.

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

The MIT grant covers this repository's Rust code. It does **not** extend to the alpha-satellite color
table transplanted into `src/color.rs`, which is derived from upstream work that carries no license file
(the header of that file says so too). See [NOTICE.md](NOTICE.md) before redistributing it or using it
commercially.
