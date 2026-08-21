# Provenance and licensing notice

## What is original here

The Rust source in `src/` is an independent implementation written for this repository and is offered
under the MIT license (see `LICENSE`).

## What is derived

`src/color.rs` contains a 218-entry table mapping alpha-satellite HOR/SF model names to RGB colors. That
table — and the `score / alignment_length ≥ 0.7` threshold, the ±10 bp near-deduplication rule, and the
overall four-track output contract — were transplanted from upstream work, not devised here.

Provenance chain:

| Component | Origin |
|---|---|
| `hmmertblout2bed.awk` (color table, threshold, BED9 contract) | Lev I. Uralsky, Institute of Molecular Genetics, Moscow — <https://github.com/enigene/hmmertblout2bed> |
| `HumAS-HMMER` (pipeline design) | Lev I. Uralsky / enigene — <https://github.com/enigene/HumAS-HMMER> |
| `HumAS-HMMER_for_AnVIL` (`hmmer-run.sh`, `hmmer-run_SF.sh`, `overlap_filter.py`, the "FEDOR" color additions: `S3C1H2-F`, `S3C1H2-E`, `S2C2H2-C`, `S3C11H2-A/B`, `S4C15H4`) | Fedor Ryabov — <https://github.com/fedorrik/HumAS-HMMER_for_AnVIL> |
| HMM profiles (`AS-HORs-*.hmm`, `AS-SFs-*.hmm`) | Not redistributed here. Obtain from `HumAS-HMMER_for_AnVIL`. |

## Licensing situation — read this before redistributing

**None of the upstream repositories carries a license file.** Verified 2026-08-21:

- `enigene/HumAS-HMMER` — no license
- `enigene/hmmertblout2bed` — no license
- `fedorrik/HumAS-HMMER_for_AnVIL` — no license
- `kmiga/alphaAnnotation` (which vendors `HumAS-HMMER_for_AnVIL`; Zenodo [10.5281/zenodo.5715445](https://doi.org/10.5281/zenodo.5715445)) — no license file; the Zenodo deposit is tagged only "Other (Open)"

The associated publication —

> Uralsky LI, Shepelev VA, Alexandrov AA, Yurov YB, Rogaev EI, Alexandrov IA (2019). *Data in Brief*
> 24:103708, [doi:10.1016/j.dib.2019.103708](https://doi.org/10.1016/j.dib.2019.103708)

— is open access under **CC BY-NC-ND 4.0**. That license covers the article and its data; it is not a
software license, and taken literally its `NC` and `ND` terms would preclude both commercial use and
derivative works such as this port. The paper points at `github.com/enigene/hmmertblout2bed` for the
color-coding script but states no terms for it.

Practical consequences:

1. The MIT license in `LICENSE` covers **this repository's Rust code only**. It cannot and does not
   relicense the upstream color table.
2. Redistribution of the color table has been widespread in the field (the T2T/AnVIL toolchain, Zenodo
   deposits, downstream Snakemake wrappers) and no upstream author has objected, but that is convention,
   not permission.
3. If you need clear rights — commercial use in particular — contact the upstream authors.

Corrections and clarifications from the upstream authors are welcome; open an issue and this notice will
be updated.
