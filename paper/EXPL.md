# FastHumAS: Extended Technical & Biological Explanation (EXPL)
## Architectural Paradigm, Concordance with Original HumAS-HMMER, and Transferability of Results

**Authors:** Marina Popova, Simona Giunta, Aleksey Komissarov  
**Reference Codebase:** [https://github.com/aglabx/fasthumas](https://github.com/aglabx/fasthumas)  
**Evaluated Assembly:** T2T-CHM13 v2.0 (`GCF_009914755.1`, 3,117,275,501 bp)  
**Profile Libraries:**
- `AS-HORs-hmmer3.3.2-120124.hmm` (1,963 monomer models, MD5: `4c642b44df914a47bed2cd3d0998f8e6`)
- `AS-SFs-hmmer3.0.290621.hmm` (39 SF models, MD5: `b7c211e2e72edfbe9f9068a5432f4025`)

---

## 1. Executive Summary & Conceptual Framework: "A Minimap2 for HMMER"

The central methodological question raised regarding FastHumAS is:
> *"FastHumAS replaces exhaustive Profile HMM alignment with a two-stage heuristic (consensus Master Profile peak detection followed by bounded banded dynamic programming). Does this speedup (23.4× wall-clock acceleration, zero intermediate scratch disk churn) come at the cost of annotation accuracy? Does it miss genuine alpha-satellite monomers, and are its outputs 100% transferable to the existing T2T CenSat ecosystem?"*

The conceptual analogy to **minimap2** (Li, 2018) is direct and instructive:
- When `minimap2` replaced quadratic `BLAST`/`BWA-MEM` by combining minimizer seeding with banded dynamic programming, it reduced alignment runtimes by orders of magnitude while preserving $>99.5\%$ mapping sensitivity and improving alignment quality across indels.
- In exactly the same way, **FastHumAS acts as an in-memory heuristic aligner for Profile HMM libraries**. Instead of quadratic, unconstrained dynamic programming scanning every genomic base against 2,000 separate Profile HMM models (the computational bottleneck of `nhmmer`), FastHumAS seeds candidate monomer endpoints using a unified 171-bp consensus Master Profile filter (Stage 1), and then performs bounded profile-anchored dynamic programming in an asymmetric window $[-14, +28]$ bp (Stage 2).

As proven empirically below across **200,337 loci across 10 finished chromosomes of the T2T-CHM13 genome**, FastHumAS achieves **99.44% sensitivity, 99.96% precision, 100.00% strand concordance, and 99.32% label concordance** against the original HumAS-HMMER pipeline, while filtering out degenerate sub-monomer noise and fixing systematic coordinate offset defects in the legacy shell pipeline.

---

## 2. The Original Baseline: What Was Evaluated?

The baseline used for all comparative evaluations is the **unmodified, official HumAS-HMMER toolchain**:
1. **Official Pipeline Repository:** [`fedorrik/HumAS-HMMER_for_AnVIL`](https://github.com/fedorrik/HumAS-HMMER_for_AnVIL) (authored by Fedor Ryabov, Nicolas Altemose, and Karen Miga).
2. **Execution Environment:** Executed on cluster `aglab0` using the official shell scripts:
   - `hmmer-run.sh` (executes `nhmmer --dna --cpu 48 --notextw --noali --tblout ...` against `AS-HORs-hmmer3.3.2-120124.hmm`)
   - `hmmer-run_SF.sh` (executes `nhmmer` against `AS-SFs-hmmer3.0.290621.hmm`)
   - `hmmertblout2bed.awk` (AWK-based tabular parsing and score threshold filtering at threshold 0.70)
   - `overlap_filter.py` (greedy non-overlapping interval selection prioritizing HOR over SF)
3. **Independent Third-Party Verification:** We additionally cross-referenced the independent published CHM13 alpha-satellite annotations from **Gao et al. 2025 (*Science*)** (71,282 active HOR units across all 23 nuclear chromosomes in `/mnt/data/t2t/gao2025_centromeres/annotations/hgsvc3_chm/`), confirming model and coordinate concordance.

---

## 3. Empirical Accuracy & Concordance Metrics

Across 10 finished nuclear chromosomes of T2T-CHM13 v2.0 (chr1, chr3, chr8, chr10, chr11, chr12, chr14, chr17, chr22, chrY), encompassing **200,337 legacy loci and 199,302 FastHumAS loci**, direct locus-by-locus benchmarking with strict reciprocal overlap ($\ge 50\%$ on both intervals) and 1:1 bipartite matching establishes:

| Metric | Measured Value | Biological Meaning & Assessment |
|---|:---:|---|
| **Sensitivity (Recall)** | **99.44%** | FastHumAS captures **199,220 of the 200,337** legacy loci, demonstrating high recovery across chromosomes. |
| **Precision** | **99.96%** | 199,220 of 199,302 FastHumAS calls are validated by legacy annotations (virtually zero spurious false positives). |
| **Strand Concordance** | **100.00%** | 199,220 / 199,220 agreement across all matched intervals (zero polarity inversion errors). |
| **Subfamily / HOR Label Agreement** | **99.32%** | Exact model name assignment (e.g. `S1C10H1L.5`, `Ka`, `Ba`) matches character-for-character in 197,874 / 199,220 pairs. |
| **Boundary Concordance ($\le 3$ bp)** | **97.59%** | 194,419 / 199,220 monomers match within 1–3 bp; boundary shifts reflect intentional heuristic junction healing. |
| **Boundary Concordance ($\le 10$ bp)** | **99.74%** | 198,711 / 199,220 monomers match within 10 bp, confirming positional stability. |
| **Exact 0-bp Boundary Match** | **41.11%** | 81,909 / 199,220 monomers match bit-for-bit to the exact single nucleotide. |

---

## 4. Forensic Analysis of Discordant Loci

A crucial test for any heuristic aligner is understanding what occurs at discordant loci:

### 4.1 Legacy-Unique Loci ($n = 1,117$, 0.56%)
Why did `nhmmer` detect 1,117 loci that FastHumAS did not output?
- **Length Distribution:** Mean length is only **105.2 bp** (median: 104.0 bp, range: 24–204 bp).
- **Score Distribution:** Mean score is **87.8**; 65.4% scored $<100$, and 94.1% scored $<150$.
- **Biological Nature:** These are **truncated sub-monomer fragments** located predominantly in degenerate pericentromeric transitions and transposable element borders. `nhmmer`'s local alignment picked up partial 30–80 bp motif matches near the 0.70 score density threshold that fail the full-monomer criteria of FastHumAS. FastHumAS omits these sub-monomers, preventing artificial fragmentation.

### 4.2 FastHumAS-Unique Loci ($n = 82$, 0.04%)
Why did FastHumAS detect 82 loci not present in the legacy output?
- **Length Distribution:** Mean length is **168.1 bp** (median: 169.0 bp), matching canonical 171-bp monomer architecture.
- **Score Distribution:** Mean score is **114.2** (range: 99.3–121.7).
- **Biological Hypothesis:** Length and coordinate properties support the hypothesis that these loci represent **canonical monomer candidates** that were discarded by the legacy pipeline during downstream `overlap_filter.py` resolution due to minus-strand coordinate collision artifacts. Definitive biological validation of discordant loci remains an open hypothesis subject to orthogonal confirmation.

---

## 5. Elimination of Systematic Legacy Defects

FastHumAS not only matches legacy sensitivity, but fixes two severe defects inherent in the legacy shell implementation:

### 5.1 The Minus-Strand 2-bp Offset Bug (~166,660 Gaps)
In `hmmertblout2bed.awk`:
```awk
if ($9 == "-") {
    aliFrom = $7 - 1;  # BUG: on minus strand, $7 is the larger coordinate!
    aliTo = $8;
}
```
Because `$7 > $8` on the minus strand, subtracting 1 from `$7` placed the adjustment on the wrong end and left `$8` unadjusted. This systematically shifted every minus-strand annotation 1 bp to the right and shortened it by 2 bp.
- **Empirical Quantification:** On the 10 evaluated chromosomes, legacy HumAS-HMMER contained 100,900 micro-gaps (1–3 bp). **68,292 of these (67.68%) were exact 2-bp gaps between adjacent minus-strand monomers**.
- **Genome-wide Extrapolation:** Extrapolating to the full CHM13 assembly ($68,292 \times \frac{488,904}{200,337}$), this bug created **~166,660 artificial minus-strand gaps**.
- **FastHumAS Fix:** FastHumAS uses symmetric interval projection:
  $$\text{start} = \min(\text{aliFrom}, \text{aliTo}) - 1, \quad \text{end} = \max(\text{aliFrom}, \text{aliTo})$$
  This eliminates 98.05% of artificial 2-bp minus-strand gaps (decreasing from 68,292 to 1,335).

### 5.2 Junction Healing
Local alignment drops terminal bases due to register ambiguity at monomer junctions. FastHumAS tests 1–3 bp gaps against model consensus and array recurrence:
- **Whole-Genome Acceptance:** FastHumAS evaluates candidate gaps in HOR and SF, accepting 35.1% and 55.7% of boundary extensions respectively. Gaps $>3$ bp and divergent transitions are strictly left open.
- **Flush Continuity:** The proportion of flush abutting junctions (gap = 0 bp) increases from **40.41% (80,955) to 93.72% (186,770)** across evaluated chromosomes. Micro-gaps (1–3 bp) decrease by 91.92% (from 100,900 to 8,152).
- **Active Dense Arrays:** In the canonical active centromeric core of chromosome 22 (`NC_060946.1:12,788,180–15,711,065`, array `hor_22_9` / `S2C14/22H1L`, spanning 2.92 Mb with 17,145 monomers), **99.59% of junctions (17,073 / 17,144)** are resolved completely flush (gap = 0 bp), compared to **0.00% (0 / 17,150)** in legacy HumAS-HMMER. Only 65 micro-gaps, 5 gaps > 3 bp, and 1 overlap (-1 bp) remain. Investigation of the 5 gaps > 3 bp reveals:
  (i) a 219-bp unannotated gap at a degenerate pericentromeric boundary (`12,789,031–12,789,250`) where legacy called two sub-100-bit fragments (`S2C9H1L.7` [score 97.9, 118 bp] and `S2C14/22H1L.3` [score 90.4, 95 bp]);
  (ii) three identical recurrent 31-bp unannotated sequences (`CAACAAAAAGTGTTTTTCAAAACTGCTGTAT`) between monomers `.7` and `.6`, where FastHumAS calls full-length 177-bp monomers while legacy split monomer `.6` into two fragments; and
  (iii) an unannotated 6-bp sequence (`CTAAAA`, `15,709,700–15,709,706`) in the assembly sequence (where legacy similarly leaves an 8-bp unannotated gap `15,709,699–15,709,707`).

---

## 6. Scoring Model Calibration (The 0.85 Scaling Factor)

Stage 2 dynamic programming performs a profile-anchored semi-global alignment. To map Gotoh match log-odds sums to HMMER's Null2-adjusted bit score scale, FastHumAS applies a linear scaling factor $0.85$.

We calibrated and validated this factor across **81,661 exact-boundary, model-matched monomers** across 10 chromosomes with a strict training/validation partition:
- **Training Split (chr1, chr3, chr8; $n = 23,506$):**
  - Mean ratio: $0.8495 \pm 0.0074$ (Sample SD)
  - SEM: $0.000048$
  - Median ratio: $0.8504$ (IQR: 0.8469–0.8534)
  - 95% empirical range: $[0.8305, 0.8599]$
- **Independent Validation Split (chr10, 11, 12, 14, 17, 22, Y; $n = 58,155$):**
  - Mean ratio: $0.8505 \pm 0.0078$ (Sample SD)
  - SEM: $0.000032$
  - Median ratio: $0.8522$ (IQR: 0.8478–0.8555)
  - 95% empirical range: $[0.8299, 0.8595]$

### Score Stratification & Null2 Dynamics
Analyzing the calibration ratio across score tiers reveals the expected thermodynamic effect of HMMER's Null2 background model:
- **High-Scoring Monomers ($\ge 150$ bits; $n = 80,341$, 98.4% of evaluated exact-boundary calibration pairs):**
  - Mean ratio: $0.8507 \pm 0.0065$ (Median: $0.8517$)
- **Mid-Scoring Monomers (100–150 bits; $n = 1,320$, 1.6% of evaluated exact-boundary calibration pairs):**
  - Mean ratio: $0.8198 \pm 0.0128$ (Median: $0.8210$)
- **Low-Scoring Monomers (< 100 bits; $n = 0$):**
  - Zero exact-boundary model-matched monomers fall below 100 bits.

Weighted mean across splits is $\mathbf{0.850209}$, matching the weighted mean across score tiers ($\mathbf{0.850210}$) within $<10^{-6}$.

For canonical alpha-satellite monomers ($\ge 150$ bits, representing 98.4% of the evaluated exact-boundary calibration pairs), the $0.85$ factor is essentially exact ($0.8507$). For lower-scoring matches near the detection threshold (100–150 bits), composite differences in HMMER's statistical model yield a lower empirical ratio ($0.8198$), meaning that multiplying by $0.85$ produces a bit score approximately ~3.7% higher than HMMER. FastHumAS includes the candidate density check $\text{Score}_{\text{scaled}} / \text{Length} \ge \text{threshold} - 0.05$ as an **upstream inclusion pre-filter** aligned with legacy AWK rounding (`format!("{:.1}", score_per_len)`), ensuring that any candidate that could round to $\ge 0.70$ is not prematurely discarded before final filtering.

---

## 7. Production Overlap Resolution (BEDOPS Parity)

FastHumAS strictly replicates the 3-stage overlap filtering pipeline of legacy HumAS-HMMER while guaranteeing coordinate order preservation for junction healing:
1. **Local Non-Maximum Suppression (`bedmap --max-element --fraction-either 0.1`):** Maps each reference interval to the highest-scoring candidate among all intervals sharing $\ge 10\%$ overlap with either element. When candidate scores are identical, FastHumAS selects candidate elements matching BEDOPS `ScoreThenGenomicCompareGreater`: preferring higher start coordinate, then higher end coordinate, and for identical coordinates preserving the earliest candidate in `sort-bed` order `(chrom, start, end, name, strand)`.
2. **Exact Deduplication (`exact_dedup`):** Eliminates duplicate records, matching legacy `awk '{if(!($0 in a)){a[$0]; print}}'`.
3. **Near-Duplicate Suppression (`near_dedup`):** Directly processes the deduplicated stream, sequentially eliminating records whose start or end coordinates lie within the half-open window $[-10, +10)$ bp (`range(prev - 10, prev + 10)`) of an adjacent higher-scoring survivor.
4. **Coordinate Order Restoration:** Coordinates are sorted prior to downstream junction healing without altering filtering membership, ensuring that `junction::candidates` evaluates all true adjacent monomer junctions along the chromosome.

---

## 8. Gotoh DP Recurrence Fix Validation

A rigorous comparison of whole-genome CHM13 annotations before and after correcting the Gotoh affine gap DP recurrence confirmed that **99.976% of records (488,636 / 488,754) are completely bit-identical**:
- **Discordant Records:** The 125 records affected by the gap fix have a mean length of **169.5 bp** (122/125 $\ge 160$ bp, mean score 164.2 bits). They represent full-length competitive monomer alignments whose optimal DP path shifted under the corrected gap open/extension penalty, rather than truncated fragments.
- **Boundary Refinements:** The 4 boundary refinement loci exhibited length adjustments between 2 and 12 bp (two $-2$ bp, one $-9$ bp, one $-12$ bp).

---

## 9. Transferability & Ecosystem Compatibility

FastHumAS ensures 100% drop-in portability for any existing T2T CenSat workflow:
1. **Standard Multi-Track BEDs:** Produces `AS-HOR+SF.bed`, `AS-HOR.bed`, `AS-SF.bed`, and `AS-strand.bed` with the official 218-color RGB palette.
2. **Downstream Pipeline Support:** Direct input compatibility with `CenStats`, `CenMAP`, `stv_row`, `StringDecomposer`, and the UCSC Genome Browser.
3. **Built-in `--legacy` Verification Mode:** FastHumAS provides a native `--legacy` flag that orchestrates the original `nhmmer` subprocess pipeline directly from Rust, enabling researchers to perform bit-level verification whenever required.

All reproducibility scripts, manifests, and logs are publicly available in the GitHub repository at [https://github.com/aglabx/fasthumas/tree/main/benchmarks](https://github.com/aglabx/fasthumas/tree/main/benchmarks).
