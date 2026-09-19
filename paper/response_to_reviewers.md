# Point-by-Point Response to Reviewer 3 (Iteration 3)

**Manuscript Title:** FastHumAS: Ultra-Fast In-Memory Annotation of Centromeric Alpha-Satellite Higher-Order Repeats in Complete Genomes  
**Authors:** Marina Popova, Simona Giunta, and Aleksey Komissarov  
**Repository:** [https://github.com/aglabx/fasthumas](https://github.com/aglabx/fasthumas)  

We once again express our sincere gratitude to Reviewer 3 for their rigorous and insightful critique. The reviewer's meticulous verification identified a subtle ordering bug between the overlap filter and junction healing, an omission in the calibration script, an erroneous mathematical sign in our narrative regarding Null2, and hardcoded summaries in the healing analysis script.

In this revision, we have addressed every single point raised by the reviewer:

1. **Healing Input Coordinate Order (`src/overlap.rs`, `src/scanner.rs`, `src/pipeline.rs`):** We resolved the coordinate reordering caused by `bedmap_max_element`. FastHumAS now explicitly enforces `sort_bed_records` after `exact_dedup` (before `near_dedup`) and again before returning from `filter_overlaps`. The reviewer's exact counterexample ($A [0, 176), B [5, 93), C [158, 329), D [330, 501)$) has been added as a unit test (`test_filter_overlaps_restores_order_for_healing`); the output is strictly sorted as $A, C, D$, and the 1-bp junction between $C$ and $D$ is properly evaluated and healed.
2. **BEDOPS Tie-Breaking Parity (`src/bed.rs`, `src/overlap.rs`):** We expanded `sort_bed_records` to sort by `(chrom, start, end, name, strand)` matching BEDOPS `sort-bed`. In `bedmap_max_element`, equal-score candidate ties are resolved by selecting the candidate appearing earliest in `sort-bed` order. A dedicated unit test (`test_bedmap_tie_break_name`) confirms that for equal-score intervals with names `B` and `A`, `A` is preserved.
3. **Scoring Calibration Statistical Consistency:** We fixed the omission of chromosome 17 (`NC_060941.1`) in `calibrate_scoring.py`. Across all 81,661 exact-boundary, model-matched monomers across 10 finished chromosomes:
   - **Training set (chr1, 3, 8; $n = 23,506$):** Mean $0.8495 \pm 0.0074$ (SEM: $0.000048$, median $0.8504$).
   - **Independent Validation set (chr10, 11, 12, 14, 17, 22, Y; $n = 58,155$):** Mean $0.8505 \pm 0.0078$ (SEM: $0.000032$, median $0.8522$).
   - **Weighted Mean Consistency:**
     $$\bar r_{\text{train/test}} = \frac{23,506 \cdot 0.849491 + 58,155 \cdot 0.850499}{81,661} = \mathbf{0.850209}$$
     $$\bar r_{\text{score tiers}} = \frac{80,341 \cdot 0.850709 + 1,320 \cdot 0.819812}{81,661} = \mathbf{0.850208}$$
     The discrepancy identified by the reviewer is completely eliminated (difference $< 10^{-6}$).
4. **Correction of the `threshold - 0.05` Narrative:** We fully acknowledge the reviewer's mathematical correction. For lower-scoring matches (100–150 bits), the empirical ratio is ~0.82, so multiplying by 0.85 over-scores relative to HMMER by ~3.7%. Lowering the threshold to `threshold - 0.05` is a relaxation in the same direction, not a compensation. We have revised the text in `paper.md`, `EXPL.md`, and this response: `threshold - 0.05` is an **upstream inclusion pre-filter** designed to match legacy AWK score density rounding (`format!("{:.1}", score_per_len)`), preventing candidates that round up to 0.70 from being prematurely eliminated.
5. **Reproducible chr22 Active Array Analysis & Characterization of Gaps > 3 bp:** We created a dedicated reproducible script ([`benchmarks/analyze_active_hor_chr22.py`](file:///Users/akomissarov/Dropbox/workspace/new/biology/humas_hmmer/benchmarks/analyze_active_hor_chr22.py)) and summary log ([`benchmarks/logs/chr22_active_hor_continuity.summary.txt`](file:///Users/akomissarov/Dropbox/workspace/new/biology/humas_hmmer/benchmarks/logs/chr22_active_hor_continuity.summary.txt)). In the active HOR core `NC_060946.1:12,788,180–15,711,065` (17,145 monomers):
   - FastHumAS achieves **99.59% flush contiguous junctions** (17,073 / 17,144) vs **0.00%** (0 / 17,150) in legacy HumAS-HMMER.
   - Exactly 5 gaps > 3 bp remain: one 219-bp gap at a degenerate boundary (where legacy had fragmented <100-bit sub-monomers), three 31-bp recurrent structural variant insertions between monomers `.7` and `.6` (where legacy fragmented monomer `.6` into two partial hits), and one genuine 6-bp assembly gap confirmed in both pipelines.
   - We explicitly clarify that 99.59% represents the end-to-end continuity of the full annotation pipeline (coordinate rectification, banded Gotoh DP, local NMS, and consensus-anchored healing).
6. **Dynamic 10-Chromosome Healing Analysis:** We completely rewrote [`benchmarks/analyze_junction_healing.py`](file:///Users/akomissarov/Dropbox/workspace/new/biology/humas_hmmer/benchmarks/analyze_junction_healing.py) to eliminate all hardcoded numbers and calculate all metrics dynamically. Across the 10 evaluated chromosomes (200,337 legacy loci):
   - Flush junctions increase from **40.41%** (80,955) in legacy to **93.72%** (186,770) in FastHumAS.
   - Micro-gaps (1–3 bp) decrease by **91.92%** (from 100,900 down to 8,152).
   - 2-bp minus-strand gaps decrease by **98.05%** (from 68,292 down to 1,335).
   - Assembly-wide extrapolation: $68,292 \times (488,904 / 200,337) = \mathbf{166,660}$ artificial minus-strand gaps genome-wide.
7. **Harmonization of Concordance Summary & CLI Ratios:** We updated [`benchmarks/compare_fasthumas_vs_nhmmer.py`](file:///Users/akomissarov/Dropbox/workspace/new/biology/humas_hmmer/benchmarks/compare_fasthumas_vs_nhmmer.py) to print exact numerators and denominators for strand agreement (199,220 / 199,220 = 100.00%), label agreement (197,874 / 199,220 = 99.32%), and boundary agreement $\le 3$ bp (194,419 / 199,220 = 97.59%). All numbers across text, tables, and logs are strictly synchronized.
8. **Inclusion of Review Materials in Repository:** All review documents—[`paper/paper.md`](file:///Users/akomissarov/Dropbox/workspace/new/biology/humas_hmmer/paper/paper.md), [`paper/paper.html`](file:///Users/akomissarov/Dropbox/workspace/new/biology/humas_hmmer/paper/paper.html), [`paper/EXPL.md`](file:///Users/akomissarov/Dropbox/workspace/new/biology/humas_hmmer/paper/EXPL.md), and [`paper/response_to_reviewers.md`](file:///Users/akomissarov/Dropbox/workspace/new/biology/humas_hmmer/paper/response_to_reviewers.md)—are now committed directly inside the repository tree under [`paper/`](file:///Users/akomissarov/Dropbox/workspace/new/biology/humas_hmmer/paper) and linked in `README.md`.

Below is our detailed, point-by-point response.

---

## 1. Restoration of Coordinate Order for Junction Healing & BEDOPS Tie-Breaking Parity

> **Reviewer 3:** *"После исправления overlap-фильтра нарушается требуемый порядок входа для healing... Новый filter_overlaps не восстанавливает координатный порядок, а следующий этап junction::candidates требует сортированный вход. Воспроизводимый пример... [0,176), [5,93), [158,329), [330,501)... Результат: C, A, D... Между C и D есть зазор длиной 1 bp... пропуская C→D... Дополнительная native-проверка обнаружила различие при одинаковых координатах и score, но разных именах: B 200.0, A 200.0... Настоящие sort-bed + bedmap сохраняют A."*

### Response:
We thank the reviewer for identifying both of these critical issues.

### 1.1 Root Cause and Algorithmic Fix for Coordinate Ordering
The reviewer's analysis of the interaction between `bedmap_max_element` and `junction::candidates` is completely correct:
- `bedmap_max_element` emits local maximum elements in the order of the reference query intervals.
- When an earlier reference interval maps to a downstream maximum element (e.g. $[0, 176)$ mapping to $[158, 329)$) while a subsequent reference interval maps to an upstream element (e.g. $[5, 93)$ mapping to $[0, 176)$), the emitted sequence is $[158, 329)$, $[0, 176)$, $[330, 501)$.
- `exact_dedup` preserves the order of first appearance. If records are not re-sorted before `near_dedup` and `junction::candidates`, the sequential pairing evaluates $[158, 329) \to [0, 176)$ and $[0, 176) \to [330, 501)$, completely missing the genuine adjacent junction $[158, 329) \to [330, 501)$ (gap = 1 bp).

**The Fix:** In [`src/overlap.rs`](file:///Users/akomissarov/Dropbox/workspace/new/biology/humas_hmmer/src/overlap.rs), `filter_overlaps` now explicitly enforces coordinate sorting:
1. After `exact_dedup` (before `near_dedup`), ensuring that `near_dedup` only compares genuinely adjacent intervals along the chromosome.
2. After `near_dedup` (before returning), ensuring that `filter_overlaps` always returns records strictly sorted by `(chrom, start, end, name, strand)`.
3. In `write_bed_file`, records are sorted prior to emitting BED output.

We added the reviewer's exact counterexample as a unit test:
```rust
#[test]
fn test_filter_overlaps_restores_order_for_healing() {
    let mut recs = vec![
        make_rec("chr1", 0, 176, "A", 150.0),
        make_rec("chr1", 5, 93, "B", 100.0),
        make_rec("chr1", 158, 329, "C", 200.0),
        make_rec("chr1", 330, 501, "D", 190.0),
    ];
    crate::bed::sort_bed_records(&mut recs);
    let filtered = filter_overlaps(&recs);
    assert_eq!(filtered.len(), 3);
    assert_eq!(filtered[0].name, "A");
    assert_eq!(filtered[0].start, 0);
    assert_eq!(filtered[1].name, "C");
    assert_eq!(filtered[1].start, 158);
    assert_eq!(filtered[2].name, "D");
    assert_eq!(filtered[2].start, 330);
}
```

### 1.2 BEDOPS Tie-Breaking Parity
BEDOPS `sort-bed` sorts by `chrom`, `start`, `end`, and then lexicographically across remaining columns (`name`, `score`, `strand`). When candidate scores are tied, `bedmap --max-element` selects the candidate that appears earliest in the sorted order.

**The Fix:**
1. In [`src/bed.rs`](file:///Users/akomissarov/Dropbox/workspace/new/biology/humas_hmmer/src/bed.rs), `sort_bed_records` now sorts by:
   ```rust
   a.chrom.cmp(&b.chrom)
       .then(a.start.cmp(&b.start))
       .then(a.end.cmp(&b.end))
       .then(a.name.cmp(&b.name))
       .then(a.strand.cmp(&b.strand))
   ```
2. In [`src/overlap.rs`](file:///Users/akomissarov/Dropbox/workspace/new/biology/humas_hmmer/src/overlap.rs), `bedmap_max_element` resolves score ties by selecting the candidate appearing earliest in `sort-bed` order:
   ```rust
   let is_better = cand.score > best.score
       || ((cand.score - best.score).abs() < 1e-6
           && (cand.start < best.start
               || (cand.start == best.start
                   && (cand.end < best.end
                       || (cand.end == best.end
                           && (cand.name < best.name
                               || (cand.name == best.name && cand.strand < best.strand)))))));
   ```
We verified this with a dedicated unit test matching the reviewer's case:
```rust
#[test]
fn test_bedmap_tie_break_name() {
    let mut recs = vec![
        make_rec("chr1", 0, 171, "B", 200.0),
        make_rec("chr1", 0, 171, "A", 200.0),
    ];
    crate::bed::sort_bed_records(&mut recs);
    let filtered = filter_overlaps(&recs);
    assert_eq!(filtered.len(), 1);
    assert_eq!(filtered[0].name, "A");
}
```
All 53 unit tests pass.

---

## 2. Scoring Calibration Arithmetic, Null2 Explanation, and Validation Consistency

> **Reviewer 3:** *"В калибровке chr17 включена в выгрузку test, но исключена из расчёта test-статистики... Опубликованные средние дополнительно противоречат друг другу... Объяснение компенсации Null2 через threshold − 0.05 имеет неправильный знак."*

### Response:
We thank the reviewer for this rigorous mathematical verification.

### 2.1 Resolution of the Statistical Discrepancy & Chromosome 17 Inclusion
In commit `8b10936`, `test_chroms` in `calibrate_scoring.py` had accidentally omitted `"NC_060941.1"` (chr17) during the print calculation, even though chr17 was assigned `split = "test"` in the TSV export.

We corrected this in [`benchmarks/calibrate_scoring.py`](file:///Users/akomissarov/Dropbox/workspace/new/biology/humas_hmmer/benchmarks/calibrate_scoring.py) and re-ran calibration across all **81,661 exact-boundary, model-matched monomers**:
- **Training Set (chr1, chr3, chr8; $n = 23,506$):**
  - Mean ratio: **$0.8495 \pm 0.0074$** (Sample SD), SEM: $0.000048$, Median: **$0.8504$**
  - IQR: $[0.8469, 0.8534]$, 95% range: $[0.8305, 0.8599]$
- **Independent Validation Set (chr10, 11, 12, 14, 17, 22, Y; $n = 58,155$):**
  - Mean ratio: **$0.8505 \pm 0.0078$** (Sample SD), SEM: $0.000032$, Median: **$0.8522$**
  - IQR: $[0.8478, 0.8555]$, 95% range: $[0.8299, 0.8595]$
- **Score-Stratified Tiers (All Chromosomes):**
  - Canonical High Score ($\ge 150$ bits, $n = 80,341$, 98.4%): Mean **$0.8507 \pm 0.0065$** (Median $0.8517$)
  - Mid Score (100–150 bits, $n = 1,320$, 1.6%): Mean **$0.8198 \pm 0.0128$** (Median $0.8210$)

**Verification of Internal Consistency:**
$$\bar r_{\text{train/test}} = \frac{23,506 \cdot 0.849491 + 58,155 \cdot 0.850499}{81,661} = \mathbf{0.850209}$$
$$\bar r_{\text{score tiers}} = \frac{80,341 \cdot 0.850709 + 1,320 \cdot 0.819812}{81,661} = \mathbf{0.850208}$$
Both partitions yield an identical overall mean of **0.85021** (difference $< 10^{-6}$).

### 2.2 Correction of the `threshold - 0.05` Explanation
The reviewer is mathematically 100% correct:
- At an empirical ratio of $0.8198$, multiplying by $0.85$ yields an approximate score that is $\frac{0.85}{0.8198} - 1 \approx +3.68\%$ **higher** than HMMER's score.
- Checking $\text{score} \ge \text{threshold} - 0.05$ lowers the acceptance bar further in the same direction; it does not "compensate" for an over-score.
- **The True Role in Code:** The check `hmmer_sc_per_len >= score_threshold - 0.05` is an **upstream inclusion pre-filter** designed to accommodate downstream AWK formatting. In legacy HumAS-HMMER, `hmmertblout2bed.awk` computes density as `sprintf("%.1f", score / length)` and filters for $\ge 0.70$. Because standard 1-decimal rounding rounds values such as $0.65$ up to $0.70$, checking $\text{threshold} - 0.05$ ensures that candidates that would satisfy the downstream rounding threshold are not prematurely pruned during Stage 2 search.
- We have completely rewritten this explanation in `paper.md`, `EXPL.md`, and this response.

---

## 3. Active HOR Array on chr22 & Investigation of Gaps > 3 bp

> **Reviewer 3:** *"Активный массив chr22 теперь выбран правильно, но новые региональные результаты не опубликованы в проверяемом виде... В коммите нет регионального BED, журнала или скрипта, воспроизводящего 17 145 / 17 073 / 65 / 6. И даже наличие шести зазоров >3 bp само по себе не устанавливает, что это шесть реальных биологических вставок. Кроме того, 99,59% — конечная непрерывность аннотации после всего pipeline... Это число нельзя целиком приписывать эффективности самого healing."*

### Response:

### 3.1 Dedicated Reproducible Script
We created [`benchmarks/analyze_active_hor_chr22.py`](file:///Users/akomissarov/Dropbox/workspace/new/biology/humas_hmmer/benchmarks/analyze_active_hor_chr22.py), which extracts `NC_060946.1:12,788,180–15,711,065` directly from the whole-genome BEDs and computes all metrics. Its output is saved in [`benchmarks/logs/chr22_active_hor_continuity.summary.txt`](file:///Users/akomissarov/Dropbox/workspace/new/biology/humas_hmmer/benchmarks/logs/chr22_active_hor_continuity.summary.txt):

```text
==========================================================================================================
Active Centromeric HOR Array Analysis: chr22 (NC_060946.1:12,788,180–15,711,065, 2.92 Mb)
CenSat Domain: hor_22_9 (S2C14/22H1L)
==========================================================================================================
Metric                                  | FastHumAS                 | Legacy HumAS-HMMER
----------------------------------------------------------------------------------------------------------
Total Annotated Monomers                | 17145                     | 17150                    
Total Inter-Monomer Junctions           | 17144                     | 17149                    
Flush Contiguous Junctions (gap = 0 bp) | 17073 ( 99.59%)           | 0 (  0.00%)
Micro-gaps (1–3 bp)                     | 65 (  0.38%)              | 16525 ( 96.36%)
  - 1-bp gaps                           | 54                        | 0                        
  - 2-bp gaps                           | 8                         | 13773                    
  - 3-bp gaps                           | 3                         | 2752                     
Gaps > 3 bp                             | 5                         | 624                      
Overlaps (< 0 bp)                       | 1                         | 0                        
==========================================================================================================
```
Sum: $17,073 + 65 + 5 + 1 = 17,144$ junctions across 17,145 monomers.

### 3.2 Sequence Characterization of Gaps > 3 bp
We investigated the exact genomic sequences at all 5 non-flush gaps > 3 bp:
1. **Locus `[12,789,031, 12,789,250)` (219 bp):** Between `S2C14/22H1L.3` and `S2C14/22H1L.2`. This is a degenerate pericentromeric boundary region where legacy HumAS-HMMER called two low-scoring sub-monomer fragments (`S2C9H1L.7` [97.9 bits, 118 bp] and `S2C14/22H1L.3` [90.4 bits, 95 bp]). FastHumAS discards these sub-threshold fragments, leaving the 219-bp boundary unannotated.
2. **Loci `[15,478,940, 15,478,971)` (31 bp), `[15,480,341, 15,480,372)` (31 bp), and `[15,488,562, 15,488,593)` (31 bp):** These three identical 31-bp gaps occur between monomer `S2C14/22H1L.7` and monomer `S2C14/22H1L.6`. In legacy HumAS-HMMER, `nhmmer` split monomer `.6` into two fragmented sub-hits: a 58-bp fragment (score 46.4) and a 145-bp fragment (score 161.8). FastHumAS's affine Gotoh DP aligns the full canonical 177-bp monomer, cleanly preserving the authentic 31-bp structural variant insertion.
3. **Locus `[15,709,700, 15,709,706)` (6 bp):** Between `S2C14/22H1L.4` and `S2C14/22H1L.3`. This is a genuine 6-bp assembly micro-gap present in both legacy (8 bp) and FastHumAS (6 bp).
4. **Overlap `[12,789,420, 12,789,421)` (-1 bp):** A 1-bp coordinate overlap between adjacent boundary monomers.

### 3.3 Pipeline Continuity vs Standalone Healing
We agree with the reviewer that **99.59%** is the cumulative continuity of the **complete FastHumAS annotation pipeline**, reflecting the combined contributions of:
1. Coordinate rectification (eliminating the 2-bp reverse-strand shift).
2. Profile-anchored banded dynamic programming (preventing artificial monomer fragmentation into sub-hits).
3. Local NMS overlap filtering (preserving non-overlapping adjacent monomers).
4. Consensus-anchored junction healing (closing 1–3 bp micro-gaps).
We have clarified this distinction explicitly in the manuscript.

---

## 4. Dynamic 10-Chromosome Healing Analysis & Assembly-Wide Extrapolation

> **Reviewer 3:** *"Десятихромосомный healing-анализ фактически не обновлён... junction_healing_analysis.summary.txt и analyze_junction_healing.py побайтно совпадают с предыдущим коммитом... 9 хромосом, без chr17... Оценка ~190 587 остаётся экстраполяцией старого девятихромосомного результата... Сам скрипт всё ещё анализирует только геометрию BED и печатает готовое биологическое заключение."*

### Response:
We have completely rewritten [`benchmarks/analyze_junction_healing.py`](file:///Users/akomissarov/Dropbox/workspace/new/biology/humas_hmmer/benchmarks/analyze_junction_healing.py) to compute all numbers, proportions, and extrapolations dynamically from the input BED files, eliminating all hardcoded strings.

Running the updated script on the 10 benchmark chromosomes yields [`benchmarks/logs/junction_healing_analysis.summary.txt`](file:///Users/akomissarov/Dropbox/workspace/new/biology/humas_hmmer/benchmarks/logs/junction_healing_analysis.summary.txt):

```text
==========================================================================================================
Inter-Monomer Gap and Healing Analysis: Legacy HumAS-HMMER vs FastHumAS
==========================================================================================================
Chromosome   | Legacy Recs  | Legacy Gap=0 | Legacy 1-3bp | Legacy 2bp(-) | Fast Gap=0   | Fast 1-3bp   | Fast 2bp(-) 
----------------------------------------------------------------------------------------------------------
NC_060925.1  | 30571        | 10543        | 17119        | 9550          | 30039        | 232          | 38          
NC_060927.1  | 15589        | 1252         | 12805        | 11056         | 14865        | 404          | 48          
NC_060932.1  | 15804        | 12518        | 2182         | 223           | 14713        | 506          | 11          
NC_060934.1  | 16915        | 10753        | 3925         | 260           | 14378        | 1153         | 149         
NC_060935.1  | 26385        | 18912        | 4168         | 609           | 22507        | 2975         | 617         
NC_060936.1  | 20107        | 599          | 16887        | 14080         | 18156        | 1282         | 286         
NC_060938.1  | 22434        | 1857         | 18476        | 14009         | 21205        | 691          | 57          
NC_060941.1  | 25412        | 23735        | 1339         | 81            | 24994        | 225          | 6           
NC_060946.1  | 24531        | 694          | 21952        | 17501         | 23717        | 487          | 53          
NC_060948.1  | 2589         | 92           | 2047         | 923           | 2196         | 197          | 70          
==========================================================================================================
TOTAL (10 chr) | 200337       | 80955        | 100900       | 68292         | 186770       | 8152         | 1335        
```

**Key Findings (Dynamically Computed across 10 Chromosomes):**
1. **Legacy Micro-Gaps:** Legacy HumAS-HMMER contains **100,900 micro-gaps (1–3 bp)** (50.37% of all monomers).
2. **Minus-Strand Coordinate Bug:** Of these micro-gaps, **68,292 (67.68%)** are exact 2-bp gaps on the minus strand caused by `aliFrom = $7 - 1; aliTo = $8`.
3. **Assembly-Wide Extrapolation:** Scaled to the 488,904 monomers across whole-genome CHM13:
   $$68,292 \times \frac{488,904}{200,337} = \mathbf{166,660}\text{ artificial gaps genome-wide}.$$
4. **FastHumAS Reduction:**
   - Flush contiguous junctions increase from **40.41%** (80,955) in legacy to **93.72%** (186,770) in FastHumAS.
   - Micro-gaps (1–3 bp) decrease by **91.92%** (from 100,900 down to 8,152).
   - Artificial 2-bp minus-strand gaps decrease by **98.05%** (from 68,292 down to 1,335).

---

## 5. Concordance Summary Synchronization & CLI Ratios

> **Reviewer 3:** *"Таблица в сообщении и опубликованная concordance-сводка содержат разные результаты... Заявление об исправленном выводе точного strand agreement тоже пока не соответствует скрипту: он считает total_strand, но печатает только округлённый процент."*

### Response:
We updated [`benchmarks/compare_fasthumas_vs_nhmmer.py`](file:///Users/akomissarov/Dropbox/workspace/new/biology/humas_hmmer/benchmarks/compare_fasthumas_vs_nhmmer.py) to explicitly print the raw counts and ratios in the standard CLI output.

Running the updated comparator on the fixed whole-genome output against all 10 finished legacy chromosomes yields [`benchmarks/logs/nhmmer_concordance.summary.txt`](file:///Users/akomissarov/Dropbox/workspace/new/biology/humas_hmmer/benchmarks/logs/nhmmer_concordance.summary.txt):

```text
================================================================================
Chrom          | FastHumAS  | Legacy   | Overlap  | Sens %  | Prec %  | <=3bp Bdry  | Label Agr  | Strand Agr
--------------------------------------------------------------------------------
NC_060925.1    | 30514      | 30571    | 30514    |  99.81% | 100.00% |  98.95% (30193) |  99.78% | 100.00%
NC_060927.1    | 15493      | 15589    | 15493    |  99.38% | 100.00% |  98.24% (15221) |  98.83% | 100.00%
NC_060932.1    | 15710      | 15804    | 15673    |  99.17% |  99.76% |  97.26% (15244) |  99.55% | 100.00%
NC_060934.1    | 16838      | 16915    | 16824    |  99.46% |  99.92% |  98.08% (16501) |  99.54% | 100.00%
NC_060935.1    | 26173      | 26385    | 26164    |  99.16% |  99.97% |  96.86% (25343) |  99.15% | 100.00%
NC_060936.1    | 19951      | 20107    | 19945    |  99.19% |  99.97% |  97.27% (19400) |  99.20% | 100.00%
NC_060938.1    | 22275      | 22434    | 22271    |  99.27% |  99.98% |  94.94% (21144) |  98.89% | 100.00%
NC_060941.1    | 25360      | 25412    | 25356    |  99.78% |  99.98% |  99.34% (25189) |  99.61% | 100.00%
NC_060946.1    | 24434      | 24531    | 24431    |  99.59% |  99.99% |  97.10% (23722) |  99.14% | 100.00%
NC_060948.1    | 2554       | 2589     | 2549     |  98.46% |  99.80% |  96.59% (2462) |  99.65% | 100.00%
================================================================================
TOTAL / MEAN   | 199302     | 200337   | 199220   |  99.44% |  99.96% |  97.59% (194419) |  99.32% | 100.00%

Exact Boundary (=0 bp): 81909 / 199220 (41.11%)
Boundary within <=3 bp: 194419 / 199220 (97.59%)
Boundary within <=10 bp: 198711 / 199220 (99.74%)
Strand Orientation Concordance: 199220 / 199220 (100.00%)
Subfamily / Label Concordance: 197874 / 199220 (99.32%)
FastHumAS unique (not in legacy): 82 / 199302 (0.04%)
Legacy unique (not in FastHumAS): 1117 / 200337 (0.56%)
```

All figures in `paper.md`, `EXPL.md`, and this response are strictly synchronized with this verified output.

---

## 6. Inclusion of Manuscript, Explanations, and Review Materials in the Repository

> **Reviewer 3:** *"В дереве коммита отсутствуют paper.md, paper.html, EXPL.md и response_to_reviewers.md. Поэтому заявленные исправления новых редакций этих документов по опубликованному репозиторию проверить невозможно."*

### Response:
We have moved all documentation directly into the repository under the [`paper/`](file:///Users/akomissarov/Dropbox/workspace/new/biology/humas_hmmer/paper) directory:
- [`paper/paper.md`](file:///Users/akomissarov/Dropbox/workspace/new/biology/humas_hmmer/paper/paper.md): Full manuscript source in GitHub Flavored Markdown.
- [`paper/paper.html`](file:///Users/akomissarov/Dropbox/workspace/new/biology/humas_hmmer/paper/paper.html): Standalone HTML version with MathJax equations.
- [`paper/EXPL.md`](file:///Users/akomissarov/Dropbox/workspace/new/biology/humas_hmmer/paper/EXPL.md): Complete algorithmic derivations, formulas, and benchmark documentation.
- [`paper/response_to_reviewers.md`](file:///Users/akomissarov/Dropbox/workspace/new/biology/humas_hmmer/paper/response_to_reviewers.md): This complete point-by-point response.

A prominent navigation section has been added to [`README.md`](file:///Users/akomissarov/Dropbox/workspace/new/biology/humas_hmmer/README.md) linking directly to these documents.
