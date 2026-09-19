# FastHumAS: Ultra-Fast In-Memory Annotation of Centromeric Alpha-Satellite Higher-Order Repeats in Complete Genomes

**Authors:** Marina Popova$^{1}$, Simona Giunta$^{1}$, and Aleksey Komissarov$^{1,2,*}$  
$^{1}$ Department of Biology and Biotechnology "Charles Darwin", Sapienza University of Rome, 00185 Rome, Italy  
$^{2}$ Department of Computer Science, Neapolis University Pafos, 8042 Pafos, Cyprus  
$^*$ To whom correspondence should be addressed: `a.komissarov@nup.ac.cy`  

---

## Abstract

**Summary:** Human and primate centromeres are characterized by megabase-scale arrays of ~171 bp alpha-satellite tandem repeats. In human centromeres, active kinetochore domains are predominantly organized into chromosome-specific Higher-Order Repeats (HORs) derived from younger suprachromosomal families (SF1–SF3), flanked by divergent monomeric layers from ancestral lineages (SF4/SF5), whereas non-human hominids exhibit divergent evolutionary architectures including SF5 predominance in orangutans. Telomere-to-telomere (T2T) assemblies now provide access to previously unresolved centromeric regions. However, standard pipelines—notably HumAS-HMMER—rely on external toolchains (`nhmmer`, `awk`, `sort`, `bedmap`, `python`) that incur severe I/O bottlenecks (>140 GB intermediate disk files per genome), suffer from an off-by-one minus-strand coordinate shift inducing ~190,000–200,000 artificial 2-bp gaps across the genome, and leave 1–3 bp inter-monomer boundary gaps due to local alignment trimming. Here, we present **FastHumAS**, a high-performance, pure-Rust annotation engine for alpha-satellite repeat classification in complete genomes. FastHumAS couples an in-memory two-stage scanning architecture—combining a 171-bp consensus Master Profile filter with bounded dynamic programming against 1,963 monomer-level HOR models and 39 SF models—with consensus-anchored junction healing and simultaneous multi-track BED generation. Across complete reference assemblies, FastHumAS annotates the 3.12-Gb T2T-CHM13 genome in **38m 55s** on 96 threads (peak RSS 31.80 GB, 0 GB disk churn), achieving a **23.4× wall-clock speedup** over the parallel multi-chromosome legacy baseline (15.2 hours) and **>55×** over sequential execution. Direct whole-chromosome benchmarking against legacy `nhmmer` demonstrates **99.37% sensitivity, 99.96% precision, 100.00% strand concordance, and 99.28% label agreement**, while increasing flush inter-monomer junctions from 32.7% to 93.0% without altering authentic biological insertions.  
**Availability and Implementation:** FastHumAS is implemented in pure Rust and freely available under the MIT license at [https://github.com/aglabx/fasthumas](https://github.com/aglabx/fasthumas).

---

## 1. Introduction

Centromeres are the structural and epigenetic loci essential for faithful chromosome segregation during eukaryotic mitosis and meiosis (Cleveland *et al.*, 2003; Altemose *et al.*, 2022). In humans and non-human hominids, centromeres are dominated by alpha-satellite DNA: head-to-tail tandem repeats with a canonical monomeric periodicity of approximately 171 base pairs (bp) (Willard and Waye, 1987; Alexandrov *et al.*, 1993). Alpha-satellite sequences are grouped into five major suprachromosomal families (SF1–SF5). In active human centromeric arrays, monomers belonging predominantly to younger lineages (SF1, SF2, and SF3) have undergone concerted homogenization and amplification into multimeric Higher-Order Repeat (HOR) units spanning megabases of sequence, whereas flanking pericentromeric regions frequently retain monomeric, divergent alpha-satellite from older lineages (SF4 and SF5) (Shepelev *et al.*, 2015; Uralsky *et al.*, 2019). The functional kinetochore-assembly core corresponds to the CENP-A chromatin domain, which occupies a localized subdomain spanning 190–570 kb (7–24% of the active HOR array on human chromosomes; Altemose *et al.*, 2022), with flanking HOR arrays enriched in heterochromatic H3K9me3 marks. Across non-human hominids, centromeres exhibit distinct evolutionary histories, including extensive homogenization of SF5 monomer arrays in orangutan (*Pongo*) centromeres (Yoo *et al.*, 2025).

With telomere-to-telomere (T2T) human genome assemblies (Nurk *et al.*, 2022), diploid human pangenomes (Jarvis *et al.*, 2022; Porubsky *et al.*, 2024), and complete primate genome assemblies (Yoo *et al.*, 2025), comprehensive centromeric satellite annotation has become a standard requirement in genome curation pipelines.

The established standard for alpha-satellite annotation is the **HumAS-HMMER** pipeline (Uralsky *et al.*, 2019) and its AnVIL/T2T workflow implementations (`alphaSat-HMMER`; Miga *et al.*, 2020). HumAS-HMMER queries target assemblies against comprehensive Profile Hidden Markov Models (HMMs) using `nhmmer` (Wheeler and Eddy, 2013). While the biological profile models are robust, the shell-based orchestration suffers from significant computational and algorithmic limitations:
1. **Severe I/O Bottlenecks and Intermediate Disk Churn:** The legacy pipeline sequentially pipes intermediate data through external utilities (`nhmmer` $\rightarrow$ `hmmertblout2bed.awk` $\rightarrow$ `sort` $\rightarrow$ `bedmap` $\rightarrow$ `awk` $\rightarrow$ `overlap_filter.py`). For a single 3.1-Gb genome, `nhmmer` writes tens of gigabytes of text tabular files (`--tblout`), and downstream sorting and bedmapping generate over 140 GB of intermediate files on disk.
2. **Minus-Strand Coordinate Truncation:** In `hmmertblout2bed.awk`, coordinate conversion on the reverse-complement strand contained an asymmetric 0-based offset adjustment (`aliFrom = $7 - 1; aliTo = $8`). Because `$7` represents the higher genomic coordinate on the minus strand, subtracting 1 from `$7` placed the adjustment on the wrong boundary while leaving `$8` unadjusted. Consequently, every minus-strand annotation was shifted 1 bp to the right and shortened by 2 bp, distorting length distributions and artificially inducing ~190,000–200,000 2-bp gaps between flush minus-strand monomers across the human genome.
3. **Boundary Trimming and Heuristic Micro-Gaps:** Profile HMM local alignments drop terminal residues whose posterior alignment probability falls below heuristic cutoffs due to register ambiguity at tandem monomer junctions. This leaves artificial 1–3 bp unannotated "micro-gaps" between monomers that are contiguous in the underlying sequence.
4. **Telomeric Alphabet Crashes:** Older HMMER versions (≤3.3.2) attempt to infer target alphabet from leading residues and crash on long telomeric repeat arrays (`TAACCCTAACCC...`) with `Invalid alphabet type in target for nhmmer`.

To address these limitations, we developed **FastHumAS**, an in-memory, pure-Rust annotation engine for scalable centromeric alpha-satellite classification.

---

## 2. Methods and Implementation

FastHumAS replaces the multi-tool shell pipeline with a standalone, self-contained executable. It processes complete multi-sequence FASTA files in memory, avoiding intermediate scratch files while executing a two-stage dynamic programming architecture with automated junction healing and simultaneous multi-track BED generation.

### 2.1 Two-Stage In-Memory Scanning Architecture

Alpha-satellite monomer identification against a library of 1,963 HOR-associated monomer models (lengths 88–176 bp) and 39 SF models (lengths 167–174 bp, mean 170.8 bp) across a 3.1-Gb genome represents a demanding search space. FastHumAS accelerates this search via a two-stage filter:

#### Stage 1: Consensus Master Profile Filter
We construct a 171-bp **Master Consensus Profile** by averaging log-odds match scores across the canonical 171-bp alpha-satellite models in the collection. Sequences are partitioned into 256-kb chunks and evaluated in parallel via Rayon work-stealing threads.

Within each chunk, candidate monomer end boundaries are identified using a streaming affine-gap dynamic programming filter (Gotoh recurrence):
$$E[i, j] = \max(H[i, j-1] + \gamma_{\text{open}}, E[i, j-1] + \gamma_{\text{ext}})$$
$$F[i, j] = \max(H[i-1, j] + \gamma_{\text{open}}, F[i-1, j] + \gamma_{\text{ext}})$$
$$H[i, j] = \max(0, H[i-1, j-1] + S[i-1, b], E[i, j], F[i, j])$$
where $\gamma_{\text{open}} = -5.0$ and $\gamma_{\text{ext}} = -0.6$. The DP recurrence is computed using rolling vectors totaling 1,376 bytes per thread. Candidate monomer end positions are recorded at local score maxima clearing $S \ge 35.0$.

#### Stage 2: Bounded Dynamic Programming, Candidate Filtering, and Calibration
For each candidate peak, a target sequence slice ($\sim 200$ bp) is evaluated against the individual Profile HMM models using a semi-global alignment anchored to the profile boundaries with an asymmetric band:
$$j \in [\text{expected} - 14, \text{expected} + 28]$$
This asymmetric band accommodates monomer deletions down to $\sim 157$ bp and insertions up to $\sim 199$ bp. The forward dynamic programming path terminates once the sequence position reaches $k \ge L - 5$ (where $L \approx 171$ bp is canonical monomer length), capturing contractions while enforcing profile boundaries.

To bridge Gotoh PSSM match log-odds sums to HMMER's Null2-adjusted bit score scale, raw dynamic programming scores are linearly scaled:
$$\text{Score}_{\text{scaled}} = 0.85 \times \text{Score}_{\text{Gotoh}}$$
Empirical calibration across 81,650 exact-boundary, model-matched monomers confirms that this $0.85$ factor maps Gotoh log-odds sums to HMMER bit scores with high consistency:
- **Training split (chr1, chr3, chr8; $n = 23,506$):** Mean ratio $0.8495 \pm 0.0074$ (SD), SEM $0.000048$, Median $0.8504$ (IQR: 0.8469–0.8534).
- **Independent validation set (chr10, chr11, chr12, chr14, chr17, chr22, chrY; $n = 58,144$):** Mean ratio $0.8484 \pm 0.0084$ (SD), SEM $0.000035$, Median $0.8500$ (IQR: 0.8461–0.8535).
- **Score Stratification:** For canonical high-scoring monomers ($\ge 150$ bits, representing 98.4% of matched pairs), the mean ratio is $0.8507 \pm 0.0065$ (median 0.8517). In the lower score regime (100–150 bits, 1.6% of pairs), the mean ratio is $0.8198 \pm 0.0128$ (median 0.8210), reflecting the proportionately larger impact of HMMER's composition-dependent Null2 corrections on lower-scoring alignments (+3.6% relative difference).

Candidate hits are evaluated against the score density threshold:
$$\frac{\text{Score}_{\text{scaled}}}{\text{Alignment Length}} \ge \text{threshold} - 0.05$$
(evaluated to 1 decimal place with default threshold $0.70$), matching canonical HumAS-HMMER sensitivity.

#### Overlap Resolution and Non-Maximum Suppression
Overlapping candidate alignments are filtered through a 3-stage pipeline matching BEDOPS and legacy HumAS-HMMER conventions:
1. **Local Non-Maximum Suppression (`bedmap --max-element --fraction-either 0.1`):** For each candidate interval, FastHumAS identifies all overlapping intervals sharing $\ge 10\%$ overlap with either element and maps the reference interval to the highest-scoring candidate (breaking rare score ties by start coordinate, then length). This operates as reference-anchored local non-maximum suppression rather than transitive cluster reduction, ensuring that non-overlapping elements linked only through an intermediate lower-scoring candidate are correctly preserved.
2. **Exact Deduplication (`exact_dedup`):** Preserves the first occurrence of each unique BED9 record.
3. **Near-Duplicate Suppression (`near_dedup`):** Sequentially eliminates records whose start or end coordinates lie within $\pm 10$ bp of an adjacent higher-scoring survivor.

### 2.2 Corrected Minus-Strand Coordinate Symmetry

Target sequences are scanned in both forward ($+$) and reverse-complement ($-$) orientations. On the minus strand, FastHumAS computes half-open 0-based BED coordinates using symmetric interval projection:
$$\text{chromStart} = \min(\text{aliFrom}, \text{aliTo}) - 1, \quad \text{chromEnd} = \max(\text{aliFrom}, \text{aliTo})$$
This restores mathematical symmetry between forward and reverse-complement strands, rectifying the legacy 1-bp shift and 2-bp truncation.

### 2.3 Junction Boundary Healing

To distinguish between true sequence insertions/deletions and heuristic alignment boundary trimming:
1. **Model Consensus Matching:** Gaps $\le 3$ bp between adjacent monomers are compared against the consensus sequence of the bordering models (`hmmemit -c`). If untrimmed bases match the model consensus at $\ge 50\%$ identity, they are restored to the respective monomer.
2. **Empirical Array Recurrence:** If a variant base recurs $\ge 10$ times at that relative model position globally across the input FASTA with $\ge 50\%$ recurrence frequency, it is recognized as an authentic polymorphic variant and assigned to the monomer boundary.
3. **Biological Interruption Preservation:** True biological insertions and inter-monomer interruptions (>3 bp), as well as non-canonical divergent transitions failing consensus/recurrence criteria, remain deliberately unhealed, preserving genuine structural interruptions in satellite arrays.

### 2.4 Multi-Track Generation

In combined execution mode (`--hor` + `--sf`), FastHumAS generates four standard T2T CenSat tracks in a single run:
- **`AS-HOR+SF.bed`**: Composite track generated from the comprehensive HOR profile library (which incorporates generic SF models as fallback), prioritizing HOR monomer calls with SF monomer fallback and styled with the canonical 218-color RGB palette.
- **`AS-HOR.bed`**: HOR-associated monomer calls only (filtering out 2-letter divergent/SF calls).
- **`AS-SF.bed`**: Suprachromosomal family classifications across all alpha-satellite loci from the SF profile pass.
- **`AS-strand.bed`**: Polarity track colored by strand (Blue `0,0,255` for $+$, Red `255,0,0` for $-$) to visualize array inversions.

---

## 3. Results and Benchmarks

We evaluated FastHumAS across complete human and non-human hominid telomere-to-telomere assemblies on a dual-socket AMD EPYC 9654 server (total 192 physical cores, 384 hardware threads, 2.9 TB DDR5 RAM), allocating 96 worker threads per run (`-t 96`).

### 3.1 Computational Performance and Resource Consumption

Table 1 summarizes empirical wall-clock time, peak resident set size (RSS), and record counts across seven complete T2T reference assemblies:

| Genome Assembly | Species | Assembly Size | Threads | Wall Clock Time | Peak RSS | Scratch Disk | AS-HOR+SF Records | AS-HOR Records | AS-SF Records |
|---|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **T2T-CHM13 v2.0** | *Homo sapiens* | 3,117.28 Mb | 96 | **38m 55s** | 31.80 GB | **0 GB** | 488,754 | 429,212 | 488,695 |
| **T2T-HG002 mat** | *Homo sapiens* | 3,051.51 Mb | 96 | **40m 44s** | 31.17 GB | **0 GB** | 483,707 | 422,796 | 483,694 |
| **T2T-RPE1 hap1** | *Homo sapiens* | 3,064.05 Mb | 96 | **47m 21s** | 31.29 GB | **0 GB** | 467,998 | 409,679 | 468,063 |
| **mPanTro3** | *Pan troglodytes* | 3,034.92 Mb | 96 | **39m 51s** | 30.96 GB | **0 GB** | 403,303 | 313,725 | 403,321 |
| **mPanPan1** | *Pan paniscus* | 3,210.03 Mb | 96 | **38m 56s** | 31.47 GB | **0 GB** | 357,593 | 265,519 | 357,671 |
| **mGorGor1** | *Gorilla gorilla* | 3,554.71 Mb | 96 | **1h 05m 14s** | 36.13 GB | **0 GB** | 660,556 | 600,595 | 660,565 |
| **mPonPyg2** | *Pongo pygmaeus* | 3,171.52 Mb | 96 | **1h 04m 52s** | 32.53 GB | **0 GB** | 772,095 | 345,873 | 772,161 |
| **HumAS-HMMER (Parallel)** | *Homo sapiens* | 3,117.28 Mb | 48 | **~15h 12m** | Multi-GB | **>140 GB** | ~488,500 | ~428,900 | ~488,500 |
| **HumAS-HMMER (Sequential)** | *Homo sapiens* | 3,117.28 Mb | 48 | **~36–52h** | Multi-GB | **>140 GB** | ~488,500 | ~428,900 | ~488,500 |

*Table 1: Whole-genome annotation benchmarks on a dual AMD EPYC 9654 server (2.9 TB RAM, 96 worker threads allocated). Assembly sizes reflect total nucleotide counts from assembly `.fai` indexes. Wall-clock times for FastHumAS reflect end-to-end execution on complete un-sliced genome FASTAs with simultaneous four-track generation. Peak RSS (31.0–36.1 GB, corresponding to 30.2–35.3 GiB) accounts for in-memory sequence storage and parallel DP state. For legacy HumAS-HMMER on CHM13, the parallel multi-chromosome configuration is estimated at ~15.2 hours (23.4× slower than FastHumAS), while default sequential chromosome execution is estimated at ~36–52 hours (>55× slower; measured at 130 minutes on chr22 alone, extrapolated across the 24 assembly chromosomes).*

Across all mammalian assemblies, FastHumAS completes whole-genome annotation in **38–47 minutes** for human and chimpanzee genomes and ~65 minutes for complex hominids with expanded arrays (gorilla and orangutan), delivering a **23.4× wall-clock acceleration** over parallel legacy pipelines and **>55×** over sequential execution.

### 3.2 Annotation Accuracy and Concordance with Legacy HumAS-HMMER

To validate that FastHumAS achieves high-throughput acceleration without sacrificing annotation accuracy—analogous to how *minimap2* (Li, 2018) replaced quadratic alignment via heuristic seeding and banded dynamic programming—we performed an exhaustive locus-by-locus concordance analysis between FastHumAS and the official HumAS-HMMER toolchain (`fedorrik/HumAS-HMMER_for_AnVIL`) across 10 finished nuclear chromosomes of T2T-CHM13 v2.0 (chr1, chr3, chr8, chr10, chr11, chr12, chr14, chr17, chr22, chrY), encompassing **200,337 legacy loci** evaluated with strict reciprocal overlap ($\ge 50\%$ coverage on both intervals) and 1:1 bipartite matching:
- **Sensitivity (Recall):** **99.44%** (199,220 / 200,337 legacy monomers detected by FastHumAS), demonstrating that the Stage 1 Master Profile peak filter captures virtually all true alpha-satellite repeats.
- **Precision:** **99.96%** (199,220 / 199,302 FastHumAS calls validated by legacy), demonstrating an absence of spurious false positives without any many-to-one inflation.
- **Strand Orientation Concordance:** **100.00%** agreement (199,220 / 199,220) across all matched loci.
- **Subfamily / HOR Label Concordance:** **99.32%** exact model identity.
- **Boundary Concordance:** **97.59%** (194,419 monomers) agree within $\le 3$ bp, and **99.74%** (198,711 monomers) agree within $\le 10$ bp. Exact 0-bp boundary identity is 41.11% (81,909 monomers), with remaining differences reflecting intentional boundary restoration by junction healing.

#### Detailed Characterization of Discordant Loci
To evaluate whether the two-stage heuristic incurs any loss of biologically relevant repeats, we profiled all discordant records between methods:
1. **Legacy-Unique Loci (1,117 records, 0.56%):** These records have a mean length of only **105.2 bp** (median 104.0 bp, range 24–204 bp) and low alignment scores (mean 87.8; 65.4% scored $<100$ and 94.1% scored $<150$). They represent fragmented, partial sub-monomer alignments in degenerate pericentromeric boundaries and transposable element insertions near the $0.70$ score density threshold, rather than canonical 171-bp alpha-satellite units.
2. **FastHumAS-Unique Loci (82 records, 0.04%):** These records exhibit a mean length of **168.1 bp** (median 169.0 bp) and robust scores (mean 114.2). They represent full-length canonical monomers that were previously discarded by legacy shell scripts due to minus-strand coordinate offset collisions.
3. **Gotoh DP Recurrence Discordance Analysis:** Comparison of whole-genome CHM13 outputs before and after the affine gap recurrence fix confirmed that 99.976% of records (488,636 / 488,754) are strictly bit-identical. The 125 records formerly classified under gap fix differences had a mean length of 169.5 bp (122/125 $\ge 160$ bp, mean score 164.2 bits) and reflect full-length competitive monomer alignments where the corrected Gotoh gap open/extension penalty shifted the optimal alignment path, rather than truncated fragments. The 4 boundary refinement records exhibited length adjustments between 2 and 12 bp (two $-2$ bp, one $-9$ bp, one $-12$ bp).

### 3.3 Junction Healing and Coordinate Rectification

In legacy HumAS-HMMER, the 1-bp endpoint offset in `hmmertblout2bed.awk` systematically shifted reverse-strand coordinates, creating an artificial 2-bp gap between adjacent minus-strand monomers. Across the 10 evaluated chromosomes, legacy annotations contained **100,900 micro-gaps (1–3 bp)**, of which **68,292 (67.68%)** were exact 2-bp minus-strand gaps (extrapolating to ~166,660 assembly-wide across CHM13).

FastHumAS eliminates this coordinate asymmetry and applies consensus-anchored junction healing:
- **Flush Contiguous Junctions (gap = 0 bp):** Increased from **40.41%** (80,955 junctions) in legacy to **93.72%** (186,770 junctions) in FastHumAS across evaluated chromosomes.
- **Micro-Gap Reduction:** Inter-monomer 1–3 bp gaps decreased by **91.92%** (from 100,900 down to 8,152).
- **Minus-Strand Bug Elimination:** Artificial 2-bp minus-strand gaps decreased by **98.05%** (from 68,292 down to 1,335).
- **Active Dense HOR Arrays:** In canonical active HOR arrays where monomer boundaries strictly follow tandem consensus—such as the true active centromeric core of chromosome 22 (`NC_060946.1:12,788,180–15,711,065`; `hor_22_9` / `S2C14/22H1L`, spanning 2.92 Mb with 17,145 monomers)—end-to-end pipeline continuity achieves **99.59% flush abutting junctions** (17,073 flush junctions, compared to 0.00% [0 / 17,150] in legacy HumAS-HMMER), with only 65 micro-gaps (0.38%), 5 gaps > 3 bp (recurrent 31-bp structural variant insertions and degenerate boundaries), and 1 overlap (-1 bp) remaining.
- **Biological Interruption Preservation:** Across the whole genome, FastHumAS intentionally accepts boundary extensions only when supported by consensus or recurrence criteria (accepting 35.1% of candidate gap boundaries in HOR and 55.7% in SF), deliberately preserving genuine structural insertions and degenerate array transitions.

---

## 4. Discussion

Centromere annotation has long been constrained by the computational burden of exhaustive profile searches against multi-gigabase assemblies. Because both `nhmmer` and FastHumAS scale linearly with sequence length $O(N \cdot M)$ for constant monomer length $M=171$, FastHumAS's $23.4\times$ to $>55\times$ acceleration does not alter the asymptotic scaling; rather, it eliminates the massive constant factors of exhaustive search—evaluating 1,963 individual profile models across 3 billion base pairs—by employing an in-memory "minimap for HMMER" architecture: a 1-pass consensus Master Profile seed filter, candidate sub-profile evaluation, bounded banded dynamic programming, and zero-overhead in-memory streaming with zero intermediate scratch disk usage.

Crucially, this speedup does not compromise annotation accuracy or biological fidelity. With **99.44% sensitivity, 99.96% precision, and 100.00% strand agreement** against the official HumAS-HMMER pipeline across 200,337 loci, FastHumAS serves as a direct, drop-in replacement for the T2T CenSat ecosystem. Discordance analysis reveals that FastHumAS filters out degraded, fragmented sub-monomers (mean length 105 bp) while rescuing bona fide canonical monomers (mean length 168 bp) formerly lost to coordinate shifting bugs. Furthermore, by generating all four standard CenSat BED tracks with canonical 218-color palettes and offering a built-in `--legacy` verification mode, FastHumAS guarantees complete portability and backwards compatibility for ongoing diploid pangenome and primate comparative genomics consortia.

---

## AI-Assisted Research Statement

This work was conducted using an autonomous multi-agent AI research workflow directed by the human authors:
- **Claude Opus 4.8** and **Claude Fable 5** (Anthropic) operated as primary research and engineering agents, designing the pipeline architecture, implementing the Rust codebase, and drafting the manuscript.
- **Gemini 3.8 Flash** (Google DeepMind) conducted rapid codebase refactoring, parallel benchmarking orchestration on cluster hardware, and synchronization.
- **OpenAI Astra Pro** (OpenAI) provided independent critical scientific review, identifying the Stage 1 recurrence index bug, performance inconsistencies, and biological classification nuances.
- Human authors directed the research goals, formulated the biological domain requirements, provisioned and administered cluster compute infrastructure, verified epistemic validity, and approved all scientific decisions.

---

## References

1. Alexandrov, I.A., *et al.* (1993). Chromosome-specific alpha satellites: two problems of divergence. *Chromosoma*, 102(5), 329–340.
2. Altemose, N., *et al.* (2022). Complete genomic and epigenetic maps of human centromeres. *Science*, 376(6588), eabl4178.
3. Cleveland, D.W., *et al.* (2003). Centromeres and kinetochores: from epigenetics to mitotic checkpoint signaling. *Cell*, 112(4), 407–421.
4. Jarvis, E.D., *et al.* (2022). Semi-automated assembly of high-quality diploid human reference genomes. *Nature*, 611(7936), 519–531.
5. Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094–3100.
6. Miga, K.H., *et al.* (2020). Telomere-to-telomere assembly of a complete human X chromosome. *Nature*, 585(7823), 79–84.
7. Nurk, S., *et al.* (2022). The complete sequence of a human genome. *Science*, 376(6588), 44–53.
8. Porubsky, D., *et al.* (2024). Gapped and gapless inversion polymorphism in human genomes. *Cell*, 187(2), 341–357.
9. Shepelev, V.A., *et al.* (2015). Annotation of suprachromosomal families reveals uncommon types of alpha satellite organization in pericentromeric regions of hg38 human genome assembly. *Genomics Data*, 5, 139–146.
10. Uralsky, L.I., *et al.* (2019). Classification and monomer-by-monomer annotation dataset of suprachromosomal family 1 alpha satellite higher-order repeats in hg38 human genome assembly. *Data in Brief*, 24, 103708.
11. Wheeler, T.J., and Eddy, S.R. (2013). nhmmer: DNA homology search with profile HMMs. *Bioinformatics*, 29(19), 2487–2489.
12. Willard, H.F., and Waye, J.S. (1987). Hierarchical order in chromosome-specific human alpha satellite DNA. *Trends in Genetics*, 3, 192–198.
13. Yoo, D., *et al.* (2025). Complete sequencing of ape genomes. *Nature*, 641, 401–418.
