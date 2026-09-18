use rayon::prelude::*;

use crate::bed::BedRecord;
use crate::color::get_color;
use crate::hmm::{HmmCollection, HmmModel};
use crate::overlap::filter_overlaps;
use crate::tblout::HmmHit;

/// Maps ASCII nucleotide to index 0..3 (A, C, G, T), or 0 for unknown.
#[inline(always)]
pub fn base_to_idx(b: u8) -> usize {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 0,
    }
}

/// Reverse complement of a DNA sequence in uppercase.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let mut rev = Vec::with_capacity(seq.len());
    for &b in seq.iter().rev() {
        let comp = match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        };
        rev.push(comp);
    }
    rev
}

/// Master consensus profile built from averaging all 171-bp models.
pub struct MasterProfile {
    pub length: usize,
    pub pssm: Vec<[f32; 4]>,
}

impl MasterProfile {
    pub fn from_collection(coll: &HmmCollection) -> Self {
        let length = 171;
        let mut count = 0usize;
        let mut avg_scores = vec![[0.0f32; 4]; length];

        for m in &coll.models {
            if m.length == length {
                count += 1;
                for (avg_row, m_row) in avg_scores.iter_mut().take(length).zip(&m.match_scores) {
                    for (avg_s, &m_s) in avg_row.iter_mut().zip(m_row) {
                        *avg_s += m_s;
                    }
                }
            }
        }

        if count > 0 {
            let inv = 1.0 / count as f32;
            for row in &mut avg_scores {
                for s in row {
                    *s *= inv;
                }
            }
        } else if !coll.models.is_empty() {
            let m0 = &coll.models[0];
            let l = m0.length.min(length);
            avg_scores[..l].copy_from_slice(&m0.match_scores[..l]);
        }

        MasterProfile {
            length,
            pssm: avg_scores,
        }
    }

    /// Fast parallel scan of the sequence using the master profile to find candidate monomer end positions (1-based).
    /// Advances a single column of the Gotoh affine Smith-Waterman recurrence.
    ///
    /// Rolling buffers `h` and `e` are updated in place for target nucleotide index `b`.
    /// Returns the updated score at the end of the profile `h[l]`.
    #[inline(always)]
    pub fn step_dp(
        &self,
        b: usize,
        h: &mut [f32],
        e: &mut [f32],
        gap_open: f32,
        gap_ext: f32,
    ) -> f32 {
        let l = self.length;
        let mut prev_h = 0.0f32;
        let mut f = -99999.0f32;

        for i in 1..=l {
            let match_sc = self.pssm[i - 1][b];
            let curr_e = (h[i] + gap_open).max(e[i] + gap_ext);
            e[i] = curr_e;

            let from_diag = prev_h + match_sc;
            let new_h = (0.0f32).max(from_diag).max(curr_e).max(f);
            prev_h = h[i];
            h[i] = new_h;

            f = (new_h + gap_open).max(f + gap_ext);
        }

        h[l]
    }

    /// Aligns a sequence against the profile using the exact production DP recurrence,
    /// returning the peak score attained at the profile end.
    #[allow(dead_code)]
    pub fn align_sequence(&self, seq_indices: &[usize], gap_open: f32, gap_ext: f32) -> f32 {
        let l = self.length;
        let mut h = vec![0.0f32; l + 1];
        let mut e = vec![-99999.0f32; l + 1];
        let mut max_sc = 0.0f32;

        for &b in seq_indices {
            let sc = self.step_dp(b, &mut h, &mut e, gap_open, gap_ext);
            if sc > max_sc {
                max_sc = sc;
            }
        }

        max_sc
    }

    #[allow(clippy::needless_range_loop)]
    pub fn find_monomer_ends(&self, seq_indices: &[usize]) -> Vec<usize> {
        let l = self.length;
        let n = seq_indices.len();
        if n < l {
            return Vec::new();
        }

        const CHUNK_SIZE: usize = 256 * 1024;
        let num_chunks = n.div_ceil(CHUNK_SIZE);

        let raw_peaks: Vec<(usize, f32)> = (0..num_chunks)
            .into_par_iter()
            .flat_map_iter(|chunk_idx| {
                let start_idx = chunk_idx * CHUNK_SIZE;
                let end_idx = ((chunk_idx + 1) * CHUNK_SIZE).min(n);

                let scan_start = start_idx.saturating_sub(200);
                let scan_end = (end_idx + 2).min(n);

                let mut h = vec![0.0f32; l + 1];
                let mut e = vec![-99999.0f32; l + 1];
                let gap_open = -5.0f32;
                let gap_ext = -0.6f32;

                let mut chunk_peaks = Vec::new();

                let mut p2_sc = 0.0f32;
                let mut p1_pos = 0usize;
                let mut p1_sc = 0.0f32;

                for pos_0 in scan_start..scan_end {
                    let b = seq_indices[pos_0];
                    let curr_sc = self.step_dp(b, &mut h, &mut e, gap_open, gap_ext);
                    let curr_pos = pos_0 + 1; // 1-based coordinate

                    if p1_sc >= 35.0
                        && p1_sc >= p2_sc
                        && p1_sc >= curr_sc
                        && p1_pos > start_idx
                        && p1_pos <= end_idx
                    {
                        chunk_peaks.push((p1_pos, p1_sc));
                    }

                    p2_sc = p1_sc;
                    p1_pos = curr_pos;
                    p1_sc = curr_sc;
                }

                if p1_sc >= 35.0 && p1_sc >= p2_sc && p1_pos > start_idx && p1_pos <= end_idx {
                    chunk_peaks.push((p1_pos, p1_sc));
                }

                chunk_peaks
            })
            .collect();

        // Merge peaks within 60 bp
        let mut filtered_peaks: Vec<usize> = Vec::new();
        let mut last_p = 0usize;
        let mut last_s = 0.0f32;

        for (p, s) in raw_peaks {
            if !filtered_peaks.is_empty() && p.saturating_sub(last_p) <= 60 {
                if s > last_s {
                    *filtered_peaks.last_mut().unwrap() = p;
                    last_p = p;
                    last_s = s;
                }
            } else {
                filtered_peaks.push(p);
                last_p = p;
                last_s = s;
            }
        }

        filtered_peaks
    }
}

/// A DP cell tracking score and start position.
#[derive(Clone, Copy, Debug)]
struct Cell {
    score: f32,
    start: u16,
}

impl Default for Cell {
    fn default() -> Self {
        Cell {
            score: -99999.0,
            start: 0,
        }
    }
}

/// High-performance aligner with preallocated zero-copy buffers.
#[derive(Clone)]
pub struct Aligner {
    pub gap_open: f32,
    pub gap_ext: f32,
    prev_m: Vec<Cell>,
    prev_i: Vec<Cell>,
    prev_d: Vec<Cell>,
    curr_m: Vec<Cell>,
    curr_i: Vec<Cell>,
    curr_d: Vec<Cell>,
}

impl Default for Aligner {
    fn default() -> Self {
        let cap = 256;
        Aligner {
            gap_open: -5.0,
            gap_ext: -0.6,
            prev_m: vec![Cell::default(); cap],
            prev_i: vec![Cell::default(); cap],
            prev_d: vec![Cell::default(); cap],
            curr_m: vec![Cell::default(); cap],
            curr_i: vec![Cell::default(); cap],
            curr_d: vec![Cell::default(); cap],
        }
    }
}

#[derive(Debug, Clone)]
pub struct AlignmentResult {
    pub score: f32,
    pub hmm_from: usize,
    pub hmm_to: usize,
    pub ali_from: usize,
    pub ali_to: usize,
}

impl Aligner {
    /// Align target_slice against `model` using a narrow band anchored around `end_in_slice`.
    /// Reuses preallocated internal buffers without heap allocations.
    pub fn align(
        &mut self,
        model: &HmmModel,
        target_slice: &[usize],
        end_in_slice: usize,
    ) -> Option<AlignmentResult> {
        let m_len = model.length;
        let t_len = target_slice.len();
        if t_len < 100 || m_len == 0 {
            return None;
        }

        let needed = t_len + 1;
        if self.prev_m.len() < needed {
            let cap = needed + 64;
            self.prev_m.resize(cap, Cell::default());
            self.prev_i.resize(cap, Cell::default());
            self.prev_d.resize(cap, Cell::default());
            self.curr_m.resize(cap, Cell::default());
            self.curr_i.resize(cap, Cell::default());
            self.curr_d.resize(cap, Cell::default());
        }

        let expected_start = end_in_slice.saturating_sub(m_len);
        let start_min = expected_start.saturating_sub(14);
        let start_max = (expected_start + 28).min(t_len.saturating_sub(100));

        let clear_start = start_min.saturating_sub(2);
        let clear_end = (expected_start + m_len + 32).min(t_len + 1);
        for j in clear_start..clear_end {
            self.prev_m[j] = Cell::default();
            self.prev_i[j] = Cell::default();
            self.prev_d[j] = Cell::default();
        }

        for j in start_min..=start_max {
            self.prev_m[j] = Cell {
                score: 0.0,
                start: (j + 1) as u16,
            };
        }

        let mut best_score = -99999.0f32;
        let mut best_start = 1usize;
        let mut best_k = m_len;
        let mut best_j = m_len.min(t_len);

        for k in 1..=m_len {
            let j_center = expected_start + k;
            let j_lo = j_center.saturating_sub(14).max(1);
            let j_hi = (j_center + 28).min(t_len);

            for j in j_lo.saturating_sub(1)..=(j_hi + 1).min(t_len) {
                self.curr_m[j] = Cell::default();
                self.curr_i[j] = Cell::default();
                self.curr_d[j] = Cell::default();
            }

            for j in j_lo..=j_hi {
                let b = target_slice[j - 1];
                let match_sc = model.match_scores[k - 1][b];

                let m_prev = self.prev_m[j - 1];
                let i_prev = self.prev_i[j - 1];
                let d_prev = self.prev_d[j - 1];

                let mut best_prev = m_prev;
                if i_prev.score > best_prev.score {
                    best_prev = i_prev;
                }
                if d_prev.score > best_prev.score {
                    best_prev = d_prev;
                }

                if best_prev.score > -90000.0 {
                    self.curr_m[j] = Cell {
                        score: best_prev.score + match_sc,
                        start: best_prev.start,
                    };
                }

                let m_ins = self.curr_m[j - 1];
                let i_ins = self.curr_i[j - 1];
                let mut best_ins = Cell {
                    score: m_ins.score + self.gap_open,
                    start: m_ins.start,
                };
                let ext_sc = i_ins.score + self.gap_ext;
                if ext_sc > best_ins.score {
                    best_ins = Cell {
                        score: ext_sc,
                        start: i_ins.start,
                    };
                }
                if best_ins.score > -90000.0 {
                    self.curr_i[j] = best_ins;
                }

                let m_del = self.prev_m[j];
                let d_del = self.prev_d[j];
                let mut best_del = Cell {
                    score: m_del.score + self.gap_open,
                    start: m_del.start,
                };
                let ext_del = d_del.score + self.gap_ext;
                if ext_del > best_del.score {
                    best_del = Cell {
                        score: ext_del,
                        start: d_del.start,
                    };
                }
                if best_del.score > -90000.0 {
                    self.curr_d[j] = best_del;
                }

                if k >= m_len.saturating_sub(5) {
                    let sc = self.curr_m[j].score.max(self.curr_d[j].score);
                    let st = if self.curr_m[j].score >= self.curr_d[j].score {
                        self.curr_m[j].start as usize
                    } else {
                        self.curr_d[j].start as usize
                    };

                    if sc > best_score {
                        best_score = sc;
                        best_start = st;
                        best_k = k;
                        best_j = j;
                    }
                }
            }

            std::mem::swap(&mut self.prev_m, &mut self.curr_m);
            std::mem::swap(&mut self.prev_i, &mut self.curr_i);
            std::mem::swap(&mut self.prev_d, &mut self.curr_d);
        }

        if best_score <= 0.0 {
            return None;
        }

        Some(AlignmentResult {
            score: best_score,
            hmm_from: 1,
            hmm_to: best_k,
            ali_from: best_start,
            ali_to: best_j,
        })
    }
}

/// Scan a sequence against all models in `HmmCollection` completely in-memory with Rayon.
pub fn scan_sequence_in_memory(
    target_name: &str,
    seq: &[u8],
    coll: &HmmCollection,
    master: &MasterProfile,
    score_threshold: f64,
) -> Vec<BedRecord> {
    let mut hits: Vec<HmmHit> = Vec::new();

    // 1. Scan plus strand
    scan_strand(
        target_name,
        seq,
        '+',
        coll,
        master,
        score_threshold,
        &mut hits,
    );

    // 2. Scan minus strand
    let rev_seq = reverse_complement(seq);
    scan_strand(
        target_name,
        &rev_seq,
        '-',
        coll,
        master,
        score_threshold,
        &mut hits,
    );

    // 3. Convert hits to BedRecords
    let mut records: Vec<BedRecord> = Vec::new();
    let seq_len = seq.len() as u64;

    for hit in hits {
        let (lo, hi) = if hit.strand == '+' {
            (hit.ali_from, hit.ali_to)
        } else {
            let start = seq_len.saturating_sub(hit.ali_to) + 1;
            let end = seq_len.saturating_sub(hit.ali_from) + 1;
            (start, end)
        };

        let start = match lo.checked_sub(1) {
            Some(s) => s,
            None => continue,
        };
        let end = hi;
        let ali_length = end.saturating_sub(start);
        if ali_length == 0 {
            continue;
        }

        let score_per_len = hit.score / ali_length as f64;
        let formatted: f64 = format!("{:.1}", score_per_len)
            .parse()
            .unwrap_or(score_per_len);
        if formatted < score_threshold {
            continue;
        }

        let color = get_color(&hit.query_name);

        records.push(BedRecord {
            chrom: hit.target_name.clone(),
            start,
            end,
            name: hit.query_name.clone(),
            score: hit.score,
            strand: hit.strand,
            thick_start: start,
            thick_end: end,
            color,
            hmm_from: hit.hmm_from,
            hmm_to: hit.hmm_to,
        });
    }

    // 4. Sort and filter overlaps
    crate::bed::sort_bed_records(&mut records);
    filter_overlaps(&records)
}

fn scan_strand(
    target_name: &str,
    strand_seq: &[u8],
    strand_char: char,
    coll: &HmmCollection,
    master: &MasterProfile,
    score_threshold: f64,
    out_hits: &mut Vec<HmmHit>,
) {
    let target_indices: Vec<usize> = strand_seq.iter().map(|&b| base_to_idx(b)).collect();
    let n = target_indices.len();
    if n < 100 {
        return;
    }

    // 1. Find monomer peaks using the master profile
    let ends = master.find_monomer_ends(&target_indices);
    let target_slice_ref = target_indices.as_slice();

    // 2. Parallel scan across windows with Rayon, reusing scratch buffers per chunk
    let hits: Vec<HmmHit> = ends
        .par_chunks(32)
        .flat_map_iter(|chunk| {
            let mut aligner = Aligner::default();

            chunk.iter().filter_map(move |&end_pos| {
                let start_pos = end_pos.saturating_sub(190);
                let end_idx = (end_pos + 5).min(n);
                let slice = &target_slice_ref[start_pos..end_idx];
                let end_in_slice = (end_pos - start_pos).min(slice.len());

                let mut best_hit: Option<HmmHit> = None;
                let mut best_score = 0.0f32;

                for model in &coll.models {
                    if let Some(res) = aligner.align(model, slice, end_in_slice) {
                        let ali_len = (res.ali_to - res.ali_from + 1) as f32;
                        let hmmer_approx_score = (res.score * 0.85) as f64;
                        let hmmer_sc_per_len = hmmer_approx_score / ali_len as f64;

                        if hmmer_sc_per_len >= score_threshold - 0.05 && res.score > best_score {
                            best_score = res.score;
                            let ali_from = (start_pos + res.ali_from) as u64;
                            let ali_to = (start_pos + res.ali_to) as u64;
                            best_hit = Some(HmmHit {
                                target_name: target_name.to_string(),
                                query_name: model.name.clone(),
                                hmm_from: res.hmm_from as u64,
                                hmm_to: res.hmm_to as u64,
                                ali_from,
                                ali_to,
                                strand: strand_char,
                                score: hmmer_approx_score,
                            });
                        }
                    }
                }

                best_hit
            })
        })
        .collect();

    out_hits.extend(hits);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ACGTN"), b"NACGT");
        assert_eq!(reverse_complement(b"AATCT"), b"AGATT");
    }

    #[test]
    fn test_base_to_idx() {
        assert_eq!(base_to_idx(b'A'), 0);
        assert_eq!(base_to_idx(b'C'), 1);
        assert_eq!(base_to_idx(b'G'), 2);
        assert_eq!(base_to_idx(b'T'), 3);
        assert_eq!(base_to_idx(b'a'), 0);
        assert_eq!(base_to_idx(b'N'), 0);
    }

    #[test]
    fn test_affine_sw_single_insertion() {
        // Reviewer's exact test case:
        // Profile: ACGT (len 4), match = +6, mismatch = -20
        // Target sequence: ACAGT (len 5, insertion of A between C and G)
        // Analytical score: 4 matches * (+6.0) + gap_open(-5.0) = 19.0 (a 1-bp insertion only incurs opening penalty)
        let l = 4;
        let mut pssm = vec![[-20.0f32; 4]; l];
        pssm[0][0] = 6.0; // A
        pssm[1][1] = 6.0; // C
        pssm[2][2] = 6.0; // G
        pssm[3][3] = 6.0; // T

        let master = MasterProfile { length: l, pssm };
        let target = b"ACAGT";
        let seq_indices: Vec<usize> = target.iter().map(|&b| base_to_idx(b)).collect();

        let gap_open = -5.0f32;
        let gap_ext = -0.6f32;

        // Directly exercises production MasterProfile::step_dp via align_sequence:
        let score = master.align_sequence(&seq_indices, gap_open, gap_ext);
        assert_eq!(
            score, 19.0f32,
            "Production DP recurrence with 1-bp insertion must yield exactly 19.0"
        );
    }

    #[test]
    fn test_find_monomer_ends_synthetic() {
        // Representative 171-bp alpha-satellite monomer consensus
        let raw_consensus = b"AATTTCAGCTGACTAAACAGTGATTTTTGTACTCTTTGCTCGAGTATTTTGGATCCCGTCTAGCTAACGCTTGTTTTGTGTGGGGTGTGAGCTTCGCTTCCCAATTCTTTATCCAGTATATTTGGACACCCATTCCAGCTACACATTATTTGCAGTGTGTATTTGGACACCTTT";
        let consensus_171 = &raw_consensus[..171];
        assert_eq!(consensus_171.len(), 171);

        let l = 171;
        let mut pssm = vec![[-10.0f32; 4]; l];
        for (i, &b) in consensus_171.iter().enumerate() {
            pssm[i][base_to_idx(b)] = 3.0;
        }

        let master = MasterProfile { length: l, pssm };

        // Construct a target sequence with two consecutive 171-bp monomers
        let mut target = Vec::with_capacity(342);
        target.extend_from_slice(consensus_171);
        target.extend_from_slice(consensus_171);
        let seq_indices: Vec<usize> = target.iter().map(|&b| base_to_idx(b)).collect();

        let peaks = master.find_monomer_ends(&seq_indices);
        assert_eq!(peaks.len(), 2, "Expected exactly two monomer end peaks");
        assert_eq!(peaks[0], 171, "First monomer peak must end at 171 bp");
        assert_eq!(peaks[1], 342, "Second monomer peak must end at 342 bp");
    }

    #[test]
    fn test_scan_test_fa_if_present() {
        let hmm_path = std::path::Path::new("data/AS-SFs-hmmer3.0.290621.hmm");
        let fa_path = std::path::Path::new("data/test/test.fa");
        if hmm_path.exists() && fa_path.exists() {
            let coll = HmmCollection::load_from_file(hmm_path).unwrap();
            let master = MasterProfile::from_collection(&coll);

            let fasta_content = std::fs::read_to_string(fa_path).unwrap();
            let mut seq = Vec::new();
            for line in fasta_content.lines() {
                if !line.starts_with('>') {
                    seq.extend_from_slice(line.trim().as_bytes());
                }
            }

            let start = std::time::Instant::now();
            let records = scan_sequence_in_memory("chr11", &seq, &coll, &master, 0.7);
            let elapsed = start.elapsed();

            println!("In-memory SF scan with Rayon finished in {:?}", elapsed);
            println!("Found {} SF records:", records.len());
            for r in &records {
                println!(
                    "  {}\t{}\t{}\t{}\t{:.1}\t{}",
                    r.chrom, r.start, r.end, r.name, r.score, r.strand
                );
            }
            assert!(!records.is_empty());
        }
    }

    #[test]
    fn test_scan_test_fa_hor_if_present() {
        let hmm_path = std::path::Path::new("data/AS-HORs-hmmer3.3.2-120124.hmm");
        let fa_path = std::path::Path::new("data/test/test.fa");
        if hmm_path.exists() && fa_path.exists() {
            let coll = HmmCollection::load_from_file(hmm_path).unwrap();
            let master = MasterProfile::from_collection(&coll);

            let fasta_content = std::fs::read_to_string(fa_path).unwrap();
            let mut seq = Vec::new();
            for line in fasta_content.lines() {
                if !line.starts_with('>') {
                    seq.extend_from_slice(line.trim().as_bytes());
                }
            }

            let start = std::time::Instant::now();
            let records = scan_sequence_in_memory("chr11", &seq, &coll, &master, 0.7);
            let elapsed = start.elapsed();

            println!(
                "In-memory HOR scan (1963 models) with Rayon finished in {:?}",
                elapsed
            );
            println!("Found {} HOR records:", records.len());
            for r in &records {
                println!(
                    "  {}\t{}\t{}\t{}\t{:.1}\t{}",
                    r.chrom, r.start, r.end, r.name, r.score, r.strand
                );
            }
            assert!(!records.is_empty());
        }
    }
}
