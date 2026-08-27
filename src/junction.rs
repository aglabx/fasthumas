use std::collections::{BTreeSet, HashMap};

use crate::bed::BedRecord;
use crate::consensus::Consensus;

/// How far a junction gap may be closed, and what makes a returned base
/// believable.
#[derive(Debug, Clone, Copy)]
pub struct GapPolicy {
    pub max_gap: u64,
    /// Fraction of returned bases that must match the model consensus.
    pub min_identity: f64,
    /// How far the trimmed model positions may exceed the gap before the
    /// junction is left alone. The excess means the genome is missing bases the
    /// model expects, which is a deletion, not something to invent.
    pub deletion_tolerance: u64,
    /// How many times a base must recur at one model position before it counts
    /// as a variant of that monomer rather than noise.
    pub min_recurrence: usize,
    /// ... and what share of the observations at that position it must hold.
    pub min_dominance: f64,
}

/// The junction pass is deliberately not tunable: monomer boundaries are a
/// property of the models and the sequence, not a knob for the caller.
pub const POLICY: GapPolicy = GapPolicy {
    max_gap: 3,
    min_identity: 0.5,
    deletion_tolerance: 2,
    min_recurrence: 10,
    min_dominance: 0.5,
};

#[derive(Debug, Default, Clone, Copy)]
pub struct GapStats {
    /// Junctions with a gap of 1..=max_gap.
    pub candidates: usize,
    /// ... of those, gaps the trimmed model positions account for.
    pub accounted: usize,
    /// ... of those, skipped because duplicate model names left the split unclear.
    pub unresolved: usize,
    /// Sides of accounted junctions, each owed one or more bases.
    pub sides: usize,
    /// ... accepted because the bases match the model consensus.
    pub by_consensus: usize,
    /// ... accepted because the same variant recurs at that model position.
    pub by_recurrence: usize,
    /// ... left alone.
    pub rejected: usize,
    /// Bases handed back to a monomer.
    pub bases: u64,
}

impl GapStats {
    pub fn merge(&mut self, o: &GapStats) {
        self.candidates += o.candidates;
        self.accounted += o.accounted;
        self.unresolved += o.unresolved;
        self.sides += o.sides;
        self.by_consensus += o.by_consensus;
        self.by_recurrence += o.by_recurrence;
        self.rejected += o.rejected;
        self.bases += o.bases;
    }
}

/// A junction whose gap the trimmed model positions account for.
#[derive(Debug, Clone, Copy)]
pub struct Candidate {
    left: usize,
    right: usize,
    pub start: u64,
    pub end: u64,
    left_takes: u64,
    right_takes: u64,
}

/// One side of one junction: bases owed to a single record, and how well they
/// agree with that model.
#[derive(Debug, Clone)]
pub struct Side {
    /// Index of the record in the track this side belongs to.
    pub record: usize,
    /// Extend the record's right edge, otherwise its left.
    pub extend_right: bool,
    pub bases: u64,
    pub model: String,
    /// First trimmed model position.
    pub model_pos: u64,
    /// The bases read in the model's own orientation, which is what makes a hit
    /// and its reverse-complement twin count towards the same variant.
    pub variant: String,
    /// Best agreement with the model consensus at those positions.
    pub identity: f64,
}

/// Find the junctions worth trying to close.
///
/// A gap between two adjacent monomers is not arbitrary: HMMER builds its
/// alignment by maximum expected accuracy and drops residues whose posterior is
/// low, which at a tandem junction is exactly where the register is ambiguous.
/// The tblout says how much was dropped — the alignment stops short of the
/// model's first or last position — so the gap can be attributed base by base
/// instead of guessed at.
///
/// One end of every hit is free of the model's length: `hmm_from - 1` bases
/// were dropped there whatever the model. The other end needs the length, and
/// where a name covers models of several lengths every one of them stays a
/// candidate.
///
/// The trimmed positions usually add up to the gap exactly. When they add up to
/// *more*, the genome is short of bases the model expects — a deletion — and
/// the gap can still be closed with the bases that are there, as long as only
/// one side is missing anything. Adding up to less is never enough to act on.
///
/// `records` must be sorted by sequence then start.
pub fn candidates(
    records: &[BedRecord],
    consensus: &Consensus,
    policy: &GapPolicy,
) -> (Vec<Candidate>, GapStats) {
    let mut out = Vec::new();
    let mut stats = GapStats::default();

    for i in 0..records.len().saturating_sub(1) {
        let (left, right) = (&records[i], &records[i + 1]);
        if left.chrom != right.chrom || right.start <= left.end {
            continue;
        }
        let gap = right.start - left.end;
        if gap > policy.max_gap {
            continue;
        }
        stats.candidates += 1;

        let lefts = trim_right_options(left, consensus);
        let rights = trim_left_options(right, consensus);
        let pairs: Vec<(u64, u64)> = lefts
            .iter()
            .flat_map(|a| rights.iter().map(move |b| (*a, *b)))
            .collect();

        let exact: Vec<(u64, u64)> = pairs
            .iter()
            .copied()
            .filter(|(a, b)| a + b == gap)
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();

        let takes = if exact.len() == 1 {
            exact[0]
        } else if !exact.is_empty() {
            stats.unresolved += 1;
            continue;
        } else {
            // Nothing adds up exactly. A deletion shows as more trimmed model
            // positions than there are bases; act only when a single side is
            // short, so there is no choice to make about who gets what.
            let over: Vec<(u64, u64)> = pairs
                .iter()
                .copied()
                .filter(|(a, b)| a + b > gap && a + b <= gap + policy.deletion_tolerance)
                .collect();
            if over.is_empty() {
                continue;
            }
            if over.iter().all(|(a, _)| *a == 0) {
                (0, gap)
            } else if over.iter().all(|(_, b)| *b == 0) {
                (gap, 0)
            } else {
                stats.unresolved += 1;
                continue;
            }
        };

        stats.accounted += 1;
        out.push(Candidate {
            left: i,
            right: i + 1,
            start: left.end,
            end: right.start,
            left_takes: takes.0,
            right_takes: takes.1,
        });
    }

    (out, stats)
}

/// Describe each side of each candidate, without changing anything yet: the
/// decision needs counts gathered over the whole run.
///
/// `bases[k]` must hold the target sequence over `candidates[k].start..end`.
pub fn sides(
    records: &[BedRecord],
    candidates: &[Candidate],
    bases: &[Vec<u8>],
    consensus: &Consensus,
) -> Vec<Side> {
    let mut out = Vec::new();

    for (c, observed) in candidates.iter().zip(bases) {
        if observed.len() as u64 != c.end - c.start {
            continue;
        }
        let observed: Vec<u8> = observed.iter().map(|b| b.to_ascii_uppercase()).collect();
        let split = c.left_takes as usize;

        if c.left_takes > 0 {
            let r = &records[c.left];
            let expected = expected_right(r, c.left_takes, consensus);
            out.push(side(
                c.left,
                r,
                true,
                c.left_takes,
                trimmed_position_right(r),
                &observed[..split],
                &expected,
            ));
        }
        if c.right_takes > 0 {
            let r = &records[c.right];
            let expected = expected_left(r, c.right_takes, consensus);
            out.push(side(
                c.right,
                r,
                false,
                c.right_takes,
                trimmed_position_left(r),
                &observed[split..],
                &expected,
            ));
        }
    }

    out
}

#[allow(clippy::too_many_arguments)]
fn side(
    record: usize,
    r: &BedRecord,
    extend_right: bool,
    bases: u64,
    model_pos: u64,
    observed: &[u8],
    expected: &[Vec<u8>],
) -> Side {
    let variant: String = if r.strand == '-' {
        observed
            .iter()
            .rev()
            .map(|b| complement(*b) as char)
            .collect()
    } else {
        observed.iter().map(|b| *b as char).collect()
    };
    Side {
        record,
        extend_right,
        bases,
        model: r.name.clone(),
        model_pos,
        variant,
        identity: best_identity(observed, expected),
    }
}

/// How often each base has been seen at each model position, over the whole
/// run. Bases are counted in the model's own orientation, so the two strands
/// pool together.
#[derive(Debug, Default)]
pub struct Tally {
    seen: HashMap<(String, u64, String), usize>,
    total: HashMap<(String, u64), usize>,
}

impl Tally {
    pub fn add(&mut self, s: &Side) {
        *self
            .seen
            .entry((s.model.clone(), s.model_pos, s.variant.clone()))
            .or_default() += 1;
        *self
            .total
            .entry((s.model.clone(), s.model_pos))
            .or_default() += 1;
    }

    fn counts(&self, s: &Side) -> (usize, usize) {
        let key = (s.model.clone(), s.model_pos, s.variant.clone());
        (
            self.seen.get(&key).copied().unwrap_or(0),
            self.total
                .get(&(s.model.clone(), s.model_pos))
                .copied()
                .unwrap_or(0),
        )
    }

    /// A base is believable if the model expects it, or if the array itself has
    /// shown the same base at the same model position often enough and
    /// consistently enough to call it a variant of that monomer.
    pub fn verdict(&self, s: &Side, policy: &GapPolicy) -> Verdict {
        if s.identity >= policy.min_identity {
            return Verdict::Consensus;
        }
        let (seen, total) = self.counts(s);
        if seen >= policy.min_recurrence && seen as f64 >= policy.min_dominance * total as f64 {
            return Verdict::Recurrent { seen, total };
        }
        Verdict::Rejected
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Verdict {
    Consensus,
    Recurrent { seen: usize, total: usize },
    Rejected,
}

impl Verdict {
    pub fn accepted(&self) -> bool {
        !matches!(self, Verdict::Rejected)
    }
}

/// Model positions trimmed off the target-right end of a hit.
fn trim_right_options(r: &BedRecord, consensus: &Consensus) -> BTreeSet<u64> {
    if r.strand == '-' {
        BTreeSet::from([r.hmm_from.saturating_sub(1)])
    } else {
        consensus
            .lengths(&r.name)
            .into_iter()
            .map(|len| len.saturating_sub(r.hmm_to))
            .collect()
    }
}

/// Model positions trimmed off the target-left end of a hit.
fn trim_left_options(r: &BedRecord, consensus: &Consensus) -> BTreeSet<u64> {
    if r.strand == '-' {
        consensus
            .lengths(&r.name)
            .into_iter()
            .map(|len| len.saturating_sub(r.hmm_to))
            .collect()
    } else {
        BTreeSet::from([r.hmm_from.saturating_sub(1)])
    }
}

/// The first model position past the target-right end of a hit.
fn trimmed_position_right(r: &BedRecord) -> u64 {
    if r.strand == '-' {
        r.hmm_from.saturating_sub(1)
    } else {
        r.hmm_to + 1
    }
}

/// The first model position past the target-left end of a hit.
fn trimmed_position_left(r: &BedRecord) -> u64 {
    if r.strand == '-' {
        r.hmm_to + 1
    } else {
        r.hmm_from.saturating_sub(1)
    }
}

/// The `k` bases the model expects immediately to the right of `r`, in target
/// left-to-right order — one reading per model that could be behind the name.
fn expected_right(r: &BedRecord, k: u64, consensus: &Consensus) -> Vec<Vec<u8>> {
    consensus
        .variants(&r.name)
        .iter()
        .filter_map(|cons| {
            if r.strand != '-' && (cons.len() as u64) < r.hmm_to + k {
                return None;
            }
            (1..=k)
                .map(|j| {
                    Some(if r.strand == '-' {
                        complement(*cons.get(model_index(r.hmm_from.checked_sub(j)?)?)?)
                    } else {
                        *cons.get(model_index(r.hmm_to + j)?)?
                    })
                })
                .collect()
        })
        .collect()
}

/// The `k` bases the model expects immediately to the left of `r`.
fn expected_left(r: &BedRecord, k: u64, consensus: &Consensus) -> Vec<Vec<u8>> {
    consensus
        .variants(&r.name)
        .iter()
        .filter_map(|cons| {
            if r.strand == '-' && (cons.len() as u64) < r.hmm_to + k {
                return None;
            }
            (1..=k)
                .rev()
                .map(|j| {
                    Some(if r.strand == '-' {
                        complement(*cons.get(model_index(r.hmm_to + j)?)?)
                    } else {
                        *cons.get(model_index(r.hmm_from.checked_sub(j)?)?)?
                    })
                })
                .collect()
        })
        .collect()
}

fn model_index(position: u64) -> Option<usize> {
    if position == 0 {
        None
    } else {
        Some((position - 1) as usize)
    }
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other,
    }
}

/// The best agreement across the readings the name allows: a name standing for
/// several models is given the benefit of the doubt, since the hit could have
/// come from any of them.
fn best_identity(observed: &[u8], expected: &[Vec<u8>]) -> f64 {
    expected
        .iter()
        .map(|e| identity(observed, e))
        .fold(0.0, f64::max)
}

fn identity(observed: &[u8], expected: &[u8]) -> f64 {
    if observed.is_empty() || observed.len() != expected.len() {
        return 0.0;
    }
    let hits = observed
        .iter()
        .zip(expected)
        .filter(|(o, e)| o == e && **o != b'N')
        .count();
    hits as f64 / observed.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bed::test_record;

    fn table() -> Consensus {
        Consensus::parse(b">m-consensus\nACGTACGTAC\n>n-consensus\nTTTTGGGGCC\n")
    }

    /// Two plus-strand monomers, each trimmed by one base, with a 2 bp gap.
    fn pair() -> Vec<BedRecord> {
        let mut left = test_record("chr1", 0, 9, "m", 50.0);
        left.hmm_from = 1;
        left.hmm_to = 9;
        let mut right = test_record("chr1", 11, 20, "n", 50.0);
        right.hmm_from = 2;
        right.hmm_to = 10;
        vec![left, right]
    }

    fn run(recs: &mut [BedRecord], gap: &[u8], policy: &GapPolicy) -> (Vec<Side>, Tally) {
        let (cands, _) = candidates(recs, &table(), policy);
        let sides = sides(recs, &cands, &[gap.to_vec()], &table());
        let mut tally = Tally::default();
        for s in &sides {
            tally.add(s);
        }
        (sides, tally)
    }

    #[test]
    fn test_split_by_trimming() {
        let recs = pair();
        let (cands, stats) = candidates(&recs, &table(), &POLICY);
        assert_eq!((stats.candidates, stats.accounted), (1, 1));
        assert_eq!((cands[0].left_takes, cands[0].right_takes), (1, 1));
    }

    #[test]
    fn test_consensus_match_is_accepted() {
        let mut recs = pair();
        // position 10 of "m" is C; position 1 of "n" is T
        let (sides, tally) = run(&mut recs, b"CT", &POLICY);
        assert_eq!(sides.len(), 2);
        for s in &sides {
            assert_eq!(tally.verdict(s, &POLICY), Verdict::Consensus);
        }
    }

    #[test]
    fn test_lone_mismatch_is_rejected() {
        let mut recs = pair();
        let (sides, tally) = run(&mut recs, b"GA", &POLICY);
        for s in &sides {
            assert_eq!(tally.verdict(s, &POLICY), Verdict::Rejected);
        }
    }

    /// The same mismatching base seen often enough at the same model position
    /// stops being noise and counts as a variant of that monomer.
    #[test]
    fn test_recurring_variant_is_accepted() {
        let mut recs = pair();
        let (sides, _) = run(&mut recs, b"GA", &POLICY);
        let left = sides[0].clone();
        let mut tally = Tally::default();
        for _ in 0..12 {
            tally.add(&left);
        }
        assert!(matches!(
            tally.verdict(&left, &POLICY),
            Verdict::Recurrent {
                seen: 12,
                total: 12
            }
        ));
        // ... but not when it is a minority at that position
        let mut other = left.clone();
        other.variant = "C".into();
        let mut mixed = Tally::default();
        for _ in 0..12 {
            mixed.add(&left);
        }
        for _ in 0..40 {
            mixed.add(&other);
        }
        assert_eq!(mixed.verdict(&left, &POLICY), Verdict::Rejected);
    }

    /// Bases are pooled in model orientation, so a minus-strand observation
    /// counts towards the same variant as its plus-strand twin.
    #[test]
    fn test_variant_is_recorded_in_model_orientation() {
        let mut recs = pair();
        recs[0].strand = '-';
        recs[0].hmm_from = 2; // one position short at the target-right end
        recs[0].hmm_to = 10;
        recs[1].hmm_from = 1; // the right record is complete
        recs[1].start = 10; // ... leaving a single base in the gap
        let (cands, _) = candidates(&recs, &table(), &POLICY);
        let sides = sides(&recs, &cands, &[b"T".to_vec()], &table());
        assert_eq!(sides.len(), 1);
        assert_eq!(sides[0].variant, "A"); // complemented back into the model
    }

    /// More trimmed positions than there are bases means the genome is missing
    /// one: still closable, as long as only one side is short.
    #[test]
    fn test_deletion_closes_from_the_short_side() {
        let mut recs = pair();
        recs[0].hmm_to = 10; // left record is complete
        recs[1].hmm_from = 3; // right record is missing two positions
        recs[1].start = 10; // ... but only one base sits in the gap
        let (cands, stats) = candidates(&recs, &table(), &POLICY);
        assert_eq!(stats.accounted, 1);
        assert_eq!((cands[0].left_takes, cands[0].right_takes), (0, 1));
    }

    /// If both sides are short there is no way to say who gets the bases.
    #[test]
    fn test_deletion_with_both_sides_short_is_skipped() {
        let mut recs = pair();
        recs[0].hmm_to = 9; // one missing on the left
        recs[1].hmm_from = 3; // two missing on the right
        recs[1].start = 10; // one base in the gap
        let (cands, stats) = candidates(&recs, &table(), &POLICY);
        assert_eq!((stats.accounted, stats.unresolved), (0, 1));
        assert!(cands.is_empty());
    }
}
