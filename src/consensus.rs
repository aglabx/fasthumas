use std::collections::{BTreeSet, HashMap};
use std::path::Path;

use crate::error::PipelineError;
use crate::hmm::HmmCollection;

/// Majority-rule consensus sequences for the models in a profile, as produced
/// by `hmmemit -c`. One residue per match state, so a consensus is as long as
/// its model.
///
/// A profile may carry several models under one name — the HOR profile has
/// 1923 models but only 1238 distinct names, and 58 of those names cover
/// models of different lengths. A tblout hit names its model and nothing else,
/// so every model sharing that name stays a candidate; the junction pass has
/// to consider all of them and only acts when they agree on an answer.
#[derive(Debug, Default)]
pub struct Consensus {
    models: HashMap<String, Vec<Vec<u8>>>,
}

impl Consensus {
    #[allow(dead_code)]
    pub fn from_profile(hmm: &Path) -> Result<Self, PipelineError> {
        let coll = HmmCollection::load_from_file(hmm)?;
        Ok(Self::from_collection(&coll))
    }

    pub fn from_collection(coll: &HmmCollection) -> Self {
        Consensus {
            models: coll.to_consensus_map(),
        }
    }

    /// Parse the FASTA `hmmemit -c` writes. Headers look like `>Ja-consensus`.
    #[allow(dead_code)]
    pub(crate) fn parse(fasta: &[u8]) -> Self {
        let mut models: HashMap<String, Vec<Vec<u8>>> = HashMap::new();
        let mut name: Option<String> = None;
        let mut seq: Vec<u8> = Vec::new();

        for line in fasta.split(|b| *b == b'\n') {
            if line.first() == Some(&b'>') {
                if let Some(previous) = name.take() {
                    add(&mut models, previous, std::mem::take(&mut seq));
                }
                let header = String::from_utf8_lossy(&line[1..]);
                let first = header.split_whitespace().next().unwrap_or("");
                name = Some(
                    first
                        .strip_suffix("-consensus")
                        .unwrap_or(first)
                        .to_string(),
                );
            } else if name.is_some() {
                seq.extend(
                    line.iter()
                        .filter(|b| b.is_ascii_alphabetic())
                        .map(|b| b.to_ascii_uppercase()),
                );
            }
        }
        if let Some(last) = name.take() {
            add(&mut models, last, seq);
        }

        Consensus { models }
    }

    /// Every distinct consensus filed under `name`.
    pub fn variants(&self, name: &str) -> &[Vec<u8>] {
        self.models.get(name).map_or(&[], |v| v.as_slice())
    }

    /// The distinct model lengths `name` could stand for.
    pub fn lengths(&self, name: &str) -> BTreeSet<u64> {
        self.variants(name).iter().map(|c| c.len() as u64).collect()
    }

    /// Distinct model names in the profile.
    pub fn names(&self) -> usize {
        self.models.len()
    }

    /// Names that stand for models of more than one length, where the model's
    /// length cannot be recovered from a hit.
    #[allow(dead_code)]
    pub fn length_ambiguous(&self) -> usize {
        self.models
            .values()
            .filter(|variants| {
                let lengths: BTreeSet<usize> = variants.iter().map(Vec::len).collect();
                lengths.len() > 1
            })
            .count()
    }
}

/// File a model under its name, keeping distinct consensuses apart.
#[allow(dead_code)]
fn add(models: &mut HashMap<String, Vec<Vec<u8>>>, name: String, seq: Vec<u8>) {
    let variants = models.entry(name).or_default();
    if !variants.contains(&seq) {
        variants.push(seq);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_hmmemit_output() {
        let out = b">Qa-consensus\nAATCTAAGT\naatg\n>Pa-consensus\nAAGCTCAG\n";
        let c = Consensus::parse(out);
        assert_eq!(c.names(), 2);
        assert_eq!(c.variants("Qa"), &[b"AATCTAAGTAATG".to_vec()]);
        assert_eq!(c.variants("Pa"), &[b"AAGCTCAG".to_vec()]);
        assert!(c.variants("Qa-consensus").is_empty());
        assert_eq!(c.length_ambiguous(), 0);
    }

    /// One name over models of different lengths: both stay candidates, and
    /// the name counts as length-ambiguous.
    #[test]
    fn test_duplicate_names_keep_every_variant() {
        let out = b">Qa-consensus\nAAAA\n>Qa-consensus\nAAAAA\n>Pa-consensus\nCCCC\n";
        let c = Consensus::parse(out);
        assert_eq!((c.names(), c.length_ambiguous()), (2, 1));
        assert_eq!(c.variants("Qa").len(), 2);
        assert_eq!(c.lengths("Qa"), BTreeSet::from([4, 5]));
        assert_eq!(c.lengths("Pa"), BTreeSet::from([4]));
    }

    /// Same length, different sequence: one length, two candidate readings.
    #[test]
    fn test_same_length_variants_are_not_length_ambiguous() {
        let c = Consensus::parse(b">Qa-consensus\nAAAA\n>Qa-consensus\nACAA\n");
        assert_eq!((c.names(), c.length_ambiguous()), (1, 0));
        assert_eq!(c.variants("Qa").len(), 2);
    }

    #[test]
    fn test_identical_duplicate_is_folded() {
        let c = Consensus::parse(b">Qa-consensus\nAAAA\n>Qa-consensus\nAAAA\n");
        assert_eq!(c.variants("Qa").len(), 1);
    }
}
