use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::error::PipelineError;

/// A single nucleotide HMM profile loaded into memory.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct HmmModel {
    pub name: String,
    pub length: usize,
    pub consensus: Vec<u8>,
    /// Match emissions log-odds in bits: [A, C, G, T] for each state 0..length.
    /// score = log2( P(residue | Mk) / P_null(residue) )
    pub match_scores: Vec<[f32; 4]>,
    /// Raw negative natural log probabilities: [A, C, G, T]
    pub match_neglog: Vec<[f32; 4]>,
    /// Transitions in nats or bits: [m->m, m->i, m->d, i->m, i->i, d->m, d->d]
    pub transitions: Vec<[f32; 7]>,
    /// Background composition in nats: [A, C, G, T]
    pub compo: [f32; 4],
}

impl HmmModel {
    /// Consensus sequence as ASCII string.
    #[allow(dead_code)]
    pub fn consensus_str(&self) -> String {
        String::from_utf8_lossy(&self.consensus).to_string()
    }
}

/// A collection of HMM models parsed from a `.hmm` file.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct HmmCollection {
    pub models: Vec<HmmModel>,
    /// Index of models by name. Note that multiple models may share the same name!
    pub by_name: HashMap<String, Vec<usize>>,
}

impl HmmCollection {
    /// Parse all models from a `.hmm` file into memory.
    pub fn load_from_file(path: &Path) -> Result<Self, PipelineError> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        Self::parse(reader)
    }

    /// Parse models from any `BufRead` stream.
    pub fn parse<R: BufRead>(mut reader: R) -> Result<Self, PipelineError> {
        let mut models = Vec::new();
        let mut by_name: HashMap<String, Vec<usize>> = HashMap::new();

        let mut line = String::new();
        let mut in_model = false;

        let mut name = String::new();
        let mut length = 0usize;
        let mut compo = [1.38629f32; 4]; // Default uniform 0.25 -> -ln(0.25) = 1.38629
        let mut match_neglog: Vec<[f32; 4]> = Vec::new();
        let mut transitions: Vec<[f32; 7]> = Vec::new();

        // ln(2) constant for converting nats to bits
        let ln2 = std::f32::consts::LN_2;

        while reader.read_line(&mut line)? > 0 {
            let trimmed = line.trim();

            if trimmed.starts_with("HMMER3") {
                in_model = true;
                name.clear();
                length = 0;
                match_neglog.clear();
                transitions.clear();
            } else if in_model {
                if let Some(rest) = trimmed.strip_prefix("NAME ") {
                    name = rest.trim().to_string();
                } else if let Some(rest) = trimmed.strip_prefix("LENG ") {
                    if let Ok(l) = rest.trim().parse::<usize>() {
                        length = l;
                        match_neglog.reserve(length);
                        transitions.reserve(length);
                    }
                } else if trimmed.starts_with("COMPO ") {
                    let parts: Vec<&str> = trimmed.split_whitespace().collect();
                    if parts.len() >= 5 {
                        for i in 0..4 {
                            compo[i] = parts[i + 1].parse::<f32>().unwrap_or(1.38629);
                        }
                    }
                } else if trimmed == "//" {
                    // End of current model
                    if !name.is_empty() && length > 0 && match_neglog.len() == length {
                        let mut consensus = Vec::with_capacity(length);
                        let mut match_scores = Vec::with_capacity(length);

                        for &nlogs in &match_neglog {
                            // argmin of negative log gives residue with highest probability
                            let mut min_val = f32::INFINITY;
                            let mut best_idx = 0;
                            for (idx, &v) in nlogs.iter().enumerate() {
                                if v < min_val {
                                    min_val = v;
                                    best_idx = idx;
                                }
                            }
                            let base = match best_idx {
                                0 => b'A',
                                1 => b'C',
                                2 => b'G',
                                _ => b'T',
                            };
                            consensus.push(base);

                            // Compute log-odds score in bits: (-val(k, base) + compo(base)) / ln(2)
                            let mut scores = [0.0f32; 4];
                            for b in 0..4 {
                                if nlogs[b] >= 900.0 {
                                    scores[b] = -50.0; // effectively zero probability
                                } else {
                                    scores[b] = (-nlogs[b] + compo[b]) / ln2;
                                }
                            }
                            match_scores.push(scores);
                        }

                        let idx = models.len();
                        by_name.entry(name.clone()).or_default().push(idx);

                        models.push(HmmModel {
                            name: std::mem::take(&mut name),
                            length,
                            consensus,
                            match_scores,
                            match_neglog: std::mem::take(&mut match_neglog),
                            transitions: std::mem::take(&mut transitions),
                            compo,
                        });
                    }

                    in_model = false;
                } else {
                    // Check if line is a match state line: starts with state number 1..length
                    // Format: <state_num> <p_A> <p_C> <p_G> <p_T> ...
                    let parts: Vec<&str> = trimmed.split_whitespace().collect();
                    if !parts.is_empty() {
                        if let Ok(state_num) = parts[0].parse::<usize>() {
                            if state_num == match_neglog.len() + 1 && parts.len() >= 5 {
                                let mut nlogs = [999.0f32; 4];
                                for i in 0..4 {
                                    if parts[i + 1] != "*" {
                                        nlogs[i] = parts[i + 1].parse::<f32>().unwrap_or(999.0);
                                    }
                                }
                                match_neglog.push(nlogs);

                                // The next line is insert emissions (skip)
                                line.clear();
                                reader.read_line(&mut line)?;

                                // The line after is transitions: m->m, m->i, m->d, i->m, i->i, d->m, d->d
                                line.clear();
                                reader.read_line(&mut line)?;
                                let trans_parts: Vec<&str> = line.split_whitespace().collect();
                                let mut tr = [999.0f32; 7];
                                for i in 0..7.min(trans_parts.len()) {
                                    if trans_parts[i] != "*" {
                                        tr[i] = trans_parts[i].parse::<f32>().unwrap_or(999.0);
                                    }
                                }
                                transitions.push(tr);
                            }
                        }
                    }
                }
            }

            line.clear();
        }

        Ok(HmmCollection { models, by_name })
    }

    /// Extract consensus mapping for `Consensus` struct directly, without calling hmmemit.
    pub fn to_consensus_map(&self) -> HashMap<String, Vec<Vec<u8>>> {
        let mut map: HashMap<String, Vec<Vec<u8>>> = HashMap::new();
        for m in &self.models {
            let list = map.entry(m.name.clone()).or_default();
            if !list.iter().any(|existing| existing == &m.consensus) {
                list.push(m.consensus.clone());
            }
        }
        map
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TINY_HMM: &str = r#"HMMER3/b [3.0 | March 2010]
NAME  TestModel
LENG  3
ALPH  DNA
HMM          A        C        G        T   
            m->m     m->i     m->d     i->m     i->i     d->m     d->d
  COMPO   1.38629  1.38629  1.38629  1.38629
          1.38629  1.38629  1.38629  1.38629
          0.04635  3.78772  3.78772  1.46634  0.26236  0.00000        *
      1   0.10000  2.50000  3.00000  4.00000      1 - -
          1.38629  1.38629  1.38629  1.38629
          0.04635  3.78772  3.78772  1.46634  0.26236  1.09861  0.40547
      2   3.00000  0.10000  3.00000  4.00000      2 - -
          1.38629  1.38629  1.38629  1.38629
          0.04635  3.78772  3.78772  1.46634  0.26236  1.09861  0.40547
      3   3.00000  3.00000  0.10000  4.00000      3 - -
          1.38629  1.38629  1.38629  1.38629
          0.04635  3.78772  3.78772  1.46634  0.26236  1.09861  0.40547
//
"#;

    #[test]
    fn test_parse_tiny_hmm() {
        let coll = HmmCollection::parse(TINY_HMM.as_bytes()).unwrap();
        assert_eq!(coll.models.len(), 1);
        let m = &coll.models[0];
        assert_eq!(m.name, "TestModel");
        assert_eq!(m.length, 3);
        assert_eq!(m.consensus, b"ACG");
        assert_eq!(m.consensus_str(), "ACG");
        assert_eq!(m.match_scores.len(), 3);
        // Position 1: A has lowest nlog (0.10000), so highest score
        assert!(m.match_scores[0][0] > m.match_scores[0][1]);
        assert!(m.match_scores[0][0] > 1.5);
    }

    #[test]
    fn test_parse_real_files_if_present() {
        let sf_path = Path::new("data/AS-SFs-hmmer3.0.290621.hmm");
        if sf_path.exists() {
            let start = std::time::Instant::now();
            let coll = HmmCollection::load_from_file(sf_path).unwrap();
            let elapsed = start.elapsed();
            println!("Loaded {} SF models in {:?}", coll.models.len(), elapsed);
            assert_eq!(coll.models.len(), 39);
        }

        let hor_path = Path::new("data/AS-HORs-hmmer3.3.2-120124.hmm");
        if hor_path.exists() {
            let start = std::time::Instant::now();
            let coll = HmmCollection::load_from_file(hor_path).unwrap();
            let elapsed = start.elapsed();
            println!("Loaded {} HOR models in {:?}", coll.models.len(), elapsed);
            assert_eq!(coll.models.len(), 1963);
        }
    }
}
