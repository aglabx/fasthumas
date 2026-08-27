use std::path::PathBuf;

#[derive(Debug, thiserror::Error)]
pub enum PipelineError {
    #[error("nhmmer failed for {path}: {message}")]
    NhmmerFailed { path: PathBuf, message: String },

    #[error("nhmmer not found in PATH. Install HMMER: http://hmmer.org/")]
    NhmmerNotFound,

    #[error("hmmemit not found in PATH. Install HMMER: http://hmmer.org/")]
    HmmemitNotFound,

    #[error("hmmemit failed for {path}: {message}")]
    HmmemitFailed { path: PathBuf, message: String },

    #[error("input FASTA not found, or is not a regular file: {0}")]
    InputNotFile(PathBuf),

    #[error("no sequences found in {0}")]
    NoSequences(PathBuf),

    #[error("HMM profile not found: {0}")]
    HmmNotFound(PathBuf),

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}
