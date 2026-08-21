use std::path::PathBuf;

#[derive(Debug, thiserror::Error)]
pub enum PipelineError {
    #[error("nhmmer failed for {path}: {message}")]
    NhmmerFailed { path: PathBuf, message: String },

    #[error("nhmmer not found in PATH. Install HMMER: http://hmmer.org/")]
    NhmmerNotFound,

    #[error("no .fa files found in {0}")]
    NoInputFiles(PathBuf),

    #[error("HMM profile not found: {0}")]
    HmmNotFound(PathBuf),

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}
