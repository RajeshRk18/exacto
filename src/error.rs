use thiserror::Error;

#[derive(Debug, Error)]
pub enum ExactoError {
    #[error("invalid parameter: {0}")]
    InvalidParam(String),

    #[error("dimension mismatch: expected {expected}, got {got}")]
    DimensionMismatch { expected: usize, got: usize },

    #[error("modulus mismatch")]
    ModulusMismatch,

    #[error("ring degree must be a power of 2, got {0}")]
    InvalidRingDegree(usize),

    #[error("decryption error: noise budget exhausted")]
    DecryptionError,

    #[error("decomposition error: {0}")]
    DecompositionError(String),

    #[error("lattice error: {0}")]
    LatticeError(String),

    #[error("key not available: {0}")]
    MissingKey(String),

    #[error("not yet implemented: {0}")]
    NotImplemented(String),
}

pub type Result<T> = std::result::Result<T, ExactoError>;
