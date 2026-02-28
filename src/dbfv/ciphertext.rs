use std::sync::Arc;
use crate::bfv::BfvCiphertext;
use crate::params::DbfvParams;

/// A dBFV ciphertext: d BFV ciphertexts, one per digit.
///
/// Represents an encryption of μ ∈ Z_p where μ = Σ μ_i · b^i.
/// Each limb encrypts μ_i ∈ Z_b under BFV with plaintext modulus b.
#[derive(Clone, Debug)]
pub struct DbfvCiphertext {
    /// The d BFV ciphertext "limbs", one per digit.
    pub limbs: Vec<BfvCiphertext>,
    /// Current polynomial degree in B (normally d, may be 2d-1 after mul before reduction).
    pub degree: usize,
    /// Number of multiplicative dBFV layers applied since fresh encryption.
    ///
    /// The current implementation supports one multiplication layer without
    /// ciphertext-level lattice reduction.
    pub mul_depth: usize,
    /// Parameters.
    pub params: Arc<DbfvParams>,
}

impl DbfvCiphertext {
    /// Number of limbs.
    pub fn num_limbs(&self) -> usize {
        self.limbs.len()
    }

    /// Check if this ciphertext needs degree reduction (degree > d).
    pub fn needs_reduction(&self) -> bool {
        self.degree > self.params.num_digits
    }
}
