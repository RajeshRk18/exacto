pub mod keygen;
pub mod encrypt;
pub mod eval;
pub mod encoding;
pub mod keyswitch;
pub mod modswitch;

pub use keygen::{SecretKey, PublicKey, RelinKey, GaloisKey};
pub use encrypt::{encrypt_pk, encrypt_sk, decrypt};
pub use eval::{bfv_add, bfv_sub, bfv_neg, bfv_mul_and_relin, bfv_plain_mul, bfv_apply_automorphism};
pub use encoding::{encode_scalar, decode_scalar, encode_simd, decode_simd};

use std::sync::Arc;
use crate::params::BfvParams;
use crate::ring::rns::RnsPoly;

/// A BFV ciphertext: (c0, c1, ..., c_k) where k=1 for fresh, k=2 after mul (before relin).
#[derive(Clone, Debug)]
pub struct BfvCiphertext {
    /// Ciphertext components. Typically 2 (c0, c1) or 3 after multiplication.
    pub c: Vec<RnsPoly>,
    /// Associated parameters.
    pub params: Arc<BfvParams>,
}

impl BfvCiphertext {
    pub fn degree(&self) -> usize {
        self.c.len() - 1
    }
}
