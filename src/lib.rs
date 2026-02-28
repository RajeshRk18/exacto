//! # Exacto: dBFV High-Precision FHE Library
//!
//! Implements the "Decomposed BFV" (dBFV) scheme from Peikert-Zarchy-Zyskind (2026).
//!
//! Core idea: represent Z_p plaintexts as d digits in base b ≈ p^{1/d},
//! encrypt each digit as a BFV ciphertext. Error growth is proportional to
//! d·p^{1/d} instead of p, enabling high-precision homomorphic computation.
//!
//! ## Quick Start
//!
//! ```no_run
//! use exacto::prelude::*;
//!
//! // Create parameters: p=256, base=16, d=2
//! let params = exacto::params::presets::compact_dbfv().unwrap();
//!
//! // Generate keys
//! let (sk, pk, rlk) = exacto::dbfv::keygen::dbfv_keygen(&params).unwrap();
//!
//! // Encrypt
//! let ct1 = exacto::dbfv::dbfv_encrypt(10, &pk, &params).unwrap();
//! let ct2 = exacto::dbfv::dbfv_encrypt(20, &pk, &params).unwrap();
//!
//! // Add homomorphically
//! let ct_sum = exacto::dbfv::dbfv_add(&ct1, &ct2).unwrap();
//!
//! // Decrypt
//! let result = exacto::dbfv::dbfv_decrypt(&ct_sum, &sk).unwrap();
//! assert_eq!(result, 30);
//! ```

pub mod error;
pub mod params;
pub mod ring;
pub mod sampling;
pub mod bfv;
pub mod dbfv;
pub mod bootstrap;

/// Convenient re-exports for common types and functions.
pub mod prelude {
    pub use crate::error::{ExactoError, Result};
    pub use crate::params::{BfvParams, BfvParamsBuilder, DbfvParams};
    pub use crate::ring::{CoeffPoly, NttPoly, RnsPoly};
    pub use crate::bfv::{
        BfvCiphertext, SecretKey, PublicKey, RelinKey, GaloisKey,
        encrypt_pk, encrypt_sk, decrypt,
        bfv_add, bfv_sub, bfv_neg, bfv_mul_and_relin, bfv_plain_mul,
        encode_scalar, decode_scalar, encode_simd, decode_simd,
    };
    pub use crate::dbfv::{
        DbfvCiphertext,
        dbfv_encrypt, dbfv_encrypt_sk, dbfv_encrypt_poly, dbfv_encrypt_poly_sk,
        dbfv_decrypt, dbfv_decrypt_poly,
        dbfv_add, dbfv_sub, dbfv_neg, dbfv_mul,
        dbfv_apply_automorphism, dbfv_div_by_base, dbfv_change_base,
    };
    pub use crate::bootstrap::{
        BootstrapKey, gen_bootstrap_key, bfv_bootstrap, dbfv_bootstrap, dbfv_mul_then_bootstrap,
        dbfv_mul_chain_then_bootstrap,
    };
}
