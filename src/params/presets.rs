use std::sync::Arc;
use crate::error::Result;
use crate::params::{BfvParams, BfvParamsBuilder, DbfvParams};

/// NTT-friendly primes: q ≡ 1 (mod 2n).
/// These are chosen to be < 62 bits to allow u64 multiplication without overflow in Barrett.

/// Primes ≡ 1 (mod 2·4096 = 8192)
pub const PRIMES_4096: &[u64] = &[
    0xFFFF_FFFF_FFE0_0001, // 2^64 - 2^17 + 1 (≡ 1 mod 8192)
    0xFFFF_FFFF_FFC0_0001, // another NTT-friendly prime
    0xFFFF_FFFF_FF60_0001,
];

/// Primes ≡ 1 (mod 2·8192 = 16384)
pub const PRIMES_8192: &[u64] = &[
    0xFFFF_FFFF_FFC0_0001,
];

/// Create compact BFV parameters for fast tests and examples.
/// n=1024, single 40-bit NTT prime, p=257, gadget_base=2^16.
/// Using 40-bit q to ensure Delta = q/p >> relin noise with gadget_base=2^16.
/// Auxiliary 50-bit NTT prime for HPS RNS multiplication.
pub fn compact_bfv() -> Result<Arc<BfvParams>> {
    BfvParamsBuilder::new()
        .ring_degree(1024)
        .plain_modulus(257)
        // 40-bit NTT prime: 1099509805057 ≡ 1 mod 2048
        .ct_moduli(vec![1099509805057])
        // 50-bit auxiliary NTT prime for HPS: 562949953443841 ≡ 1 mod 2048
        // Need QP > 2·n·(q/2)² ≈ 2^89, QP ≈ 2^90. ✓
        .aux_moduli(vec![562949953443841])
        .sigma(3.2)
        .build()
}

/// Create small BFV parameters (≈80-bit security).
/// n=4096, p=65537.
pub fn small_bfv() -> Result<Arc<BfvParams>> {
    // q ≡ 1 (mod 2·4096 = 8192)
    // 4295049217 = 4096·1048576 + 1 ≡ 1 mod 8192
    // Actually, let's use: 4611686018326724609 ≡ 1 (mod 8192)
    // For now, use a known good prime:
    // 576460752308273153 is prime and ≡ 1 mod 8192
    BfvParamsBuilder::new()
        .ring_degree(4096)
        .plain_modulus(65537)
        .ct_moduli(vec![576460752308273153])
        .sigma(3.2)
        .build()
}

/// Create dBFV parameters for u64 arithmetic (~128-bit security).
///
/// p = 2^64 (represented as 0), base = 256 (2^8), d = 8 digits.
/// n = 4096, single 60-bit NTT prime Q, two auxiliary NTT primes for HPS.
///
/// BFV plaintext modulus: 1040407 (prime > 2*d*(b-1)^2 = 1040400).
/// Gadget base: 2^8 = 256. With Q≈2^60: Delta=Q/t≈1.1T, gadget_digits=8,
/// relin noise ~420k per mul, 8 muls summed ~1.2M << Delta/2 ≈ 554G.
pub fn u64_dbfv() -> Result<Arc<DbfvParams>> {
    let bfv = BfvParamsBuilder::new()
        .ring_degree(4096)
        .plain_modulus(1040407)
        // Single 60-bit NTT prime ≡ 1 mod 8192
        .ct_moduli(vec![1152921504606830593])
        // Two auxiliary NTT primes for HPS RNS multiplication.
        // P = product ≈ 2^109, need P > n·Q/2 ≈ 2^71. ✓
        .aux_moduli(vec![18014398509998081, 36028797018972161])
        .gadget_base(256) // 2^8: small enough for Delta/2 noise budget
        .sigma(3.2)
        .build()?;

    DbfvParams::new(bfv, 256, 8, 0) // p = 2^64 (0 sentinel = full u64 range)
}

/// Create compact dBFV parameters for fast tests and examples.
/// p=256, base=16, d=2 (so 16^2 = 256 = p).
///
/// BFV plaintext modulus: 929 (prime > 2·d·(b-1)^2 = 900).
/// This ensures convolution products don't overflow, and that centering
/// (interpreting values > t/2 as negative) is unambiguous for both
/// positive (≤ 450) and negative (≥ 914) values.
///
/// Using 40-bit q for sufficient noise budget with gadget_base=2^16.
pub fn compact_dbfv() -> Result<Arc<DbfvParams>> {
    let bfv = BfvParamsBuilder::new()
        .ring_degree(1024)
        .plain_modulus(929) // Prime > 2*d*(b-1)^2 = 900, for unambiguous centering
        // 40-bit NTT prime: 1099509805057 ≡ 1 mod 2048
        .ct_moduli(vec![1099509805057])
        // 50-bit auxiliary NTT prime for HPS multiplication
        .aux_moduli(vec![562949953443841])
        .sigma(3.2)
        .build()?;

    DbfvParams::new(bfv, 16, 2, 256)
}
