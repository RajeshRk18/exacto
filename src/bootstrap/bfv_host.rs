use std::collections::HashMap;
use std::sync::Arc;
use rayon::prelude::*;

use crate::error::{ExactoError, Result};
use crate::params::BfvParams;
use crate::ring::poly::CoeffPoly;
use crate::ring::rns::RnsPoly;
use crate::bfv::{BfvCiphertext, keygen::{SecretKey, RelinKey, GaloisKey}};
use crate::bfv::eval::{bfv_add, bfv_monomial_mul, bfv_plain_mul};
use crate::bootstrap::digit_extract::{
    trivial_encrypt_poly, eval_poly_homomorphic, compute_rounding_poly,
};
use crate::bootstrap::coeffs_to_slots::{coeffs_to_slots, gen_trace_galois_keys};
use crate::dbfv::ciphertext::DbfvCiphertext;
use crate::dbfv::eval::dbfv_mul;

/// Bootstrap key: encrypts the original secret key under the bootstrap scheme.
///
/// The bootstrap key is a BFV ciphertext (under boot_params) encrypting
/// the polynomial s(X) where s is the original secret key.
/// Uses circular security assumption (same s for both schemes).
pub struct BootstrapKey {
    /// BFV ciphertext encrypting s under boot_params.
    pub bsk: BfvCiphertext,
    /// Bootstrap scheme parameters.
    pub boot_params: Arc<BfvParams>,
    /// Relinearization key for the boot scheme.
    pub boot_rlk: RelinKey,
    /// Galois keys for CoeffsToSlots (full ring bootstrap).
    pub galois_keys: HashMap<usize, GaloisKey>,
    /// Precomputed rounding polynomial coefficients in Z_{t_boot}.
    pub rounding_poly: Vec<u64>,
    /// Original plaintext modulus.
    pub t_orig: u64,
    /// Modulus switch target (= t_boot).
    pub q_prime: u64,
}

/// Generate a bootstrap key.
///
/// Creates a BFV ciphertext that encrypts the secret key polynomial s(X)
/// under the bootstrap scheme parameters. Assumes circular security.
///
/// Parameters:
/// - `sk`: original secret key
/// - `boot_params`: BFV parameters for the bootstrap scheme (larger q, t_boot ≥ q')
/// - `q_prime`: modulus switch target (should equal boot_params.plain_modulus)
/// - `t_orig`: original plaintext modulus
pub fn gen_bootstrap_key<R: rand::Rng>(
    sk: &SecretKey,
    boot_params: &Arc<BfvParams>,
    q_prime: u64,
    t_orig: u64,
    rng: &mut R,
) -> Result<BootstrapKey> {
    let n = sk.params.ring_degree;
    if boot_params.ring_degree != n {
        return Err(ExactoError::InvalidParam(
            "boot params must have same ring degree".into()
        ));
    }

    // Get s in coefficient form from original secret key
    let s_coeffs = sk.poly.components[0].to_coeff_poly();

    // s is ternary: coefficients in {0, 1, q-1} (representing {0, 1, -1}).
    // Convert to Z_{t_boot}: 0→0, 1→1, q-1→t_boot-1.
    let q_orig = sk.params.ct_basis.moduli[0];
    let t_boot = boot_params.plain_modulus;
    let mut s_plaintext_coeffs = vec![0u64; n];
    for i in 0..n {
        let c = s_coeffs.coeffs[i];
        if c == 0 {
            s_plaintext_coeffs[i] = 0;
        } else if c == 1 {
            s_plaintext_coeffs[i] = 1;
        } else if c == q_orig - 1 {
            s_plaintext_coeffs[i] = t_boot - 1; // -1 mod t_boot
        } else {
            // General case: c mod t_boot (centered)
            let half_q = q_orig / 2;
            if c > half_q {
                let neg = q_orig - c;
                s_plaintext_coeffs[i] = t_boot - (neg % t_boot);
            } else {
                s_plaintext_coeffs[i] = c % t_boot;
            }
        }
    }

    let s_pt = CoeffPoly {
        coeffs: s_plaintext_coeffs,
        modulus: t_boot,
    };

    // Create secret key for boot scheme (same polynomial s, different NTT domain)
    let boot_sk = create_boot_sk(sk, boot_params)?;

    // Encrypt s as plaintext under boot scheme
    let bsk = crate::bfv::encrypt::encrypt_sk_with_rng(&s_pt, &boot_sk, boot_params, rng)?;

    // Generate relin key for boot scheme
    let boot_rlk = crate::bfv::keygen::gen_relin_key_with_rng(&boot_sk, rng)?;

    // Generate the minimal trace-key set for CoeffsToSlots.
    let galois_keys = gen_trace_galois_keys(&boot_sk, rng)?;

    // Precompute rounding polynomial
    let rounding_poly = compute_rounding_poly(t_orig, q_prime, t_boot);

    Ok(BootstrapKey {
        bsk,
        boot_params: boot_params.clone(),
        boot_rlk,
        galois_keys,
        rounding_poly,
        t_orig,
        q_prime,
    })
}

/// Bootstrap a single BFV ciphertext: refresh noise.
///
/// Algorithm:
/// 1. Modulus switch ct from q to q' (small).
/// 2. Re-encrypt phase: ct_phase = TrivialEnc(c0') + PlainMul(bsk, c1').
/// 3. For general ciphertexts: use CoeffsToSlots to extract each coefficient
///    of the phase polynomial, evaluate the rounding polynomial on each,
///    then pack back with SlotsToCoeffs.
///
/// Supports both trivial ciphertexts (fast path, no CoeffsToSlots needed)
/// and general ciphertexts (full ring bootstrap with CoeffsToSlots).
pub fn bfv_bootstrap(
    ct: &BfvCiphertext,
    bsk: &BootstrapKey,
) -> Result<BfvCiphertext> {
    let q = ct.params.ct_basis.moduli[0];
    let q_prime = bsk.q_prime;
    let boot_params = &bsk.boot_params;
    let n = ct.params.ring_degree;

    if ct.c.len() != 2 {
        return Err(ExactoError::InvalidParam(
            "bootstrap requires degree-1 ciphertext".into()
        ));
    }

    // Step 1: Modulus switch from q to q'
    let c0_coeffs = ct.c[0].to_coeff_poly(&ct.params.ct_basis);
    let c1_coeffs = ct.c[1].to_coeff_poly(&ct.params.ct_basis);

    let mut c0_prime = vec![0u64; n];
    let mut c1_prime = vec![0u64; n];
    for i in 0..n {
        c0_prime[i] = ((q_prime as u128 * c0_coeffs.coeffs[i] as u128 + q as u128 / 2) / q as u128) as u64;
        c1_prime[i] = ((q_prime as u128 * c1_coeffs.coeffs[i] as u128 + q as u128 / 2) / q as u128) as u64;
        c0_prime[i] %= q_prime;
        c1_prime[i] %= q_prime;
    }

    let t_boot = boot_params.plain_modulus;

    // Step 2: Re-encrypt the phase under boot scheme
    let c0_pt = CoeffPoly {
        coeffs: c0_prime.iter().map(|&c| c % t_boot).collect(),
        modulus: t_boot,
    };
    let c1_pt = CoeffPoly {
        coeffs: c1_prime.iter().map(|&c| c % t_boot).collect(),
        modulus: t_boot,
    };

    let ct_c0 = trivial_encrypt_poly(&c0_pt, boot_params)?;
    let ct_c1s = bfv_plain_mul(&bsk.bsk, &c1_pt)?;
    let ct_phase = bfv_add(&ct_c0, &ct_c1s)?;

    // Step 3: Check if we can skip CoeffsToSlots (trivial ciphertext fast path)
    let is_trivial = c1_coeffs.coeffs.iter().all(|&c| c == 0);

    if is_trivial {
        // Fast path: phase is just c0 (constant-ish), evaluate directly
        let ct_result = eval_poly_homomorphic(&ct_phase, &bsk.rounding_poly, &bsk.boot_rlk)?;
        return Ok(ct_result);
    }

    // Step 3b: Full ring bootstrap with CoeffsToSlots
    // 3b.1: Extract each coefficient of the encrypted phase polynomial
    let coeff_cts = coeffs_to_slots(&ct_phase, &bsk.galois_keys)?;

    // 3b.2: Evaluate rounding polynomial on each coefficient independently
    let rounded_cts: Vec<BfvCiphertext> = coeff_cts.iter()
        .map(|ct_j| eval_poly_homomorphic(ct_j, &bsk.rounding_poly, &bsk.boot_rlk))
        .collect::<Result<Vec<_>>>()?;

    // 3b.3: Pack back with SlotsToCoeffs: result = Σ X^j · rounded_j
    let mut result = rounded_cts[0].clone(); // j=0: no monomial multiply needed
    for j in 1..n {
        let shifted = bfv_monomial_mul(&rounded_cts[j], j)?;
        result = bfv_add(&result, &shifted)?;
    }

    Ok(result)
}

/// dBFV bootstrapping: refresh noise in all d BFV limbs independently.
///
/// Algorithm (paper §5):
/// 1. Each dBFV limb is already a BFV ciphertext (one per digit).
/// 2. Bootstrap each BFV ciphertext independently (parallelizable with rayon).
/// 3. Return the refreshed dBFV ciphertext.
pub fn dbfv_bootstrap(
    ct: &DbfvCiphertext,
    bsk: &BootstrapKey,
) -> Result<DbfvCiphertext> {
    // Bootstrapping changes the underlying BFV parameter set. Keep dBFV metadata
    // (base, digits, plaintext modulus) and swap in the new BFV params.
    let refreshed_params = crate::params::DbfvParams::new(
        bsk.boot_params.clone(),
        ct.params.base,
        ct.params.num_digits,
        ct.params.plain_modulus,
    )?;

    let limbs: Vec<BfvCiphertext> = ct.limbs.par_iter()
        .map(|limb| bfv_bootstrap(limb, bsk))
        .collect::<Result<Vec<_>>>()?;

    Ok(DbfvCiphertext {
        limbs,
        degree: ct.degree,
        // Noise is refreshed; restart depth accounting.
        mul_depth: 0,
        params: refreshed_params,
    })
}

/// Multiply two dBFV ciphertexts and immediately refresh with bootstrapping.
///
/// This is the recommended way to extend multiplicative chains beyond the
/// single-depth safety bound in `dbfv_mul`.
pub fn dbfv_mul_then_bootstrap(
    ct1: &DbfvCiphertext,
    ct2: &DbfvCiphertext,
    rlk: &RelinKey,
    bsk: &BootstrapKey,
) -> Result<DbfvCiphertext> {
    let ct_mul = dbfv_mul(ct1, ct2, rlk)?;
    dbfv_bootstrap(&ct_mul, bsk)
}

/// Create a SecretKey for the boot scheme from the original secret key.
/// Same polynomial s, but in the boot NTT domain.
fn create_boot_sk(
    sk: &SecretKey,
    boot_params: &Arc<BfvParams>,
) -> Result<SecretKey> {
    let s_coeffs = sk.poly.components[0].to_coeff_poly();
    let q_orig = sk.params.ct_basis.moduli[0];
    let q_boot = boot_params.ct_basis.moduli[0];

    // Convert ternary coefficients from mod q_orig to mod q_boot
    let mut boot_coeffs = vec![0u64; s_coeffs.len()];
    for i in 0..s_coeffs.len() {
        let c = s_coeffs.coeffs[i];
        if c == 0 {
            boot_coeffs[i] = 0;
        } else if c <= q_orig / 2 {
            boot_coeffs[i] = c % q_boot;
        } else {
            // Negative value: c represents c - q_orig
            let neg = q_orig - c;
            boot_coeffs[i] = q_boot - (neg % q_boot);
        }
    }

    let boot_poly = CoeffPoly {
        coeffs: boot_coeffs,
        modulus: q_boot,
    };
    let poly = RnsPoly::from_coeff_poly(&boot_poly, &boot_params.ct_basis)?;

    Ok(SecretKey {
        poly,
        params: boot_params.clone(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bfv::keygen::*;
    use crate::bfv::encrypt::*;
    use crate::bfv::encoding::*;
    use crate::dbfv::encrypt::dbfv_encrypt_poly_sk_with_rng;
    use crate::dbfv::decrypt::{dbfv_decrypt, dbfv_decrypt_poly};
    use crate::dbfv::eval::dbfv_mul;
    use crate::params::{BfvParamsBuilder, DbfvParams};
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    /// Bootstrap test parameters.
    ///
    /// Original: n=16, q=65537, t=5.
    /// Boot: n=16, Q_boot=1125899906842817, t_boot=29, gadget_base=8.
    /// Switch modulus: q' = 25 (< t_boot).
    fn bootstrap_test_params() -> Result<(Arc<BfvParams>, Arc<BfvParams>, u64)> {
        let orig_params = BfvParamsBuilder::new()
            .ring_degree(16)
            .plain_modulus(5)
            .ct_moduli(vec![65537]) // Fermat prime, ≡ 1 mod 32
            .sigma(3.2)
            .build()?;

        let boot_params = BfvParamsBuilder::new()
            .ring_degree(16)
            .plain_modulus(29) // Prime > q'=25
            .ct_moduli(vec![1125899906842817]) // 51-bit NTT prime ≡ 1 mod 32
            .sigma(3.2)
            .gadget_base(8)
            .build()?;

        let q_prime = 25u64;

        Ok((orig_params, boot_params, q_prime))
    }

    fn dbfv_bootstrap_test_params() -> Result<(Arc<DbfvParams>, Arc<BfvParams>, u64)> {
        // Small dbfv instance for end-to-end API testing.
        let orig_bfv = BfvParamsBuilder::new()
            .ring_degree(16)
            .plain_modulus(97) // > 2*d*(b-1)^2 with b=4,d=2 => 36
            .ct_moduli(vec![65537])
            .sigma(3.2)
            .gadget_base(8)
            .build()?;
        let dbfv = DbfvParams::new(orig_bfv, 4, 2, 16)?;

        let boot_bfv = BfvParamsBuilder::new()
            .ring_degree(16)
            .plain_modulus(257) // prime > q'
            .ct_moduli(vec![1125899906842817]) // 51-bit NTT prime ≡ 1 mod 32
            .sigma(3.2)
            .gadget_base(8)
            .build()?;

        Ok((dbfv, boot_bfv, 64))
    }

    /// Test bootstrap pipeline using trivial ciphertexts (c1=0).
    ///
    /// With c1=0 the phase is just c0 (a constant polynomial after modswitching),
    /// so the rounding polynomial evaluates on a scalar value in [0, q').
    /// Full ring bootstrapping of general ciphertexts requires CoeffsToSlots
    /// (Chen-Han 2018) which uses Galois automorphisms to extract individual
    /// coefficients into SIMD slots before polynomial evaluation.
    #[test]
    fn test_bootstrap_single() {
        let (orig_params, boot_params, q_prime) = bootstrap_test_params().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let t_orig = orig_params.plain_modulus; // 5

        // Generate original secret key (needed for boot key generation)
        let sk = gen_secret_key_with_rng(&orig_params, &mut rng).unwrap();

        // Generate bootstrap key
        let bsk = gen_bootstrap_key(&sk, &boot_params, q_prime, t_orig, &mut rng).unwrap();
        let boot_sk = create_boot_sk(&sk, &boot_params).unwrap();

        // Test with trivial encryptions (c1=0): phase = c0 = Δ*m (constant poly).
        // After modswitching, phase is a single scalar in [0, q'), so the
        // rounding polynomial correctly extracts m.
        for m in [0u64, 1, 2, 3, 4] {
            let ct = crate::bootstrap::digit_extract::trivial_encrypt(m, &orig_params).unwrap();

            // Bootstrap
            let ct_boot = bfv_bootstrap(&ct, &bsk).unwrap();

            // Decrypt under boot scheme
            let dec_boot = decrypt(&ct_boot, &boot_sk).unwrap();
            let result = decode_scalar(&dec_boot) % t_orig;

            assert_eq!(result, m, "bootstrap failed for m={m}: got {result}");
        }
    }

    /// Test full ring bootstrapping with a real (non-trivial) ciphertext.
    ///
    /// This uses CoeffsToSlots to extract individual polynomial coefficients,
    /// evaluates the rounding polynomial on each, and packs them back.
    #[test]
    fn test_bootstrap_ring() {
        let (orig_params, boot_params, q_prime) = bootstrap_test_params().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let t_orig = orig_params.plain_modulus; // 5

        let sk = gen_secret_key_with_rng(&orig_params, &mut rng).unwrap();
        let bsk = gen_bootstrap_key(&sk, &boot_params, q_prime, t_orig, &mut rng).unwrap();
        let boot_sk = create_boot_sk(&sk, &boot_params).unwrap();

        // Encrypt scalar m=3 with secret key encryption (c1 ≠ 0, general ciphertext)
        let pt = encode_scalar(3, &orig_params).unwrap();
        let ct = encrypt_sk_with_rng(&pt, &sk, &orig_params, &mut rng).unwrap();

        // Verify it decrypts correctly before bootstrap
        let pre_dec = decrypt(&ct, &sk).unwrap();
        assert_eq!(decode_scalar(&pre_dec), 3, "pre-bootstrap decrypt failed");

        // Bootstrap (full ring path — c1 ≠ 0)
        let ct_boot = bfv_bootstrap(&ct, &bsk).unwrap();

        // Decrypt under boot scheme
        let dec_boot = decrypt(&ct_boot, &boot_sk).unwrap();
        let result = decode_scalar(&dec_boot) % t_orig;
        assert_eq!(result, 3, "ring bootstrap failed: expected 3, got {result}");
    }

    #[test]
    fn test_dbfv_mul_then_bootstrap_allows_next_mul() {
        let (dbfv_params, boot_params, q_prime) = dbfv_bootstrap_test_params().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(777);

        let sk = gen_secret_key_with_rng(&dbfv_params.bfv_params, &mut rng).unwrap();
        let rlk = gen_relin_key_with_rng(&sk, &mut rng).unwrap();

        let bsk = gen_bootstrap_key(
            &sk,
            &boot_params,
            q_prime,
            dbfv_params.bfv_params.plain_modulus,
            &mut rng,
        ).unwrap();
        let boot_sk = create_boot_sk(&sk, &boot_params).unwrap();
        let boot_dbfv_params = DbfvParams::new(
            boot_params.clone(),
            dbfv_params.base,
            dbfv_params.num_digits,
            dbfv_params.plain_modulus,
        ).unwrap();

        // First multiply under original params, then refresh.
        let mut c1 = vec![0u64; dbfv_params.bfv_params.ring_degree];
        c1[0] = 3;
        c1[1] = 1;
        let pt1 = CoeffPoly {
            coeffs: c1,
            modulus: dbfv_params.plain_modulus,
        };
        let mut c2 = vec![0u64; dbfv_params.bfv_params.ring_degree];
        c2[0] = 2;
        let pt2 = CoeffPoly {
            coeffs: c2,
            modulus: dbfv_params.plain_modulus,
        };
        let ct1 = dbfv_encrypt_poly_sk_with_rng(&pt1, &sk, &dbfv_params, &mut rng).unwrap();
        let ct2 = dbfv_encrypt_poly_sk_with_rng(&pt2, &sk, &dbfv_params, &mut rng).unwrap();
        let ct_refreshed = dbfv_mul_then_bootstrap(&ct1, &ct2, &rlk, &bsk).unwrap();
        assert_eq!(ct_refreshed.mul_depth, 0);
        assert_eq!(ct_refreshed.params.bfv_params.plain_modulus, boot_params.plain_modulus);

        // Multiply again under boot params (using boot key material).
        let mut c3 = vec![0u64; boot_dbfv_params.bfv_params.ring_degree];
        c3[0] = 2;
        let pt3 = CoeffPoly {
            coeffs: c3,
            modulus: boot_dbfv_params.plain_modulus,
        };
        let ct3 = dbfv_encrypt_poly_sk_with_rng(&pt3, &boot_sk, &boot_dbfv_params, &mut rng).unwrap();
        let ct_next = dbfv_mul(&ct_refreshed, &ct3, &bsk.boot_rlk).unwrap();
        assert_eq!(ct_next.mul_depth, 1);

        // Smoke check that refreshed ciphertext remains decryptable under boot key
        // through both poly and scalar compatibility APIs.
        let dec_poly = dbfv_decrypt_poly(&ct_next, &boot_sk).unwrap();
        assert_eq!(dec_poly.modulus, boot_dbfv_params.plain_modulus);
        assert_eq!(dec_poly.coeffs.len(), boot_dbfv_params.bfv_params.ring_degree);
        let _ = dbfv_decrypt(&ct_next, &boot_sk).unwrap();
    }
}
