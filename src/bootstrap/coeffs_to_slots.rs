use std::collections::HashMap;

use crate::error::{ExactoError, Result};
use crate::ring::modular::mod_inv;
use crate::bfv::{BfvCiphertext, keygen::GaloisKey};
use crate::bfv::eval::{
    bfv_add, bfv_apply_automorphism, bfv_monomial_mul, bfv_plain_mul,
};
use crate::ring::poly::CoeffPoly;

/// Extract coefficient j from an encrypted polynomial.
///
/// Given ct encrypting m(X) = Σ a_i · X^i, returns a ciphertext encrypting
/// the scalar a_j (as a constant polynomial).
///
/// Algorithm:
/// 1. Multiply by X^{-j} = X^{2n-j} to shift coefficient j to position 0.
/// 2. Apply the homomorphic trace Tr = Σ_{k ∈ Gal} σ_k to isolate position 0.
///    Since Tr(X^0) = n and Tr(X^i) = 0 for i ≠ 0, we get n · a_j.
/// 3. Scale by n^{-1} mod t to recover a_j.
pub fn extract_coefficient(
    ct: &BfvCiphertext,
    j: usize,
    galois_keys: &HashMap<usize, GaloisKey>,
) -> Result<BfvCiphertext> {
    let params = &ct.params;
    let n = params.ring_degree;
    let t = params.plain_modulus;

    // Step 1: Multiply by X^{-j} = X^{2n-j} mod (X^n+1)
    let shifted = if j == 0 {
        ct.clone()
    } else {
        bfv_monomial_mul(ct, 2 * n - j)?
    };

    // Step 2: Compute trace Tr(shifted) = Σ_{k odd, 1≤k<2n} σ_k(shifted)
    let result = shifted_trace(&shifted, n, galois_keys)?;

    // Step 3: Scale by n^{-1} mod t
    let n_inv = mod_inv(n as u64 % t, t)
        .ok_or_else(|| ExactoError::InvalidParam("n not invertible mod t".into()))?;
    let mut scale_coeffs = vec![0u64; n];
    scale_coeffs[0] = n_inv;
    let scale_poly = CoeffPoly { coeffs: scale_coeffs, modulus: t };
    let result = bfv_plain_mul(&result, &scale_poly)?;

    Ok(result)
}

/// Compute the full trace: Tr(ct) = Σ_{k odd, 1≤k<2n} σ_k(ct).
///
/// Uses a doubling approach for efficiency when possible:
/// For n = 2^l, we can compute the trace in log2(n) steps.
fn shifted_trace(
    ct: &BfvCiphertext,
    n: usize,
    galois_keys: &HashMap<usize, GaloisKey>,
) -> Result<BfvCiphertext> {
    // Small rings and non-power-of-2 rings use the generic summation path.
    if n <= 32 || !n.is_power_of_two() {
        return naive_trace(ct, n, galois_keys);
    }

    // For n=2^ℓ, compute full trace with ℓ relative-trace steps:
    // Tr = Π_{s ∈ {n, n/2, ..., 2}} (1 + σ_{s+1}).
    // Each step doubles the subgroup size, giving O(log n) automorphisms.
    let mut result = ct.clone();
    for k in required_trace_elements(n) {
        let gk = galois_keys.get(&k).ok_or_else(|| {
            ExactoError::InvalidParam(format!("missing Galois key for element {k}"))
        })?;
        let auto = bfv_apply_automorphism(&result, gk)?;
        result = bfv_add(&result, &auto)?;
    }
    Ok(result)
}

/// Naive trace computation: directly sum σ_k(ct) for all odd k in [1, 2n).
fn naive_trace(
    ct: &BfvCiphertext,
    n: usize,
    galois_keys: &HashMap<usize, GaloisKey>,
) -> Result<BfvCiphertext> {
    let mut result = ct.clone(); // σ_1(ct) = ct (identity)

    for k in (3..2 * n).step_by(2) {
        let gk = galois_keys.get(&k).ok_or_else(|| {
            ExactoError::InvalidParam(format!("missing Galois key for element {k}"))
        })?;
        let auto_ct = bfv_apply_automorphism(ct, gk)?;
        result = bfv_add(&result, &auto_ct)?;
    }

    Ok(result)
}

/// CoeffsToSlots: extract all n coefficients from an encrypted polynomial.
///
/// Returns n ciphertexts, where the j-th ciphertext encrypts the j-th
/// coefficient of the original plaintext polynomial (as a constant polynomial).
///
/// This is the main building block for full ring bootstrapping (§5 of the paper).
pub fn coeffs_to_slots(
    ct: &BfvCiphertext,
    galois_keys: &HashMap<usize, GaloisKey>,
) -> Result<Vec<BfvCiphertext>> {
    let n = ct.params.ring_degree;
    let mut slots = Vec::with_capacity(n);

    for j in 0..n {
        let slot_j = extract_coefficient(ct, j, galois_keys)?;
        slots.push(slot_j);
    }

    Ok(slots)
}

/// SlotsToCoeffs: pack n scalar ciphertexts back into a polynomial ciphertext.
///
/// Given ct_0, ..., ct_{n-1} where ct_j encrypts a_j (as constant polynomial),
/// returns a ciphertext encrypting a_0 + a_1·X + ... + a_{n-1}·X^{n-1}.
pub fn slots_to_coeffs(
    slots: &[BfvCiphertext],
) -> Result<BfvCiphertext> {
    if slots.is_empty() {
        return Err(ExactoError::InvalidParam("empty slots".into()));
    }
    let n = slots[0].params.ring_degree;
    if slots.len() != n {
        return Err(ExactoError::InvalidParam(
            format!("expected {} slots, got {}", n, slots.len())
        ));
    }

    // result = Σ X^j · ct_j
    let mut result = slots[0].clone(); // X^0 · ct_0 = ct_0
    for j in 1..n {
        let shifted = bfv_monomial_mul(&slots[j], j)?;
        result = bfv_add(&result, &shifted)?;
    }

    Ok(result)
}

/// Generate all Galois keys needed for CoeffsToSlots.
///
/// Generates keys for all odd elements k in [3, 2n-1] (the non-identity
/// elements of the Galois group for Z[X]/(X^n+1)).
pub fn gen_all_galois_keys<R: rand::Rng>(
    sk: &crate::bfv::keygen::SecretKey,
    rng: &mut R,
) -> Result<HashMap<usize, GaloisKey>> {
    let n = sk.params.ring_degree;
    let mut keys = HashMap::new();

    for k in (3..2 * n).step_by(2) {
        let gk = crate::bfv::keygen::gen_galois_key_with_rng(sk, k, rng)?;
        keys.insert(k, gk);
    }

    Ok(keys)
}

/// Galois elements needed for coefficient-trace extraction.
///
/// For n=2^ℓ, this is the minimal relative-trace chain: {n+1, n/2+1, ..., 3}.
/// For non-power-of-2, falls back to all non-identity odd automorphisms.
pub fn required_trace_elements(n: usize) -> Vec<usize> {
    if n <= 32 {
        (3..2 * n).step_by(2).collect()
    } else if n.is_power_of_two() {
        let mut elems = Vec::new();
        let mut step = n;
        while step >= 2 {
            elems.push(step + 1);
            step >>= 1;
        }
        elems
    } else {
        (3..2 * n).step_by(2).collect()
    }
}

/// Generate only the Galois keys needed for trace-based CoeffsToSlots.
pub fn gen_trace_galois_keys<R: rand::Rng>(
    sk: &crate::bfv::keygen::SecretKey,
    rng: &mut R,
) -> Result<HashMap<usize, GaloisKey>> {
    let n = sk.params.ring_degree;
    let mut keys = HashMap::new();
    for k in required_trace_elements(n) {
        let gk = crate::bfv::keygen::gen_galois_key_with_rng(sk, k, rng)?;
        keys.insert(k, gk);
    }
    Ok(keys)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use crate::bfv::keygen::*;
    use crate::bfv::encrypt::*;
    use crate::params::{BfvParams, BfvParamsBuilder};
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    fn small_test_params() -> Result<Arc<BfvParams>> {
        // Small params for CoeffsToSlots testing: n=16, large q for noise budget
        BfvParamsBuilder::new()
            .ring_degree(16)
            .plain_modulus(97) // Prime, and gcd(16, 97) = 1 so n^{-1} exists
            .ct_moduli(vec![1125899906842817]) // 51-bit NTT prime ≡ 1 mod 32
            .sigma(3.2)
            .gadget_base(8)
            .build()
    }

    #[test]
    fn test_extract_constant_coefficient() {
        let params = small_test_params().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let n = params.ring_degree;

        let sk = gen_secret_key_with_rng(&params, &mut rng).unwrap();
        let galois_keys = gen_all_galois_keys(&sk, &mut rng).unwrap();

        // Encrypt polynomial m(X) = 42 (constant)
        let mut pt_coeffs = vec![0u64; n];
        pt_coeffs[0] = 42;
        let pt = CoeffPoly { coeffs: pt_coeffs, modulus: params.plain_modulus };
        let ct = encrypt_sk_with_rng(&pt, &sk, &params, &mut rng).unwrap();

        // Extract coefficient 0 — should be 42
        let ct_coeff0 = extract_coefficient(&ct, 0, &galois_keys).unwrap();
        let dec = decrypt(&ct_coeff0, &sk).unwrap();
        assert_eq!(dec.coeffs[0], 42, "expected coeff[0]=42, got {}", dec.coeffs[0]);
    }

    #[test]
    fn test_extract_nonzero_coefficient() {
        let params = small_test_params().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let n = params.ring_degree;

        let sk = gen_secret_key_with_rng(&params, &mut rng).unwrap();
        let galois_keys = gen_all_galois_keys(&sk, &mut rng).unwrap();

        // Encrypt polynomial m(X) = 5 + 10·X + 15·X^2
        let mut pt_coeffs = vec![0u64; n];
        pt_coeffs[0] = 5;
        pt_coeffs[1] = 10;
        pt_coeffs[2] = 15;
        let pt = CoeffPoly { coeffs: pt_coeffs, modulus: params.plain_modulus };
        let ct = encrypt_sk_with_rng(&pt, &sk, &params, &mut rng).unwrap();

        // Extract coefficient 0 — should be 5
        let ct_c0 = extract_coefficient(&ct, 0, &galois_keys).unwrap();
        let dec0 = decrypt(&ct_c0, &sk).unwrap();
        assert_eq!(dec0.coeffs[0], 5, "coeff[0]: expected 5, got {}", dec0.coeffs[0]);

        // Extract coefficient 1 — should be 10
        let ct_c1 = extract_coefficient(&ct, 1, &galois_keys).unwrap();
        let dec1 = decrypt(&ct_c1, &sk).unwrap();
        assert_eq!(dec1.coeffs[0], 10, "coeff[1]: expected 10, got {}", dec1.coeffs[0]);

        // Extract coefficient 2 — should be 15
        let ct_c2 = extract_coefficient(&ct, 2, &galois_keys).unwrap();
        let dec2 = decrypt(&ct_c2, &sk).unwrap();
        assert_eq!(dec2.coeffs[0], 15, "coeff[2]: expected 15, got {}", dec2.coeffs[0]);
    }

    #[test]
    fn test_coeffs_to_slots_roundtrip() {
        let params = small_test_params().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let n = params.ring_degree;
        let t = params.plain_modulus;

        let sk = gen_secret_key_with_rng(&params, &mut rng).unwrap();
        let galois_keys = gen_all_galois_keys(&sk, &mut rng).unwrap();

        // Encrypt m(X) = 1 + 2X + 3X^2 + ... + 16X^15
        let pt_coeffs: Vec<u64> = (1..=n as u64).map(|i| i % t).collect();
        let pt = CoeffPoly { coeffs: pt_coeffs.clone(), modulus: t };
        let ct = encrypt_sk_with_rng(&pt, &sk, &params, &mut rng).unwrap();

        // CoeffsToSlots
        let slots = coeffs_to_slots(&ct, &galois_keys).unwrap();
        assert_eq!(slots.len(), n);

        // Verify each slot
        for (j, slot_ct) in slots.iter().enumerate() {
            let dec = decrypt(slot_ct, &sk).unwrap();
            let expected = pt_coeffs[j];
            assert_eq!(
                dec.coeffs[0], expected,
                "slot[{j}]: expected {expected}, got {}", dec.coeffs[0]
            );
        }

        // SlotsToCoeffs
        let ct_recon = slots_to_coeffs(&slots).unwrap();
        let dec_recon = decrypt(&ct_recon, &sk).unwrap();

        for j in 0..n {
            assert_eq!(
                dec_recon.coeffs[j], pt_coeffs[j],
                "reconstructed coeff[{j}]: expected {}, got {}",
                pt_coeffs[j], dec_recon.coeffs[j]
            );
        }
    }

    #[test]
    fn test_required_trace_elements_power_of_two() {
        // Small n keeps naive trace keys.
        assert_eq!(required_trace_elements(8), vec![3, 5, 7, 9, 11, 13, 15]);
        // Large power-of-two n uses relative-trace chain.
        assert_eq!(required_trace_elements(64), vec![65, 33, 17, 9, 5, 3]);
    }
}
