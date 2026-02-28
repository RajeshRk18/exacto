use std::sync::Arc;
use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

use crate::error::{ExactoError, Result};
use crate::params::BfvParams;
use crate::ring::modular::{mod_inv, mod_mul};
use crate::ring::poly::CoeffPoly;
use crate::ring::rns::RnsPoly;
use crate::bfv::{BfvCiphertext, keygen::{SecretKey, PublicKey}};
use crate::sampling::{sample_uniform_poly, sample_gaussian_poly, sample_binary_poly};

/// Encrypt a plaintext polynomial using the public key.
///
/// ct = (pk0·u + e1 + Δ·m, pk1·u + e2)
/// where Δ = ⌊q/p⌋, u is binary, e1,e2 are Gaussian errors.
pub fn encrypt_pk(
    plaintext: &CoeffPoly,
    pk: &PublicKey,
    params: &Arc<BfvParams>,
) -> Result<BfvCiphertext> {
    let mut rng = ChaCha20Rng::from_os_rng();
    encrypt_pk_with_rng(plaintext, pk, params, &mut rng)
}

/// Encrypt with provided RNG (for deterministic testing).
pub fn encrypt_pk_with_rng<R: rand::Rng>(
    plaintext: &CoeffPoly,
    pk: &PublicKey,
    params: &Arc<BfvParams>,
    rng: &mut R,
) -> Result<BfvCiphertext> {
    let basis = &params.ct_basis;
    let n = params.ring_degree;

    // Scale plaintext: Δ·m mod each q_i
    let delta_m = scale_plaintext(plaintext, params)?;

    // Sample u ← binary
    let u_coeffs = sample_binary_poly(n, basis.moduli[0], rng);
    let u = RnsPoly::from_coeff_poly(&u_coeffs, basis)?;

    // Sample errors e1, e2
    let e1_coeffs = sample_gaussian_poly(n, basis.moduli[0], params.sigma, rng);
    let e1 = RnsPoly::from_coeff_poly(&e1_coeffs, basis)?;

    let e2_coeffs = sample_gaussian_poly(n, basis.moduli[0], params.sigma, rng);
    let e2 = RnsPoly::from_coeff_poly(&e2_coeffs, basis)?;

    // c0 = pk0·u + e1 + Δ·m
    let pk0_u = pk.pk0.mul(&u)?;
    let c0 = pk0_u.add(&e1)?.add(&delta_m)?;

    // c1 = pk1·u + e2
    let pk1_u = pk.pk1.mul(&u)?;
    let c1 = pk1_u.add(&e2)?;

    Ok(BfvCiphertext {
        c: vec![c0, c1],
        params: params.clone(),
    })
}

/// Encrypt using the secret key (symmetric encryption).
///
/// ct = (-a·s + e + Δ·m, a)
pub fn encrypt_sk(
    plaintext: &CoeffPoly,
    sk: &SecretKey,
    params: &Arc<BfvParams>,
) -> Result<BfvCiphertext> {
    let mut rng = ChaCha20Rng::from_os_rng();
    encrypt_sk_with_rng(plaintext, sk, params, &mut rng)
}

/// Encrypt with secret key and provided RNG.
pub fn encrypt_sk_with_rng<R: rand::Rng>(
    plaintext: &CoeffPoly,
    sk: &SecretKey,
    params: &Arc<BfvParams>,
    rng: &mut R,
) -> Result<BfvCiphertext> {
    let basis = &params.ct_basis;
    let n = params.ring_degree;

    let delta_m = scale_plaintext(plaintext, params)?;

    // Sample uniform a
    let a_coeffs = sample_uniform_poly(n, basis.moduli[0], rng);
    let a = RnsPoly::from_coeff_poly(&a_coeffs, basis)?;

    // Sample error e
    let e_coeffs = sample_gaussian_poly(n, basis.moduli[0], params.sigma, rng);
    let e = RnsPoly::from_coeff_poly(&e_coeffs, basis)?;

    // c0 = -a·s + e + Δ·m
    let a_s = a.mul(&sk.poly)?;
    let c0 = a_s.neg().add(&e)?.add(&delta_m)?;

    Ok(BfvCiphertext {
        c: vec![c0, a],
        params: params.clone(),
    })
}

/// Decrypt a BFV ciphertext.
///
/// For degree-1 ct (c0, c1): m = ⌊p · (c0 + c1·s) / q⌉ mod p
pub fn decrypt(ct: &BfvCiphertext, sk: &SecretKey) -> Result<CoeffPoly> {
    let params = &ct.params;
    let basis = &params.ct_basis;

    // Compute the "phase": c0 + c1·s + c2·s^2 + ...
    let mut phase = ct.c[0].clone();
    let mut s_power = sk.poly.clone();
    for i in 1..ct.c.len() {
        let c_i_s = ct.c[i].mul(&s_power)?;
        phase = phase.add(&c_i_s)?;
        if i < ct.c.len() - 1 {
            s_power = s_power.mul(&sk.poly)?;
        }
    }

    // Convert phase from RNS-NTT to coefficient form
    // Scale down: m_i = round(p * phase_i / Q) mod p
    // Works for both single-prime and multi-prime Q via CRT reconstruction.
    let p = params.plain_modulus;
    let p_big = BigUint::from(p);

    let coeff_components: Vec<CoeffPoly> = phase.components.iter()
        .map(|c| c.to_coeff_poly())
        .collect();

    let mut q_big = BigUint::one();
    for &qi in &basis.moduli {
        q_big *= BigUint::from(qi);
    }
    let half_q = &q_big >> 1;

    // Precompute CRT terms: (Q / q_i) * (Q / q_i)^{-1} mod q_i
    let mut crt_terms = Vec::with_capacity(basis.moduli.len());
    for &qi in &basis.moduli {
        let qi_big = BigUint::from(qi);
        let q_star = &q_big / &qi_big;
        let q_star_mod_qi = (&q_star % &qi_big)
            .to_u64()
            .ok_or_else(|| ExactoError::InvalidParam(
                "failed to reduce CRT factor modulo q_i".into()
            ))?;
        let inv = mod_inv(q_star_mod_qi, qi)
            .ok_or_else(|| ExactoError::InvalidParam("non-coprime ciphertext moduli".into()))?;
        crt_terms.push(q_star * BigUint::from(inv));
    }

    let mut result = vec![0u64; params.ring_degree];
    for coeff_idx in 0..params.ring_degree {
        let mut x = BigUint::zero();
        for (i, coeffs_i) in coeff_components.iter().enumerate() {
            x += &crt_terms[i] * BigUint::from(coeffs_i.coeffs[coeff_idx]);
        }
        x %= &q_big;

        let scaled = (&x * &p_big + &half_q) / &q_big;
        let reduced: BigUint = &scaled % &p_big;
        result[coeff_idx] = reduced
            .to_u64()
            .ok_or_else(|| ExactoError::InvalidParam(
                "decryption coefficient does not fit in u64".into()
            ))?;
    }

    Ok(CoeffPoly {
        coeffs: result,
        modulus: p,
    })
}

/// Scale plaintext by Δ = ⌊Q/p⌋ (Q = ∏q_i) and convert to RNS-NTT form.
pub(crate) fn scale_plaintext(plaintext: &CoeffPoly, params: &BfvParams) -> Result<RnsPoly> {
    let basis = &params.ct_basis;
    let delta_residues = compute_delta_residues(params)?;
    let mut components = Vec::with_capacity(basis.num_moduli());

    for (idx, (&qi, plan)) in basis.moduli.iter().zip(basis.plans.iter()).enumerate() {
        let delta_i = delta_residues[idx];
        let bk = basis.barrett_ks[idx];

        let coeffs: Vec<u64> = plaintext.coeffs.iter()
            .map(|&m| mod_mul(m % qi, delta_i, qi, bk))
            .collect();

        let cp = CoeffPoly { coeffs, modulus: qi };
        components.push(crate::ring::ntt::NttPoly::from_coeff_poly(&cp, plan.clone())?);
    }

    Ok(RnsPoly {
        components,
        ring_degree: params.ring_degree,
    })
}

/// Compute Δ residues for BFV encoding where Δ = floor(Q / p) and Q = ∏ q_i.
fn compute_delta_residues(params: &BfvParams) -> Result<Vec<u64>> {
    let basis = &params.ct_basis;
    let p_big = BigUint::from(params.plain_modulus);

    let mut q_big = BigUint::one();
    for &qi in &basis.moduli {
        q_big *= BigUint::from(qi);
    }

    let delta = &q_big / &p_big;
    if delta.is_zero() {
        return Err(ExactoError::InvalidParam(
            "ciphertext modulus product Q must be >= plaintext modulus p".into(),
        ));
    }

    basis.moduli.iter().map(|&qi| {
        let qi_big = BigUint::from(qi);
        (&delta % &qi_big)
            .to_u64()
            .ok_or_else(|| ExactoError::InvalidParam(
                "failed to reduce Δ modulo q_i".into()
            ))
    }).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bfv::keygen::*;
    use crate::bfv::encoding::*;
    use crate::params::BfvParamsBuilder;
    use crate::params::presets::compact_bfv;
    use num_bigint::BigUint;
    use num_traits::{One, ToPrimitive};
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_encrypt_decrypt_sk() {
        let params = compact_bfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let sk = gen_secret_key_with_rng(&params, &mut rng).unwrap();
        let pt = encode_scalar(42, &params).unwrap();
        let ct = encrypt_sk_with_rng(&pt, &sk, &params, &mut rng).unwrap();
        let decrypted = decrypt(&ct, &sk).unwrap();

        assert_eq!(decode_scalar(&decrypted), 42);
    }

    #[test]
    fn test_encrypt_decrypt_pk() {
        let params = compact_bfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let sk = gen_secret_key_with_rng(&params, &mut rng).unwrap();
        let pk = gen_public_key_with_rng(&sk, &mut rng).unwrap();

        let pt = encode_scalar(100, &params).unwrap();
        let ct = encrypt_pk_with_rng(&pt, &pk, &params, &mut rng).unwrap();
        let decrypted = decrypt(&ct, &sk).unwrap();

        assert_eq!(decode_scalar(&decrypted), 100);
    }

    #[test]
    fn test_encrypt_decrypt_zero() {
        let params = compact_bfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let sk = gen_secret_key_with_rng(&params, &mut rng).unwrap();
        let pt = encode_scalar(0, &params).unwrap();
        let ct = encrypt_sk_with_rng(&pt, &sk, &params, &mut rng).unwrap();
        let decrypted = decrypt(&ct, &sk).unwrap();

        assert_eq!(decode_scalar(&decrypted), 0);
    }

    #[test]
    fn test_encrypt_decrypt_poly() {
        let params = compact_bfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let sk = gen_secret_key_with_rng(&params, &mut rng).unwrap();
        let values = vec![1, 2, 3, 4, 5, 0, 0, 0];
        let pt = encode_simd(&values, &params).unwrap();
        let ct = encrypt_sk_with_rng(&pt, &sk, &params, &mut rng).unwrap();
        let decrypted = decrypt(&ct, &sk).unwrap();

        let decoded = decode_simd(&decrypted, 5);
        assert_eq!(decoded, vec![1, 2, 3, 4, 5]);
    }

    #[test]
    fn test_scale_plaintext_uses_full_q_over_p() {
        let params = BfvParamsBuilder::new()
            .ring_degree(16)
            .plain_modulus(257)
            .ct_moduli(vec![1099509805057, 562949953443841])
            .sigma(3.2)
            .build()
            .unwrap();

        let mut coeffs = vec![0u64; params.ring_degree];
        coeffs[0] = 42;
        coeffs[3] = 19;
        let pt = CoeffPoly { coeffs: coeffs.clone(), modulus: params.plain_modulus };

        let scaled = scale_plaintext(&pt, &params).unwrap();

        let mut q_big = BigUint::one();
        for &qi in &params.ct_basis.moduli {
            q_big *= BigUint::from(qi);
        }
        let delta = &q_big / BigUint::from(params.plain_modulus);

        for (i, &qi) in params.ct_basis.moduli.iter().enumerate() {
            let delta_i = (&delta % BigUint::from(qi)).to_u64().unwrap();
            let coeffs_i = scaled.components[i].to_coeff_poly();
            for (j, &m) in coeffs.iter().enumerate() {
                let expected = ((m as u128 * delta_i as u128) % qi as u128) as u64;
                assert_eq!(coeffs_i.coeffs[j], expected, "mismatch at component {i}, coeff {j}");
            }
        }
    }

    #[test]
    fn test_encrypt_decrypt_sk_multi_modulus() {
        let params = BfvParamsBuilder::new()
            .ring_degree(16)
            .plain_modulus(257)
            .ct_moduli(vec![1099509805057, 562949953443841])
            .sigma(3.2)
            .build()
            .unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(31415);

        let sk = gen_secret_key_with_rng(&params, &mut rng).unwrap();
        for value in [0u64, 1, 7, 42, 128, 256] {
            let pt = encode_scalar(value, &params).unwrap();
            let ct = encrypt_sk_with_rng(&pt, &sk, &params, &mut rng).unwrap();
            let decrypted = decrypt(&ct, &sk).unwrap();
            assert_eq!(decode_scalar(&decrypted), value);
        }
    }
}
