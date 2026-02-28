use crate::bfv::BfvCiphertext;
use crate::bfv::eval::{bfv_add, bfv_apply_automorphism, bfv_plain_mul, bfv_sub};
use crate::bfv::keygen::GaloisKey;
use crate::dbfv::ciphertext::DbfvCiphertext;
use crate::dbfv::decomposition::digit_decompose;
use crate::error::{ExactoError, Result};
use crate::params::DbfvParams;
use crate::ring::modular::mod_inv;
use crate::ring::poly::CoeffPoly;

/// Apply a BFV automorphism to all dBFV limbs.
///
/// This realizes dBFV automorphisms (§4.4.3) by applying τ coefficient-wise
/// to each BFV limb and key-switching back via the provided Galois key.
pub fn dbfv_apply_automorphism(
    ct: &DbfvCiphertext,
    gk: &GaloisKey,
) -> Result<DbfvCiphertext> {
    let limbs: Vec<BfvCiphertext> = ct.limbs.iter()
        .map(|limb| bfv_apply_automorphism(limb, gk))
        .collect::<Result<Vec<_>>>()?;

    Ok(DbfvCiphertext {
        limbs,
        degree: ct.degree,
        mul_depth: ct.mul_depth,
        params: ct.params.clone(),
    })
}

/// Divide by the dBFV base b using the φ_b map from §4.4.4.
///
/// For c(B)=c0 + B*c_tilde(B), computes:
///   φ_b(c) = c0/b + c_tilde(B)
/// This divides plaintext and modulus by b when divisibility conditions hold.
pub fn dbfv_div_by_base(
    ct: &DbfvCiphertext,
) -> Result<DbfvCiphertext> {
    let d = ct.params.num_digits;
    if d == 0 || ct.limbs.is_empty() {
        return Err(ExactoError::InvalidParam("empty dBFV ciphertext".into()));
    }

    let base = ct.params.base;
    let t = ct.params.bfv_params.plain_modulus;
    let base_inv_t = mod_inv(base % t, t).ok_or_else(|| {
        ExactoError::InvalidParam("base not invertible modulo BFV plaintext modulus".into())
    })?;

    let old_p = if ct.params.plain_modulus == 0 {
        1u128 << 64
    } else {
        ct.params.plain_modulus as u128
    };
    if old_p % base as u128 != 0 {
        return Err(ExactoError::InvalidParam(
            format!("plaintext modulus {} is not divisible by base {}", old_p, base)
        ));
    }
    let new_p = old_p / base as u128;
    let new_p_u64 = if new_p == (1u128 << 64) { 0 } else { new_p as u64 };

    let c0_div = bfv_scalar_plain_mul(&ct.limbs[0], base_inv_t)?;
    let zero = bfv_sub(&ct.limbs[d - 1], &ct.limbs[d - 1])?;

    let mut limbs = vec![zero.clone(); d];
    limbs[0] = if d >= 2 {
        bfv_add(&ct.limbs[1], &c0_div)?
    } else {
        c0_div
    };
    for i in 1..d {
        limbs[i] = if i + 1 < d {
            ct.limbs[i + 1].clone()
        } else {
            zero.clone()
        };
    }

    let new_params = DbfvParams::new(
        ct.params.bfv_params.clone(),
        ct.params.base,
        ct.params.num_digits,
        new_p_u64,
    )?;

    Ok(DbfvCiphertext {
        limbs,
        degree: ct.degree.saturating_sub(1).max(1),
        mul_depth: ct.mul_depth,
        params: new_params,
    })
}

/// Homomorphic change of base from b to b' (§4.5), for scalar Z_p plaintexts.
///
/// Implements the linear transform form (Lemma 4.8) by expressing each old
/// basis power b^i mod p as a digit vector in the new base b'.
pub fn dbfv_change_base(
    ct: &DbfvCiphertext,
    new_base: u64,
    new_num_digits: usize,
) -> Result<DbfvCiphertext> {
    if new_base < 2 {
        return Err(ExactoError::InvalidParam("new base must be >= 2".into()));
    }
    if new_num_digits == 0 {
        return Err(ExactoError::InvalidParam("new_num_digits must be >= 1".into()));
    }

    let old_base = ct.params.base as u128;
    let old_d = ct.params.num_digits;
    let p = if ct.params.plain_modulus == 0 {
        1u128 << 64
    } else {
        ct.params.plain_modulus as u128
    };

    // Build linear transform D where each column i is the base-b' digit
    // decomposition of (b^i mod p).
    let mut transform = vec![vec![0u64; old_d]; new_num_digits];
    let mut b_pow = 1u128;
    for i in 0..old_d {
        let repr = (b_pow % p) as u64;
        let digits = digit_decompose(repr, new_base, new_num_digits);
        for j in 0..new_num_digits {
            transform[j][i] = digits[j];
        }
        b_pow = (b_pow * old_base) % p;
    }

    let zero = bfv_sub(&ct.limbs[0], &ct.limbs[0])?;
    let mut new_limbs = vec![zero.clone(); new_num_digits];
    for j in 0..new_num_digits {
        let mut acc = zero.clone();
        for i in 0..old_d {
            let coeff = transform[j][i];
            if coeff == 0 {
                continue;
            }
            let term = bfv_scalar_plain_mul(&ct.limbs[i], coeff)?;
            acc = bfv_add(&acc, &term)?;
        }
        new_limbs[j] = acc;
    }

    let new_params = DbfvParams::new(
        ct.params.bfv_params.clone(),
        new_base,
        new_num_digits,
        ct.params.plain_modulus,
    )?;

    Ok(DbfvCiphertext {
        limbs: new_limbs,
        degree: new_num_digits,
        mul_depth: ct.mul_depth,
        params: new_params,
    })
}

fn bfv_scalar_plain_mul(ct: &BfvCiphertext, scalar: u64) -> Result<BfvCiphertext> {
    let params = &ct.params;
    let mut coeffs = vec![0u64; params.ring_degree];
    coeffs[0] = scalar % params.plain_modulus;
    let pt = CoeffPoly {
        coeffs,
        modulus: params.plain_modulus,
    };
    bfv_plain_mul(ct, &pt)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bfv::keygen::{gen_galois_key_with_rng, gen_secret_key_with_rng};
    use crate::dbfv::decrypt::dbfv_decrypt;
    use crate::dbfv::encrypt::dbfv_encrypt_sk_with_rng;
    use crate::params::presets::compact_dbfv;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_dbfv_apply_automorphism_scalar() {
        let params = compact_dbfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let sk = gen_secret_key_with_rng(&params.bfv_params, &mut rng).unwrap();
        let gk = gen_galois_key_with_rng(&sk, 3, &mut rng).unwrap();

        let ct = dbfv_encrypt_sk_with_rng(42, &sk, &params, &mut rng).unwrap();
        let ct_auto = dbfv_apply_automorphism(&ct, &gk).unwrap();
        let dec = dbfv_decrypt(&ct_auto, &sk).unwrap();
        assert_eq!(dec, 42);
    }

    #[test]
    fn test_dbfv_div_by_base() {
        let params = compact_dbfv().unwrap(); // base=16, p=256
        let mut rng = ChaCha20Rng::seed_from_u64(43);
        let sk = gen_secret_key_with_rng(&params.bfv_params, &mut rng).unwrap();

        let ct = dbfv_encrypt_sk_with_rng(48, &sk, &params, &mut rng).unwrap();
        let ct_div = dbfv_div_by_base(&ct).unwrap();
        let dec = dbfv_decrypt(&ct_div, &sk).unwrap();
        assert_eq!(dec, 3);
        assert_eq!(ct_div.params.plain_modulus, 16);
    }

    #[test]
    fn test_dbfv_change_base() {
        let params = compact_dbfv().unwrap(); // base=16, d=2, p=256
        let mut rng = ChaCha20Rng::seed_from_u64(44);
        let sk = gen_secret_key_with_rng(&params.bfv_params, &mut rng).unwrap();

        for value in [0u64, 1, 15, 42, 127, 255] {
            let ct = dbfv_encrypt_sk_with_rng(value, &sk, &params, &mut rng).unwrap();
            let ct_b4 = dbfv_change_base(&ct, 4, 4).unwrap();
            let dec = dbfv_decrypt(&ct_b4, &sk).unwrap();
            assert_eq!(dec, value);
        }
    }
}
