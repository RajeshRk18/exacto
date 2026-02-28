use crate::error::{ExactoError, Result};
use crate::bfv::keygen::SecretKey;
use crate::bfv::encrypt::decrypt as bfv_decrypt;
use crate::bfv::encoding::decode_scalar;
use crate::dbfv::ciphertext::DbfvCiphertext;
use crate::dbfv::decomposition::{digit_recompose_signed, poly_digit_recompose_signed};
use crate::ring::poly::CoeffPoly;

/// Decrypt a dBFV ciphertext.
///
/// 1. Decrypt each BFV limb to get digit μ_i ∈ Z_t (BFV plaintext modulus).
/// 2. Center each digit: values > t/2 are interpreted as negative.
/// 3. Recompose: μ = Σ centered(μ_i) · b^i mod p.
///
/// The signed recomposition is essential because homomorphic operations
/// (especially subtraction) produce "negative" digit values that wrap mod t.
/// After multiplication, digits may also exceed b but stay within [-t/2, t/2),
/// and the signed recomposition combined with mod-p reduction gives the
/// correct result.
pub fn dbfv_decrypt(ct: &DbfvCiphertext, sk: &SecretKey) -> Result<u64> {
    if ct.params.plain_modulus != 0 {
        let pt = dbfv_decrypt_poly(ct, sk)?;
        return Ok(decode_scalar(&pt));
    }

    let params = &ct.params;
    let bfv_plain_mod = params.bfv_params.plain_modulus;

    let mut digits = Vec::with_capacity(ct.num_limbs());
    for limb in &ct.limbs {
        let pt = bfv_decrypt(limb, sk)?;
        let digit = decode_scalar(&pt);
        digits.push(digit);
    }

    let num_digits_to_use = params.num_digits.min(digits.len());
    let result = digit_recompose_signed(
        &digits[..num_digits_to_use],
        params.base,
        params.plain_modulus,
        bfv_plain_mod,
    );

    Ok(result)
}

/// Decrypt a dBFV ciphertext to a ring plaintext in Z_p[X]/(X^n+1).
///
/// Each limb decrypts to a polynomial of base-b digits for all coefficients.
/// Signed recomposition is then applied coefficient-wise.
pub fn dbfv_decrypt_poly(ct: &DbfvCiphertext, sk: &SecretKey) -> Result<CoeffPoly> {
    let params = &ct.params;
    if params.plain_modulus == 0 {
        return Err(ExactoError::InvalidParam(
            "polynomial dBFV decrypt requires finite plain_modulus (plain_modulus=0 is scalar-only)"
                .into(),
        ));
    }

    let bfv_plain_mod = params.bfv_params.plain_modulus;
    let num_digits_to_use = params.num_digits.min(ct.num_limbs());
    let mut digit_polys = Vec::with_capacity(num_digits_to_use);
    for limb in ct.limbs.iter().take(num_digits_to_use) {
        let pt = bfv_decrypt(limb, sk)?;
        digit_polys.push(pt.coeffs);
    }

    let coeffs = poly_digit_recompose_signed(
        &digit_polys,
        params.base,
        params.plain_modulus,
        bfv_plain_mod,
    );

    Ok(CoeffPoly {
        coeffs,
        modulus: params.plain_modulus,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bfv::keygen::*;
    use crate::dbfv::encrypt::*;
    use crate::params::presets::compact_dbfv;
    use crate::ring::poly::CoeffPoly;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_dbfv_encrypt_decrypt() {
        let params = compact_dbfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let sk = gen_secret_key_with_rng(&params.bfv_params, &mut rng).unwrap();

        for value in [0u64, 1, 42, 100, 255] {
            let ct = dbfv_encrypt_sk_with_rng(value, &sk, &params, &mut rng).unwrap();
            let decrypted = dbfv_decrypt(&ct, &sk).unwrap();
            assert_eq!(decrypted, value, "failed for value={value}");
        }
    }

    #[test]
    fn test_dbfv_encrypt_decrypt_pk() {
        let params = compact_dbfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let sk = gen_secret_key_with_rng(&params.bfv_params, &mut rng).unwrap();
        let pk = gen_public_key_with_rng(&sk, &mut rng).unwrap();

        for value in [0u64, 1, 42, 200] {
            let ct = dbfv_encrypt_with_rng(value, &pk, &params, &mut rng).unwrap();
            let decrypted = dbfv_decrypt(&ct, &sk).unwrap();
            assert_eq!(decrypted, value, "failed for value={value}");
        }
    }

    #[test]
    fn test_dbfv_encrypt_decrypt_poly() {
        let params = compact_dbfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(45);

        let sk = gen_secret_key_with_rng(&params.bfv_params, &mut rng).unwrap();

        let mut coeffs = vec![0u64; params.bfv_params.ring_degree];
        coeffs[0] = 1;
        coeffs[1] = 42;
        coeffs[2] = 100;
        coeffs[3] = 255;
        let pt = CoeffPoly {
            coeffs: coeffs.clone(),
            modulus: params.plain_modulus,
        };

        let ct = dbfv_encrypt_poly_sk_with_rng(&pt, &sk, &params, &mut rng).unwrap();
        let dec = dbfv_decrypt_poly(&ct, &sk).unwrap();
        assert_eq!(dec.coeffs[..4], coeffs[..4]);
    }
}
