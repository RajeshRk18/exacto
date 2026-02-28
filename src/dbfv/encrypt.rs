use std::sync::Arc;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

use crate::error::{ExactoError, Result};
use crate::params::DbfvParams;
use crate::bfv::encrypt::{encrypt_pk_with_rng, encrypt_sk_with_rng};
use crate::bfv::keygen::{SecretKey, PublicKey};
use crate::dbfv::ciphertext::DbfvCiphertext;
use crate::dbfv::decomposition::{digit_decompose, poly_digit_decompose};
use crate::ring::poly::CoeffPoly;

/// Encrypt a plaintext μ ∈ Z_p using dBFV with the public key.
///
/// Method 1: decompose μ into digits, encrypt each digit as a BFV ciphertext.
/// Each BFV ciphertext encrypts μ_i ∈ Z_b with plaintext modulus b.
pub fn dbfv_encrypt(
    plaintext: u64,
    pk: &PublicKey,
    params: &Arc<DbfvParams>,
) -> Result<DbfvCiphertext> {
    let mut rng = ChaCha20Rng::from_os_rng();
    dbfv_encrypt_with_rng(plaintext, pk, params, &mut rng)
}

/// Encrypt with provided RNG.
pub fn dbfv_encrypt_with_rng<R: rand::Rng>(
    plaintext: u64,
    pk: &PublicKey,
    params: &Arc<DbfvParams>,
    rng: &mut R,
) -> Result<DbfvCiphertext> {
    let digits = digit_decompose_scalar(plaintext, params);
    encrypt_pk_digit_polys(&digits, pk, params, rng)
}

/// Encrypt a ring plaintext polynomial in Z_p[X]/(X^n+1) with the public key.
///
/// Each coefficient is decomposed into base-b digits; limb i encrypts the
/// polynomial of i-th digits across all coefficients.
pub fn dbfv_encrypt_poly(
    plaintext: &CoeffPoly,
    pk: &PublicKey,
    params: &Arc<DbfvParams>,
) -> Result<DbfvCiphertext> {
    let mut rng = ChaCha20Rng::from_os_rng();
    dbfv_encrypt_poly_with_rng(plaintext, pk, params, &mut rng)
}

/// Encrypt a ring plaintext polynomial with provided RNG.
pub fn dbfv_encrypt_poly_with_rng<R: rand::Rng>(
    plaintext: &CoeffPoly,
    pk: &PublicKey,
    params: &Arc<DbfvParams>,
    rng: &mut R,
) -> Result<DbfvCiphertext> {
    let reduced = normalize_poly_plaintext(plaintext, params)?;
    let digit_polys = poly_digit_decompose(&reduced, params.base, params.num_digits);
    encrypt_pk_digit_polys(&digit_polys, pk, params, rng)
}

/// Encrypt using the secret key (symmetric).
pub fn dbfv_encrypt_sk(
    plaintext: u64,
    sk: &SecretKey,
    params: &Arc<DbfvParams>,
) -> Result<DbfvCiphertext> {
    let mut rng = ChaCha20Rng::from_os_rng();
    dbfv_encrypt_sk_with_rng(plaintext, sk, params, &mut rng)
}

/// Encrypt with secret key and provided RNG.
pub fn dbfv_encrypt_sk_with_rng<R: rand::Rng>(
    plaintext: u64,
    sk: &SecretKey,
    params: &Arc<DbfvParams>,
    rng: &mut R,
) -> Result<DbfvCiphertext> {
    let digits = digit_decompose_scalar(plaintext, params);
    encrypt_sk_digit_polys(&digits, sk, params, rng)
}

/// Encrypt a ring plaintext polynomial in Z_p[X]/(X^n+1) using secret key.
pub fn dbfv_encrypt_poly_sk(
    plaintext: &CoeffPoly,
    sk: &SecretKey,
    params: &Arc<DbfvParams>,
) -> Result<DbfvCiphertext> {
    let mut rng = ChaCha20Rng::from_os_rng();
    dbfv_encrypt_poly_sk_with_rng(plaintext, sk, params, &mut rng)
}

/// Encrypt a ring plaintext polynomial with secret key and provided RNG.
pub fn dbfv_encrypt_poly_sk_with_rng<R: rand::Rng>(
    plaintext: &CoeffPoly,
    sk: &SecretKey,
    params: &Arc<DbfvParams>,
    rng: &mut R,
) -> Result<DbfvCiphertext> {
    let reduced = normalize_poly_plaintext(plaintext, params)?;
    let digit_polys = poly_digit_decompose(&reduced, params.base, params.num_digits);
    encrypt_sk_digit_polys(&digit_polys, sk, params, rng)
}

fn digit_decompose_scalar(plaintext: u64, params: &DbfvParams) -> Vec<Vec<u64>> {
    let reduced = if params.plain_modulus == 0 {
        plaintext
    } else {
        plaintext % params.plain_modulus
    };
    let digits = digit_decompose(reduced, params.base, params.num_digits);
    let n = params.bfv_params.ring_degree;
    digits.into_iter().map(|digit| {
        let mut coeffs = vec![0u64; n];
        coeffs[0] = digit;
        coeffs
    }).collect()
}

fn normalize_poly_plaintext(plaintext: &CoeffPoly, params: &DbfvParams) -> Result<Vec<u64>> {
    if params.plain_modulus == 0 {
        return Err(ExactoError::InvalidParam(
            "polynomial dBFV plaintext requires finite plain_modulus (plain_modulus=0 is scalar-only)"
                .into(),
        ));
    }
    if plaintext.len() != params.bfv_params.ring_degree {
        return Err(ExactoError::DimensionMismatch {
            expected: params.bfv_params.ring_degree,
            got: plaintext.len(),
        });
    }
    if plaintext.modulus != params.plain_modulus {
        return Err(ExactoError::ModulusMismatch);
    }
    Ok(plaintext
        .coeffs
        .iter()
        .map(|&c| c % params.plain_modulus)
        .collect())
}

fn validate_digit_polys(digit_polys: &[Vec<u64>], params: &DbfvParams) -> Result<()> {
    if digit_polys.len() != params.num_digits {
        return Err(ExactoError::DimensionMismatch {
            expected: params.num_digits,
            got: digit_polys.len(),
        });
    }
    let n = params.bfv_params.ring_degree;
    let t = params.bfv_params.plain_modulus;
    if params.base > t {
        return Err(ExactoError::InvalidParam(
            format!(
                "dBFV base {} must be <= BFV plaintext modulus {}",
                params.base, t
            ),
        ));
    }

    for digit_poly in digit_polys {
        if digit_poly.len() != n {
            return Err(ExactoError::DimensionMismatch {
                expected: n,
                got: digit_poly.len(),
            });
        }
        if digit_poly.iter().any(|&c| c >= t) {
            return Err(ExactoError::InvalidParam(
                "digit polynomial coefficient is outside BFV plaintext modulus".into(),
            ));
        }
    }
    Ok(())
}

fn encrypt_pk_digit_polys<R: rand::Rng>(
    digit_polys: &[Vec<u64>],
    pk: &PublicKey,
    params: &Arc<DbfvParams>,
    rng: &mut R,
) -> Result<DbfvCiphertext> {
    validate_digit_polys(digit_polys, params)?;
    let t = params.bfv_params.plain_modulus;

    let mut limbs = Vec::with_capacity(params.num_digits);
    for coeffs in digit_polys {
        let pt = CoeffPoly {
            coeffs: coeffs.clone(),
            modulus: t,
        };
        let ct = encrypt_pk_with_rng(&pt, pk, &params.bfv_params, rng)?;
        limbs.push(ct);
    }

    Ok(DbfvCiphertext {
        limbs,
        degree: params.num_digits,
        mul_depth: 0,
        params: params.clone(),
    })
}

fn encrypt_sk_digit_polys<R: rand::Rng>(
    digit_polys: &[Vec<u64>],
    sk: &SecretKey,
    params: &Arc<DbfvParams>,
    rng: &mut R,
) -> Result<DbfvCiphertext> {
    validate_digit_polys(digit_polys, params)?;
    let t = params.bfv_params.plain_modulus;

    let mut limbs = Vec::with_capacity(params.num_digits);
    for coeffs in digit_polys {
        let pt = CoeffPoly {
            coeffs: coeffs.clone(),
            modulus: t,
        };
        let ct = encrypt_sk_with_rng(&pt, sk, &params.bfv_params, rng)?;
        limbs.push(ct);
    }

    Ok(DbfvCiphertext {
        limbs,
        degree: params.num_digits,
        mul_depth: 0,
        params: params.clone(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bfv::keygen::*;
    use crate::ring::poly::CoeffPoly;
    use crate::params::presets::compact_dbfv;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_dbfv_encrypt() {
        let params = compact_dbfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let sk = gen_secret_key_with_rng(&params.bfv_params, &mut rng).unwrap();
        let pk = gen_public_key_with_rng(&sk, &mut rng).unwrap();

        let ct = dbfv_encrypt_with_rng(42, &pk, &params, &mut rng).unwrap();
        assert_eq!(ct.num_limbs(), 2); // d=2
        assert_eq!(ct.degree, 2);
    }

    #[test]
    fn test_dbfv_encrypt_poly() {
        let params = compact_dbfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(43);

        let sk = gen_secret_key_with_rng(&params.bfv_params, &mut rng).unwrap();
        let pk = gen_public_key_with_rng(&sk, &mut rng).unwrap();

        let mut coeffs = vec![0u64; params.bfv_params.ring_degree];
        coeffs[0] = 42;
        coeffs[1] = 100;
        coeffs[2] = 255;
        let pt = CoeffPoly {
            coeffs,
            modulus: params.plain_modulus,
        };

        let ct = dbfv_encrypt_poly_with_rng(&pt, &pk, &params, &mut rng).unwrap();
        assert_eq!(ct.num_limbs(), params.num_digits);
        assert_eq!(ct.degree, params.num_digits);
    }
}
