use rayon::prelude::*;

use crate::error::{ExactoError, Result};
use crate::bfv::{BfvCiphertext, eval as bfv_eval, keygen::RelinKey};
use crate::dbfv::ciphertext::DbfvCiphertext;
use crate::dbfv::reduction;

/// Homomorphic addition of two dBFV ciphertexts.
///
/// Limb-wise BFV addition: (ct1_i + ct2_i) for each digit i.
pub fn dbfv_add(
    ct1: &DbfvCiphertext,
    ct2: &DbfvCiphertext,
) -> Result<DbfvCiphertext> {
    if ct1.num_limbs() != ct2.num_limbs() {
        return Err(ExactoError::DimensionMismatch {
            expected: ct1.num_limbs(),
            got: ct2.num_limbs(),
        });
    }

    let limbs: Vec<BfvCiphertext> = ct1.limbs.iter()
        .zip(ct2.limbs.iter())
        .map(|(a, b)| bfv_eval::bfv_add(a, b))
        .collect::<Result<Vec<_>>>()?;

    Ok(DbfvCiphertext {
        limbs,
        degree: ct1.degree.max(ct2.degree),
        mul_depth: ct1.mul_depth.max(ct2.mul_depth),
        params: ct1.params.clone(),
    })
}

/// Homomorphic subtraction of two dBFV ciphertexts.
pub fn dbfv_sub(
    ct1: &DbfvCiphertext,
    ct2: &DbfvCiphertext,
) -> Result<DbfvCiphertext> {
    if ct1.num_limbs() != ct2.num_limbs() {
        return Err(ExactoError::DimensionMismatch {
            expected: ct1.num_limbs(),
            got: ct2.num_limbs(),
        });
    }

    let limbs: Vec<BfvCiphertext> = ct1.limbs.iter()
        .zip(ct2.limbs.iter())
        .map(|(a, b)| bfv_eval::bfv_sub(a, b))
        .collect::<Result<Vec<_>>>()?;

    Ok(DbfvCiphertext {
        limbs,
        degree: ct1.degree.max(ct2.degree),
        mul_depth: ct1.mul_depth.max(ct2.mul_depth),
        params: ct1.params.clone(),
    })
}

/// Negate a dBFV ciphertext.
pub fn dbfv_neg(ct: &DbfvCiphertext) -> DbfvCiphertext {
    let limbs = ct.limbs.iter()
        .map(|l| bfv_eval::bfv_neg(l))
        .collect();
    DbfvCiphertext {
        limbs,
        degree: ct.degree,
        mul_depth: ct.mul_depth,
        params: ct.params.clone(),
    }
}

/// Homomorphic multiplication of two dBFV ciphertexts.
///
/// This is the core dBFV operation:
/// 1. Compute polynomial convolution in B: result_k = Σ_{i+j=k} BfvMul(limb_i, limb_j)
///    producing 2d-1 limbs.
/// 2. Relinearize each BFV product.
/// 3. Apply degree reduction + lattice reduction.
///
/// The d² BFV multiplications are parallelized with rayon.
pub fn dbfv_mul(
    ct1: &DbfvCiphertext,
    ct2: &DbfvCiphertext,
    rlk: &RelinKey,
) -> Result<DbfvCiphertext> {
    let params = &ct1.params;
    let d = params.num_digits;

    if ct1.num_limbs() != d || ct2.num_limbs() != d {
        return Err(ExactoError::InvalidParam(
            "multiplication requires d-limb ciphertexts".into()
        ));
    }

    let next_depth = ct1.mul_depth.max(ct2.mul_depth) + 1;
    if next_depth > 1 {
        return Err(ExactoError::NotImplemented(
            "chained dBFV multiplication requires ciphertext-level lattice reduction (paper §4.6.2)"
                .into(),
        ));
    }

    // Step 1: Convolution in B — compute result[k] = Σ_{i+j=k} limb_i * limb_j
    // Total of d² BFV multiplications producing 2d-1 output limbs.
    let result_len = 2 * d - 1;

    // Collect all (i, j) pairs grouped by output index k
    let mut work_items: Vec<(usize, usize, usize)> = Vec::with_capacity(d * d);
    for i in 0..d {
        for j in 0..d {
            work_items.push((i, j, i + j));
        }
    }

    // Perform all d² BFV multiplications in parallel
    let products: Vec<(usize, BfvCiphertext)> = work_items.par_iter()
        .map(|&(i, j, k)| {
            let prod = bfv_eval::bfv_mul_and_relin(&ct1.limbs[i], &ct2.limbs[j], rlk)?;
            Ok((k, prod))
        })
        .collect::<Result<Vec<_>>>()?;

    // Sum products for each output limb k
    let mut result_limbs: Vec<Option<BfvCiphertext>> = vec![None; result_len];
    for (k, prod) in products {
        if let Some(ref existing) = result_limbs[k] {
            result_limbs[k] = Some(bfv_eval::bfv_add(existing, &prod)?);
        } else {
            result_limbs[k] = Some(prod);
        }
    }

    let limbs: Vec<BfvCiphertext> = result_limbs.into_iter()
        .map(|opt| opt.expect("missing limb"))
        .collect();

    let mut result = DbfvCiphertext {
        limbs,
        degree: result_len,
        mul_depth: next_depth,
        params: params.clone(),
    };

    // Step 2: Degree reduction + lattice reduction
    result = reduction::reduce(&result, rlk)?;

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bfv::keygen::*;
    use crate::dbfv::encrypt::*;
    use crate::dbfv::decrypt::*;
    use crate::params::presets::{compact_dbfv, u64_dbfv};
    use crate::ring::poly::CoeffPoly;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_dbfv_add() {
        let params = compact_dbfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let sk = gen_secret_key_with_rng(&params.bfv_params, &mut rng).unwrap();

        let ct1 = dbfv_encrypt_sk_with_rng(10, &sk, &params, &mut rng).unwrap();
        let ct2 = dbfv_encrypt_sk_with_rng(20, &sk, &params, &mut rng).unwrap();

        let ct_sum = dbfv_add(&ct1, &ct2).unwrap();
        let result = dbfv_decrypt(&ct_sum, &sk).unwrap();
        assert_eq!(result, 30);
    }

    #[test]
    fn test_dbfv_sub() {
        let params = compact_dbfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let sk = gen_secret_key_with_rng(&params.bfv_params, &mut rng).unwrap();

        // 50 - 20 = 30. This involves digit borrow:
        // 50 = [2, 3], 20 = [4, 1] → sub = [-2, 2] in signed representation.
        // With signed recomposition: -2 + 2*16 = 30. ✓
        let ct1 = dbfv_encrypt_sk_with_rng(50, &sk, &params, &mut rng).unwrap();
        let ct2 = dbfv_encrypt_sk_with_rng(20, &sk, &params, &mut rng).unwrap();

        let ct_diff = dbfv_sub(&ct1, &ct2).unwrap();
        let result = dbfv_decrypt(&ct_diff, &sk).unwrap();
        assert_eq!(result, 30);
    }

    #[test]
    fn test_dbfv_poly_add() {
        let params = compact_dbfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(52);

        let sk = gen_secret_key_with_rng(&params.bfv_params, &mut rng).unwrap();
        let mut a_coeffs = vec![0u64; params.bfv_params.ring_degree];
        let mut b_coeffs = vec![0u64; params.bfv_params.ring_degree];
        a_coeffs[..4].copy_from_slice(&[1, 2, 3, 4]);
        b_coeffs[..4].copy_from_slice(&[5, 6, 7, 8]);

        let pt_a = CoeffPoly {
            coeffs: a_coeffs,
            modulus: params.plain_modulus,
        };
        let pt_b = CoeffPoly {
            coeffs: b_coeffs,
            modulus: params.plain_modulus,
        };

        let ct_a = dbfv_encrypt_poly_sk_with_rng(&pt_a, &sk, &params, &mut rng).unwrap();
        let ct_b = dbfv_encrypt_poly_sk_with_rng(&pt_b, &sk, &params, &mut rng).unwrap();
        let ct_sum = dbfv_add(&ct_a, &ct_b).unwrap();
        let dec = dbfv_decrypt_poly(&ct_sum, &sk).unwrap();

        assert_eq!(&dec.coeffs[..4], &[6, 8, 10, 12]);
    }

    #[test]
    fn test_dbfv_mul() {
        let params = compact_dbfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let sk = gen_secret_key_with_rng(&params.bfv_params, &mut rng).unwrap();
        let rlk = gen_relin_key_with_rng(&sk, &mut rng).unwrap();

        // 3 * 7 = 21
        let ct1 = dbfv_encrypt_sk_with_rng(3, &sk, &params, &mut rng).unwrap();
        let ct2 = dbfv_encrypt_sk_with_rng(7, &sk, &params, &mut rng).unwrap();

        let ct_prod = dbfv_mul(&ct1, &ct2, &rlk).unwrap();
        let result = dbfv_decrypt(&ct_prod, &sk).unwrap();
        assert_eq!(result, 21);
    }

    #[test]
    fn test_dbfv_poly_mul() {
        let params = compact_dbfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(53);

        let sk = gen_secret_key_with_rng(&params.bfv_params, &mut rng).unwrap();
        let rlk = gen_relin_key_with_rng(&sk, &mut rng).unwrap();

        // (3 + X) * (2 + X) = 6 + 5X + X^2
        let mut a_coeffs = vec![0u64; params.bfv_params.ring_degree];
        let mut b_coeffs = vec![0u64; params.bfv_params.ring_degree];
        a_coeffs[0] = 3;
        a_coeffs[1] = 1;
        b_coeffs[0] = 2;
        b_coeffs[1] = 1;

        let pt_a = CoeffPoly {
            coeffs: a_coeffs,
            modulus: params.plain_modulus,
        };
        let pt_b = CoeffPoly {
            coeffs: b_coeffs,
            modulus: params.plain_modulus,
        };

        let ct_a = dbfv_encrypt_poly_sk_with_rng(&pt_a, &sk, &params, &mut rng).unwrap();
        let ct_b = dbfv_encrypt_poly_sk_with_rng(&pt_b, &sk, &params, &mut rng).unwrap();
        let ct_prod = dbfv_mul(&ct_a, &ct_b, &rlk).unwrap();
        let dec = dbfv_decrypt_poly(&ct_prod, &sk).unwrap();
        assert_eq!(&dec.coeffs[..3], &[6, 5, 1]);
    }

    #[test]
    fn test_dbfv_mul_larger() {
        let params = compact_dbfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let sk = gen_secret_key_with_rng(&params.bfv_params, &mut rng).unwrap();
        let rlk = gen_relin_key_with_rng(&sk, &mut rng).unwrap();

        // Test products that exceed base b=16 (verifying carry handling)
        for (a, b_val) in [(15, 15), (10, 20), (12, 12)] {
            let ct1 = dbfv_encrypt_sk_with_rng(a, &sk, &params, &mut rng).unwrap();
            let ct2 = dbfv_encrypt_sk_with_rng(b_val, &sk, &params, &mut rng).unwrap();

            let ct_prod = dbfv_mul(&ct1, &ct2, &rlk).unwrap();
            let result = dbfv_decrypt(&ct_prod, &sk).unwrap();
            let expected = (a as u64 * b_val as u64) % 256;
            assert_eq!(result, expected, "{a} * {b_val}: expected {expected}, got {result}");
        }
    }

    #[test]
    fn test_dbfv_mul_rejects_depth_above_one() {
        let params = compact_dbfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(777);

        let sk = gen_secret_key_with_rng(&params.bfv_params, &mut rng).unwrap();
        let rlk = gen_relin_key_with_rng(&sk, &mut rng).unwrap();

        let ct1 = dbfv_encrypt_sk_with_rng(3, &sk, &params, &mut rng).unwrap();
        let ct2 = dbfv_encrypt_sk_with_rng(7, &sk, &params, &mut rng).unwrap();
        let ct3 = dbfv_encrypt_sk_with_rng(5, &sk, &params, &mut rng).unwrap();

        let ct12 = dbfv_mul(&ct1, &ct2, &rlk).unwrap();
        assert_eq!(ct12.mul_depth, 1);

        let err = dbfv_mul(&ct12, &ct3, &rlk).unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains("chained dBFV multiplication requires ciphertext-level lattice reduction"),
            "unexpected depth-limit error: {msg}"
        );
    }

    #[test]
    fn test_u64_dbfv_encrypt_decrypt() {
        let params = u64_dbfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(99);

        let sk = gen_secret_key_with_rng(&params.bfv_params, &mut rng).unwrap();

        for val in [0u64, 1, 42, 255, 256, 65535, 1_000_000, u64::MAX - 1, u64::MAX] {
            let ct = dbfv_encrypt_sk_with_rng(val, &sk, &params, &mut rng).unwrap();
            let result = dbfv_decrypt(&ct, &sk).unwrap();
            assert_eq!(result, val, "encrypt/decrypt roundtrip failed for {val}");
        }
    }

    #[test]
    fn test_u64_dbfv_add() {
        let params = u64_dbfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(100);

        let sk = gen_secret_key_with_rng(&params.bfv_params, &mut rng).unwrap();

        let ct1 = dbfv_encrypt_sk_with_rng(1_000_000, &sk, &params, &mut rng).unwrap();
        let ct2 = dbfv_encrypt_sk_with_rng(2_000_000, &sk, &params, &mut rng).unwrap();

        let ct_sum = dbfv_add(&ct1, &ct2).unwrap();
        let result = dbfv_decrypt(&ct_sum, &sk).unwrap();
        assert_eq!(result, 3_000_000);
    }

    #[test]
    fn test_u64_bfv_mul_direct() {
        // Test BFV mul with u64_dbfv's BFV parameters directly
        use crate::bfv::encrypt::{encrypt_sk_with_rng, decrypt as bfv_decrypt_fn};
        use crate::bfv::encoding::{encode_scalar, decode_scalar};
        use crate::bfv::eval::{bfv_mul_and_relin, bfv_mul_no_relin};

        let params = u64_dbfv().unwrap();
        let bfv = &params.bfv_params;
        let mut rng = ChaCha20Rng::seed_from_u64(101);

        let sk = gen_secret_key_with_rng(bfv, &mut rng).unwrap();
        let rlk = gen_relin_key_with_rng(&sk, &mut rng).unwrap();

        // First test: mul without relin to isolate HPS scaling
        let pt_a = encode_scalar(3u64, bfv).unwrap();
        let pt_b = encode_scalar(7u64, bfv).unwrap();
        let ct_a = encrypt_sk_with_rng(&pt_a, &sk, bfv, &mut rng).unwrap();
        let ct_b = encrypt_sk_with_rng(&pt_b, &sk, bfv, &mut rng).unwrap();

        let ct_nr = bfv_mul_no_relin(&ct_a, &ct_b).unwrap();
        eprintln!("mul_no_relin: {} components", ct_nr.c.len());
        // Decrypt degree-2 ciphertext manually: phase = c0 + c1*s + c2*s^2
        // Just check relin doesn't change result
        let ct_prod = bfv_mul_and_relin(&ct_a, &ct_b, &rlk).unwrap();
        let result = decode_scalar(&bfv_decrypt_fn(&ct_prod, &sk).unwrap());
        eprintln!("3 * 7 = {result} (expected 21)");

        for (a, b, expected) in [(3u64, 7, 21), (0, 5, 0), (10, 20, 200), (100, 100, 10000)] {
            let pt_a = encode_scalar(a, bfv).unwrap();
            let pt_b = encode_scalar(b, bfv).unwrap();
            let ct_a = encrypt_sk_with_rng(&pt_a, &sk, bfv, &mut rng).unwrap();
            let ct_b = encrypt_sk_with_rng(&pt_b, &sk, bfv, &mut rng).unwrap();

            let ct_prod = bfv_mul_and_relin(&ct_a, &ct_b, &rlk).unwrap();
            let result = decode_scalar(&bfv_decrypt_fn(&ct_prod, &sk).unwrap());
            assert_eq!(result, expected, "BFV mul {a} * {b}: expected {expected}, got {result}");
        }
    }

    #[test]
    fn test_u64_bfv_mul_single_aux_rejects_unsupported_p() {
        use crate::bfv::encrypt::encrypt_sk_with_rng;
        use crate::bfv::encoding::encode_scalar;
        use crate::bfv::eval::bfv_mul_and_relin;
        use crate::params::BfvParamsBuilder;

        let bfv = BfvParamsBuilder::new()
            .ring_degree(4096)
            .plain_modulus(1040407)
            .ct_moduli(vec![18014398509506561])
            .aux_moduli(vec![36028797018972161])
            .gadget_base(256)
            .sigma(3.2)
            .build()
            .unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(101);

        let sk = gen_secret_key_with_rng(&bfv, &mut rng).unwrap();
        let rlk = gen_relin_key_with_rng(&sk, &mut rng).unwrap();

        for (a, b) in [(3u64, 7), (0, 5), (10, 20)] {
            let pt_a = encode_scalar(a, &bfv).unwrap();
            let pt_b = encode_scalar(b, &bfv).unwrap();
            let ct_a = encrypt_sk_with_rng(&pt_a, &sk, &bfv, &mut rng).unwrap();
            let ct_b = encrypt_sk_with_rng(&pt_b, &sk, &bfv, &mut rng).unwrap();

            let err = bfv_mul_and_relin(&ct_a, &ct_b, &rlk).unwrap_err();
            let msg = err.to_string();
            assert!(
                msg.contains("single aux prime too small"),
                "unexpected error for single-aux BFV mul {a} * {b}: {msg}"
            );
        }
    }

    #[test]
    fn test_u64_bfv_mul_schoolbook_rejects_overflow_params() {
        use crate::bfv::encrypt::encrypt_sk_with_rng;
        use crate::bfv::encoding::encode_scalar;
        use crate::bfv::eval::bfv_mul_and_relin;
        use crate::params::BfvParamsBuilder;

        let bfv = BfvParamsBuilder::new()
            .ring_degree(4096)
            .plain_modulus(1040407)
            .ct_moduli(vec![18014398509506561])
            .gadget_base(256)
            .sigma(3.2)
            .build()
            .unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(101);

        let sk = gen_secret_key_with_rng(&bfv, &mut rng).unwrap();
        let rlk = gen_relin_key_with_rng(&sk, &mut rng).unwrap();

        for (a, b) in [(3u64, 7), (0, 5), (10, 20)] {
            let pt_a = encode_scalar(a, &bfv).unwrap();
            let pt_b = encode_scalar(b, &bfv).unwrap();
            let ct_a = encrypt_sk_with_rng(&pt_a, &sk, &bfv, &mut rng).unwrap();
            let ct_b = encrypt_sk_with_rng(&pt_b, &sk, &bfv, &mut rng).unwrap();

            let err = bfv_mul_and_relin(&ct_a, &ct_b, &rlk).unwrap_err();
            let msg = err.to_string();
            assert!(
                msg.contains("schoolbook BFV multiplication can overflow i128"),
                "unexpected schoolbook error for BFV mul {a} * {b}: {msg}"
            );
        }
    }

    #[test]
    fn test_u64_bfv_mul_zero_noise() {
        // Test noise of BFV mul with zero encryptions (u64 params)
        use crate::bfv::encrypt::{encrypt_sk_with_rng, decrypt as bfv_decrypt_fn};
        use crate::bfv::encoding::{encode_scalar, decode_scalar};
        use crate::bfv::eval::bfv_mul_and_relin;

        let params = u64_dbfv().unwrap();
        let bfv = &params.bfv_params;
        let mut rng = ChaCha20Rng::seed_from_u64(101);

        let sk = gen_secret_key_with_rng(bfv, &mut rng).unwrap();
        let rlk = gen_relin_key_with_rng(&sk, &mut rng).unwrap();

        let q = bfv.ct_basis.moduli[0];
        let p = bfv.plain_modulus;
        let delta = q / p;
        let half_q = q / 2;
        eprintln!("gadget_base={}, gadget_digits={}", bfv.gadget_base, bfv.gadget_digits);
        eprintln!("q={q}, p={p}, Delta={delta}, Delta/2={}", delta / 2);

        // Measure phase noise for individual products and their sum
        let pt_zero = encode_scalar(0u64, bfv).unwrap();
        let mut sum_ct: Option<crate::bfv::BfvCiphertext> = None;
        for trial in 0..8 {
            let ct_a = encrypt_sk_with_rng(&pt_zero, &sk, bfv, &mut rng).unwrap();
            let ct_b = encrypt_sk_with_rng(&pt_zero, &sk, bfv, &mut rng).unwrap();
            let prod = bfv_mul_and_relin(&ct_a, &ct_b, &rlk).unwrap();

            // Measure phase noise at coefficient 0
            let mut phase = prod.c[0].clone();
            let c1s = prod.c[1].mul(&sk.poly).unwrap();
            phase = phase.add(&c1s).unwrap();
            let phase_coeffs = phase.to_coeff_poly(&bfv.ct_basis);
            let c0 = phase_coeffs.coeffs[0];
            let noise = if c0 > half_q { q - c0 } else { c0 };
            eprintln!("  trial {trial}: noise[0]={noise} (Delta/2={})", delta/2);

            sum_ct = Some(match sum_ct {
                Some(acc) => bfv_eval::bfv_add(&acc, &prod).unwrap(),
                None => prod,
            });
        }

        // Measure sum noise
        let sum = sum_ct.unwrap();
        let mut phase = sum.c[0].clone();
        let c1s = sum.c[1].mul(&sk.poly).unwrap();
        phase = phase.add(&c1s).unwrap();
        let phase_coeffs = phase.to_coeff_poly(&bfv.ct_basis);
        let c0 = phase_coeffs.coeffs[0];
        let noise = if c0 > half_q { q - c0 } else { c0 };
        let sum_digit = decode_scalar(&bfv_decrypt_fn(&sum, &sk).unwrap());
        eprintln!("  sum noise[0]={noise} (Delta/2={})", delta/2);
        eprintln!("  sum decrypts to: {sum_digit}");

        // Also check noise across ALL coefficients
        let mut max_noise = 0u64;
        for &c in &phase_coeffs.coeffs {
            let n = if c > half_q { q - c } else { c };
            if n > max_noise { max_noise = n; }
        }
        eprintln!("  max noise across all {} coeffs: {max_noise}", phase_coeffs.coeffs.len());
    }

    #[test]
    fn test_u64_dbfv_mul_small() {
        use crate::bfv::encrypt::decrypt as bfv_decrypt_fn;
        use crate::bfv::encoding::decode_scalar;

        let params = u64_dbfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(101);

        let sk = gen_secret_key_with_rng(&params.bfv_params, &mut rng).unwrap();
        let rlk = gen_relin_key_with_rng(&sk, &mut rng).unwrap();

        // 3 * 7 = 21 (single digit, simplest case)
        let ct1 = dbfv_encrypt_sk_with_rng(3, &sk, &params, &mut rng).unwrap();
        let ct2 = dbfv_encrypt_sk_with_rng(7, &sk, &params, &mut rng).unwrap();

        let ct_prod = dbfv_mul(&ct1, &ct2, &rlk).unwrap();

        // Debug: print individual limb values
        eprintln!("ct_prod has {} limbs, degree={}", ct_prod.limbs.len(), ct_prod.degree);
        for (i, limb) in ct_prod.limbs.iter().enumerate() {
            let pt = bfv_decrypt_fn(limb, &sk).unwrap();
            let digit = decode_scalar(&pt);
            eprintln!("  limb[{i}] = {digit} (t={})", params.bfv_params.plain_modulus);
        }

        let result = dbfv_decrypt(&ct_prod, &sk).unwrap();
        assert_eq!(result, 21, "3 * 7: expected 21, got {result}");
    }

    #[test]
    fn test_u64_dbfv_mul() {
        let params = u64_dbfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(101);

        let sk = gen_secret_key_with_rng(&params.bfv_params, &mut rng).unwrap();
        let rlk = gen_relin_key_with_rng(&sk, &mut rng).unwrap();

        // 1000 * 2000 = 2_000_000
        let ct1 = dbfv_encrypt_sk_with_rng(1000, &sk, &params, &mut rng).unwrap();
        let ct2 = dbfv_encrypt_sk_with_rng(2000, &sk, &params, &mut rng).unwrap();

        let ct_prod = dbfv_mul(&ct1, &ct2, &rlk).unwrap();
        let result = dbfv_decrypt(&ct_prod, &sk).unwrap();
        assert_eq!(result, 2_000_000, "1000 * 2000: expected 2000000, got {result}");
    }
}
