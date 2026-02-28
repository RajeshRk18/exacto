use proptest::prelude::*;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

use exacto::bfv::encoding::{decode_scalar, encode_scalar};
use exacto::bfv::encrypt::{decrypt, encrypt_sk_with_rng};
use exacto::bfv::keygen::{gen_relin_key_with_rng, gen_secret_key_with_rng};
use exacto::bfv::eval::bfv_mul_and_relin;
use exacto::dbfv::decrypt::{dbfv_decrypt, dbfv_decrypt_poly};
use exacto::dbfv::encrypt::{dbfv_encrypt_poly_sk_with_rng, dbfv_encrypt_sk_with_rng};
use exacto::dbfv::eval::{dbfv_add, dbfv_mul};
use exacto::params::presets::{compact_bfv, compact_dbfv};
use exacto::ring::poly::CoeffPoly;

fn terms_to_poly(terms: &[(usize, u64)], n: usize, modulus: u64) -> CoeffPoly {
    let mut coeffs = vec![0u64; n];
    for &(idx, value) in terms {
        let i = idx % n;
        coeffs[i] = (coeffs[i] + value % modulus) % modulus;
    }
    CoeffPoly { coeffs, modulus }
}

fn sparse_negacyclic_mul(
    a_terms: &[(usize, u64)],
    b_terms: &[(usize, u64)],
    n: usize,
    modulus: u64,
) -> Vec<u64> {
    let mut result = vec![0u64; n];
    for &(i, a) in a_terms {
        let a_mod = a % modulus;
        if a_mod == 0 {
            continue;
        }
        for &(j, b) in b_terms {
            let b_mod = b % modulus;
            if b_mod == 0 {
                continue;
            }
            let prod = ((a_mod as u128 * b_mod as u128) % modulus as u128) as u64;
            let idx = i + j;
            if idx < n {
                result[idx] = (result[idx] + prod) % modulus;
            } else {
                let wrap_idx = idx - n;
                result[wrap_idx] = (result[wrap_idx] + modulus - prod) % modulus;
            }
        }
    }
    result
}

proptest! {
    #[test]
    fn prop_bfv_roundtrip_scalar(m in 0u64..257u64, seed in any::<u64>()) {
        let params = compact_bfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(seed);
        let sk = gen_secret_key_with_rng(&params, &mut rng).unwrap();

        let pt = encode_scalar(m, &params).unwrap();
        let ct = encrypt_sk_with_rng(&pt, &sk, &params, &mut rng).unwrap();
        let dec = decrypt(&ct, &sk).unwrap();
        prop_assert_eq!(decode_scalar(&dec), m);
    }

    #[test]
    fn prop_bfv_mul_scalar(a in 0u64..64u64, b in 0u64..64u64, seed in any::<u64>()) {
        let params = compact_bfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(seed);
        let sk = gen_secret_key_with_rng(&params, &mut rng).unwrap();
        let rlk = gen_relin_key_with_rng(&sk, &mut rng).unwrap();

        let pt_a = encode_scalar(a, &params).unwrap();
        let pt_b = encode_scalar(b, &params).unwrap();
        let ct_a = encrypt_sk_with_rng(&pt_a, &sk, &params, &mut rng).unwrap();
        let ct_b = encrypt_sk_with_rng(&pt_b, &sk, &params, &mut rng).unwrap();
        let ct_prod = bfv_mul_and_relin(&ct_a, &ct_b, &rlk).unwrap();
        let dec = decrypt(&ct_prod, &sk).unwrap();
        prop_assert_eq!(decode_scalar(&dec), (a * b) % params.plain_modulus);
    }

    #[test]
    fn prop_dbfv_add_mul(a in 0u64..256u64, b in 0u64..256u64, seed in any::<u64>()) {
        let params = compact_dbfv().unwrap();
        let p = if params.plain_modulus == 0 { u64::MAX } else { params.plain_modulus };
        let mut rng = ChaCha20Rng::seed_from_u64(seed);
        let sk = gen_secret_key_with_rng(&params.bfv_params, &mut rng).unwrap();
        let rlk = gen_relin_key_with_rng(&sk, &mut rng).unwrap();

        let ct_a = dbfv_encrypt_sk_with_rng(a, &sk, &params, &mut rng).unwrap();
        let ct_b = dbfv_encrypt_sk_with_rng(b, &sk, &params, &mut rng).unwrap();

        let ct_sum = dbfv_add(&ct_a, &ct_b).unwrap();
        let dec_sum = dbfv_decrypt(&ct_sum, &sk).unwrap();
        prop_assert_eq!(dec_sum, (a + b) % p);

        // Single multiplication layer (matching compact_dbfv safety envelope).
        let ct_prod = dbfv_mul(&ct_a, &ct_b, &rlk).unwrap();
        let dec_prod = dbfv_decrypt(&ct_prod, &sk).unwrap();
        prop_assert_eq!(dec_prod, (a * b) % p);
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn prop_dbfv_poly_add_sparse(
        a_terms in proptest::collection::vec((0usize..1024usize, 0u64..16u64), 0..6),
        b_terms in proptest::collection::vec((0usize..1024usize, 0u64..16u64), 0..6),
        seed in any::<u64>()
    ) {
        let params = compact_dbfv().unwrap();
        let p = params.plain_modulus;
        let n = params.bfv_params.ring_degree;
        let mut rng = ChaCha20Rng::seed_from_u64(seed);
        let sk = gen_secret_key_with_rng(&params.bfv_params, &mut rng).unwrap();

        let pt_a = terms_to_poly(&a_terms, n, p);
        let pt_b = terms_to_poly(&b_terms, n, p);
        let expected = pt_a.add(&pt_b).unwrap();

        let ct_a = dbfv_encrypt_poly_sk_with_rng(&pt_a, &sk, &params, &mut rng).unwrap();
        let ct_b = dbfv_encrypt_poly_sk_with_rng(&pt_b, &sk, &params, &mut rng).unwrap();
        let ct_sum = dbfv_add(&ct_a, &ct_b).unwrap();
        let dec_sum = dbfv_decrypt_poly(&ct_sum, &sk).unwrap();

        prop_assert_eq!(dec_sum.coeffs, expected.coeffs);
    }

    #[test]
    fn prop_dbfv_poly_mul_sparse(
        a_terms in proptest::collection::vec((0usize..1024usize, 0u64..8u64), 0..5),
        b_terms in proptest::collection::vec((0usize..1024usize, 0u64..8u64), 0..5),
        seed in any::<u64>()
    ) {
        let params = compact_dbfv().unwrap();
        let p = params.plain_modulus;
        let n = params.bfv_params.ring_degree;
        let mut rng = ChaCha20Rng::seed_from_u64(seed);
        let sk = gen_secret_key_with_rng(&params.bfv_params, &mut rng).unwrap();
        let rlk = gen_relin_key_with_rng(&sk, &mut rng).unwrap();

        let pt_a = terms_to_poly(&a_terms, n, p);
        let pt_b = terms_to_poly(&b_terms, n, p);
        let expected = sparse_negacyclic_mul(&a_terms, &b_terms, n, p);

        let ct_a = dbfv_encrypt_poly_sk_with_rng(&pt_a, &sk, &params, &mut rng).unwrap();
        let ct_b = dbfv_encrypt_poly_sk_with_rng(&pt_b, &sk, &params, &mut rng).unwrap();
        let ct_prod = dbfv_mul(&ct_a, &ct_b, &rlk).unwrap();
        let dec_prod = dbfv_decrypt_poly(&ct_prod, &sk).unwrap();

        prop_assert_eq!(dec_prod.coeffs, expected);
    }
}
