use std::sync::Arc;

use crate::error::Result;
use crate::ring::modular::{mod_mul, mod_add, mod_inv, mod_neg, barrett_constant};
use crate::ring::poly::CoeffPoly;
use crate::ring::rns::RnsPoly;
use crate::bfv::{BfvCiphertext, keygen::RelinKey};
use crate::bfv::encrypt::scale_plaintext;
use crate::bfv::eval::{bfv_add, bfv_mul_and_relin, bfv_plain_mul};
use crate::params::BfvParams;

/// Compute the rounding polynomial used in BFV bootstrapping.
///
/// Defines g(x) = round(t_orig · (x mod q') / q') mod t_orig for all
/// x ∈ [0, t_boot). The polynomial is periodic with period q' and is
/// interpolated over Z_{t_boot} via Lagrange.
///
/// Returns polynomial coefficients in Z_{t_boot}[x] (degree t_boot - 1).
pub fn compute_rounding_poly(t_orig: u64, q_prime: u64, t_boot: u64) -> Vec<u64> {
    let values: Vec<u64> = (0..t_boot)
        .map(|x| {
            let reduced = x % q_prime;
            // round(t_orig * reduced / q_prime) mod t_orig
            let scaled = (t_orig as u128 * reduced as u128 + q_prime as u128 / 2) / q_prime as u128;
            (scaled % t_orig as u128) as u64
        })
        .collect();
    lagrange_interpolate(&values, t_boot)
}

/// Lagrange interpolation at consecutive integer points 0, 1, ..., n-1.
///
/// Given values[i] = f(i) for i = 0,...,n-1, returns polynomial coefficients
/// [c_0, c_1, ..., c_{n-1}] such that f(x) = Σ c_k * x^k in Z_p[x].
///
/// Requires p prime (for modular inverses).
pub fn lagrange_interpolate(values: &[u64], p: u64) -> Vec<u64> {
    let n = values.len();
    if n == 0 { return vec![]; }
    if n == 1 { return vec![values[0] % p]; }

    let bk = barrett_constant(p);
    let mut result = vec![0u64; n];

    for j in 0..n {
        if values[j] % p == 0 { continue; }

        // Compute Lagrange basis L_j(x) = Π_{k≠j} (x - k) / (j - k)
        // Build numerator polynomial: Π_{k≠j} (x - k)
        let mut num = vec![0u64; n];
        num[0] = 1;
        let mut deg = 0;

        for k in 0..n {
            if k == j { continue; }
            let neg_k = mod_neg((k as u64) % p, p);
            let mut new_num = vec![0u64; n];
            for d in 0..=deg {
                if d + 1 < n {
                    new_num[d + 1] = mod_add(new_num[d + 1], num[d], p);
                }
                new_num[d] = mod_add(new_num[d], mod_mul(num[d], neg_k, p, bk), p);
            }
            num = new_num;
            deg += 1;
        }

        // Compute denominator: Π_{k≠j} (j - k)
        let mut denom = 1u64;
        for k in 0..n {
            if k == j { continue; }
            let diff = if j >= k {
                (j - k) as u64 % p
            } else {
                p - ((k - j) as u64 % p)
            };
            denom = mod_mul(denom, diff, p, bk);
        }
        let denom_inv = mod_inv(denom, p).expect("points must be distinct mod p");

        // Scale by values[j] / denominator
        let scale = mod_mul(values[j] % p, denom_inv, p, bk);

        for d in 0..n {
            result[d] = mod_add(result[d], mod_mul(num[d], scale, p, bk), p);
        }
    }

    result
}

/// Evaluate a polynomial homomorphically on an encrypted value using
/// Paterson-Stockmeyer (baby-step/giant-step) algorithm.
///
/// Given ct encrypting x ∈ Z_t and polynomial f of degree d,
/// computes ct' encrypting f(x).
///
/// Baby-step size k = ceil(sqrt(d+1)) gives O(sqrt(d)) multiplications
/// at multiplicative depth O(log(k) + d/k).
pub fn eval_poly_homomorphic(
    ct_x: &BfvCiphertext,
    poly_coeffs: &[u64],
    rlk: &RelinKey,
) -> Result<BfvCiphertext> {
    let params = &ct_x.params;
    let d = poly_coeffs.len().saturating_sub(1);

    if d == 0 {
        return trivial_encrypt(poly_coeffs[0], params);
    }

    // Choose baby-step size k
    let k = ((d as f64 + 1.0).sqrt().ceil() as usize).max(2);

    // Baby steps: compute x^1, x^2, ..., x^k
    let mut baby = Vec::with_capacity(k + 1);
    baby.push(trivial_encrypt(1, params)?); // x^0 = 1
    baby.push(ct_x.clone()); // x^1
    for i in 2..=k {
        // Use optimal multiplication tree
        let half = i / 2;
        let other = i - half;
        let prod = bfv_mul_and_relin(&baby[half], &baby[other], rlk)?;
        baby.push(prod);
    }

    // Compute x^k (already in baby[k])
    // Split polynomial into giant-step groups: f(x) = g_0(x) + x^k*g_1(x) + x^{2k}*g_2(x) + ...
    let num_groups = (d + k) / k; // ceil((d+1)/k)

    // Evaluate each group g_i(x) = Σ_{j=0}^{k-1} a_{i*k+j} * x^j
    let mut groups: Vec<BfvCiphertext> = Vec::with_capacity(num_groups);
    for i in 0..num_groups {
        let start = i * k;
        let mut group_ct = trivial_encrypt(0, params)?;
        for j in 0..k {
            let idx = start + j;
            if idx >= poly_coeffs.len() { break; }
            let coeff = poly_coeffs[idx];
            if coeff == 0 { continue; }
            // coeff * x^j
            let term = bfv_scalar_mul(&baby[j], coeff, params)?;
            group_ct = bfv_add(&group_ct, &term)?;
        }
        groups.push(group_ct);
    }

    // Giant-step Horner: result = g_{m-1}, then result = result * x^k + g_{m-2}, ...
    let mut result = groups.pop().unwrap();
    let x_k = &baby[k];
    while let Some(g) = groups.pop() {
        result = bfv_mul_and_relin(&result, x_k, rlk)?;
        result = bfv_add(&result, &g)?;
    }

    Ok(result)
}

/// Create a trivial encryption of a scalar m (zero noise).
/// ct = (Δ*m, 0) where Δ = floor(q/t).
pub fn trivial_encrypt(m: u64, params: &Arc<BfvParams>) -> Result<BfvCiphertext> {
    let basis = &params.ct_basis;
    let mut pt_coeffs = vec![0u64; params.ring_degree];
    pt_coeffs[0] = m % params.plain_modulus;
    let pt = CoeffPoly {
        coeffs: pt_coeffs,
        modulus: params.plain_modulus,
    };

    let c0 = scale_plaintext(&pt, params)?;
    let c1 = RnsPoly::zero(basis);

    Ok(BfvCiphertext {
        c: vec![c0, c1],
        params: params.clone(),
    })
}

/// Trivial encryption of a polynomial (zero noise).
pub fn trivial_encrypt_poly(poly: &CoeffPoly, params: &Arc<BfvParams>) -> Result<BfvCiphertext> {
    let basis = &params.ct_basis;
    let c0 = scale_plaintext(poly, params)?;
    let c1 = RnsPoly::zero(basis);

    Ok(BfvCiphertext {
        c: vec![c0, c1],
        params: params.clone(),
    })
}

/// Multiply a ciphertext by a scalar plaintext value.
fn bfv_scalar_mul(ct: &BfvCiphertext, scalar: u64, params: &Arc<BfvParams>) -> Result<BfvCiphertext> {
    let mut coeffs = vec![0u64; params.ring_degree];
    coeffs[0] = scalar % params.plain_modulus;
    let pt = CoeffPoly { coeffs, modulus: params.plain_modulus };
    bfv_plain_mul(ct, &pt)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ring::modular::{mod_mul, barrett_constant};

    #[test]
    fn test_lagrange_simple() {
        // f(x) = x over Z_7: f(0)=0, f(1)=1, f(2)=2
        let values = vec![0, 1, 2];
        let coeffs = lagrange_interpolate(&values, 7);
        // Should give [0, 1, 0] (i.e., f(x) = x)
        assert_eq!(coeffs[0], 0);
        assert_eq!(coeffs[1], 1);
        assert_eq!(coeffs[2], 0);
    }

    #[test]
    fn test_lagrange_quadratic() {
        // f(x) = x^2 over Z_7: f(0)=0, f(1)=1, f(2)=4, f(3)=2 (9 mod 7)
        let values = vec![0, 1, 4, 2];
        let coeffs = lagrange_interpolate(&values, 7);
        // Should give [0, 0, 1, 0] (f(x) = x^2)
        assert_eq!(coeffs[0], 0);
        assert_eq!(coeffs[1], 0);
        assert_eq!(coeffs[2], 1);
        assert_eq!(coeffs[3], 0);
    }

    #[test]
    fn test_lagrange_eval() {
        // Interpolate random values, then evaluate to verify
        let p = 29u64;
        let values: Vec<u64> = (0..10).map(|i| (i * i + 3 * i + 7) % p).collect();
        let coeffs = lagrange_interpolate(&values, p);
        let bk = barrett_constant(p);

        for (x, &expected) in values.iter().enumerate() {
            let mut result = 0u64;
            let mut x_pow = 1u64;
            for &c in &coeffs {
                result = mod_add(result, mod_mul(c, x_pow, p, bk), p);
                x_pow = mod_mul(x_pow, x as u64, p, bk);
            }
            assert_eq!(result, expected, "mismatch at x={x}");
        }
    }

    #[test]
    fn test_rounding_poly() {
        // t=5, q'=25, t_boot=29 (bootstrap test parameters)
        let t = 5u64;
        let q_prime = 25u64;
        let t_boot = 29u64;
        let poly = compute_rounding_poly(t, q_prime, t_boot);
        let bk = barrett_constant(t_boot);

        // Verify for ALL x in [0, t_boot): g(x) = f(x mod q')
        for x in 0..t_boot {
            let reduced = x % q_prime;
            let expected = ((t as u128 * reduced as u128 + q_prime as u128 / 2)
                / q_prime as u128 % t as u128) as u64;
            let mut result = 0u64;
            let mut x_pow = 1u64;
            for &c in &poly {
                result = mod_add(result, mod_mul(c, x_pow, t_boot, bk), t_boot);
                x_pow = mod_mul(x_pow, x % t_boot, t_boot, bk);
            }
            assert_eq!(result % t, expected, "mismatch at x={x}: got {result}, expected {expected}");
        }
    }

    #[test]
    fn test_trivial_encrypt_decrypt() {
        use crate::bfv::keygen::*;
        use crate::bfv::encrypt::decrypt;
        use crate::bfv::encoding::decode_scalar;
        use crate::params::presets::compact_bfv;
        use rand::SeedableRng;
        use rand_chacha::ChaCha20Rng;

        let params = compact_bfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let sk = gen_secret_key_with_rng(&params, &mut rng).unwrap();

        for m in [0u64, 1, 42, 100, 256] {
            let ct = trivial_encrypt(m, &params).unwrap();
            let decrypted = decrypt(&ct, &sk).unwrap();
            assert_eq!(decode_scalar(&decrypted), m, "failed for m={m}");
        }
    }
}
