use std::sync::Arc;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use zeroize::Zeroize;

use crate::error::Result;
use crate::params::BfvParams;
use crate::ring::rns::RnsPoly;
use crate::ring::poly::CoeffPoly;
use crate::sampling::{sample_gaussian_poly, sample_uniform_poly, sample_ternary_poly};

/// BFV secret key: s ∈ R_q (ternary polynomial, stored in NTT/RNS form).
pub struct SecretKey {
    /// s in RNS-NTT form.
    pub poly: RnsPoly,
    pub params: Arc<BfvParams>,
}

impl Drop for SecretKey {
    fn drop(&mut self) {
        // Zero out secret key material
        for comp in &mut self.poly.components {
            comp.evals.zeroize();
        }
    }
}

/// BFV public key: pk = (pk0, pk1) where pk0 = -(a·s + e), pk1 = a.
#[derive(Clone, Debug)]
pub struct PublicKey {
    pub pk0: RnsPoly,
    pub pk1: RnsPoly,
    pub params: Arc<BfvParams>,
}

/// Relinearization key for reducing degree-2 ciphertexts back to degree-1.
/// Uses gadget decomposition.
#[derive(Clone, Debug)]
pub struct RelinKey {
    /// Each entry (rlk0_i, rlk1_i) satisfies:
    /// rlk0_i + rlk1_i · s ≈ g_i · s^2 + e_i
    /// where g_i = base^i is the gadget component.
    pub keys: Vec<(RnsPoly, RnsPoly)>,
    pub params: Arc<BfvParams>,
}

/// Galois key for automorphism X → X^k (for rotations).
#[derive(Clone, Debug)]
pub struct GaloisKey {
    /// Key-switch key from s(X^k) to s(X).
    pub keys: Vec<(RnsPoly, RnsPoly)>,
    /// The Galois element k.
    pub element: usize,
    pub params: Arc<BfvParams>,
}

/// Generate a secret key (ternary distribution).
pub fn gen_secret_key(params: &Arc<BfvParams>) -> Result<SecretKey> {
    let mut rng = ChaCha20Rng::from_os_rng();
    gen_secret_key_with_rng(params, &mut rng)
}

/// Generate a secret key with a provided RNG.
pub fn gen_secret_key_with_rng<R: rand::Rng>(
    params: &Arc<BfvParams>,
    rng: &mut R,
) -> Result<SecretKey> {
    let basis = &params.ct_basis;
    // Sample ternary polynomial
    let s_coeffs = sample_ternary_poly(params.ring_degree, basis.moduli[0], rng);

    // Convert to RNS-NTT form
    let poly = RnsPoly::from_coeff_poly(&s_coeffs, basis)?;

    Ok(SecretKey {
        poly,
        params: params.clone(),
    })
}

/// Generate a public key from a secret key.
pub fn gen_public_key(sk: &SecretKey) -> Result<PublicKey> {
    let mut rng = ChaCha20Rng::from_os_rng();
    gen_public_key_with_rng(sk, &mut rng)
}

/// Generate a public key with a provided RNG.
pub fn gen_public_key_with_rng<R: rand::Rng>(
    sk: &SecretKey,
    rng: &mut R,
) -> Result<PublicKey> {
    let params = &sk.params;
    let basis = &params.ct_basis;

    // Sample uniform a ∈ R_q
    let a_coeffs = sample_uniform_poly(params.ring_degree, basis.moduli[0], rng);
    let a = RnsPoly::from_coeff_poly(&a_coeffs, basis)?;

    // Sample error e
    let e_coeffs = sample_gaussian_poly(params.ring_degree, basis.moduli[0], params.sigma, rng);
    let e = RnsPoly::from_coeff_poly(&e_coeffs, basis)?;

    // pk0 = -(a·s + e) mod q
    let a_s = a.mul(&sk.poly)?;
    let a_s_e = a_s.add(&e)?;
    let pk0 = a_s_e.neg();

    Ok(PublicKey {
        pk0,
        pk1: a,
        params: params.clone(),
    })
}

/// Generate a relinearization key.
/// rlk_i = (-(a_i·s + e_i) + g_i·s², a_i) where g_i = base^i.
pub fn gen_relin_key(sk: &SecretKey) -> Result<RelinKey> {
    let mut rng = ChaCha20Rng::from_os_rng();
    gen_relin_key_with_rng(sk, &mut rng)
}

/// Generate relinearization key with provided RNG.
pub fn gen_relin_key_with_rng<R: rand::Rng>(
    sk: &SecretKey,
    rng: &mut R,
) -> Result<RelinKey> {
    let params = &sk.params;
    let basis = &params.ct_basis;

    // Compute s^2
    let s_sq = sk.poly.mul(&sk.poly)?;

    let mut keys = Vec::with_capacity(params.gadget_digits);
    let mut gadget_s_sq = s_sq.clone(); // base^0 * s^2

    for i in 0..params.gadget_digits {
        // Sample uniform a_i
        let a_coeffs = sample_uniform_poly(params.ring_degree, basis.moduli[0], rng);
        let a = RnsPoly::from_coeff_poly(&a_coeffs, basis)?;

        // Sample error e_i
        let e_coeffs = sample_gaussian_poly(params.ring_degree, basis.moduli[0], params.sigma, rng);
        let e = RnsPoly::from_coeff_poly(&e_coeffs, basis)?;

        // rlk0_i = -(a_i·s + e_i) + g_i·s²
        let a_s = a.mul(&sk.poly)?;
        let a_s_e = a_s.add(&e)?;
        let neg_a_s_e = a_s_e.neg();
        let rlk0 = neg_a_s_e.add(&gadget_s_sq)?;

        keys.push((rlk0, a));

        if i + 1 < params.gadget_digits {
            gadget_s_sq = gadget_s_sq.scalar_mul(params.gadget_base);
        }
    }

    Ok(RelinKey {
        keys,
        params: params.clone(),
    })
}

/// Generate a Galois key for automorphism X → X^element.
pub fn gen_galois_key(sk: &SecretKey, element: usize) -> Result<GaloisKey> {
    let mut rng = ChaCha20Rng::from_os_rng();
    gen_galois_key_with_rng(sk, element, &mut rng)
}

/// Generate a Galois key with provided RNG.
pub fn gen_galois_key_with_rng<R: rand::Rng>(
    sk: &SecretKey,
    element: usize,
    rng: &mut R,
) -> Result<GaloisKey> {
    let params = &sk.params;
    let basis = &params.ct_basis;

    // Compute s(X^element): apply automorphism to secret key
    let s_coeffs = sk.poly.components[0].to_coeff_poly();
    let s_auto = apply_automorphism(&s_coeffs, element);
    let s_auto_rns = RnsPoly::from_coeff_poly(&s_auto, basis)?;

    // Generate key-switch key from s(X^element) to s(X)
    let mut keys = Vec::with_capacity(params.gadget_digits);
    let mut gadget_s_auto = s_auto_rns.clone(); // base^0 * s_auto

    for i in 0..params.gadget_digits {
        let a_coeffs = sample_uniform_poly(params.ring_degree, basis.moduli[0], rng);
        let a = RnsPoly::from_coeff_poly(&a_coeffs, basis)?;

        let e_coeffs = sample_gaussian_poly(params.ring_degree, basis.moduli[0], params.sigma, rng);
        let e = RnsPoly::from_coeff_poly(&e_coeffs, basis)?;

        let a_s = a.mul(&sk.poly)?;
        let a_s_e = a_s.add(&e)?;
        let neg_a_s_e = a_s_e.neg();
        let ks0 = neg_a_s_e.add(&gadget_s_auto)?;

        keys.push((ks0, a));
        if i + 1 < params.gadget_digits {
            gadget_s_auto = gadget_s_auto.scalar_mul(params.gadget_base);
        }
    }

    Ok(GaloisKey {
        keys,
        element,
        params: params.clone(),
    })
}

/// Apply the automorphism X → X^k to a polynomial in Z_q[X]/(X^n+1).
///
/// The map σ_k sends X^i to X^{ik} reduced modulo X^n+1. Since X^n = -1,
/// this is a signed permutation of the polynomial coefficients.
/// For k odd and coprime to 2n, σ_k is a ring automorphism.
pub fn apply_automorphism(poly: &CoeffPoly, k: usize) -> CoeffPoly {
    let n = poly.len();
    let q = poly.modulus;
    let mut result = vec![0u64; n];

    for (i, &c) in poly.coeffs.iter().enumerate() {
        if c == 0 {
            continue;
        }
        // X^i → X^{ik}, reduced modulo X^n + 1
        let new_exp = (i * k) % (2 * n);
        if new_exp < n {
            result[new_exp] = crate::ring::modular::mod_add(result[new_exp], c, q);
        } else {
            // X^{n+j} = -X^j in X^n+1
            let j = new_exp - n;
            result[j] = crate::ring::modular::mod_sub(result[j], c, q);
        }
    }

    CoeffPoly { coeffs: result, modulus: q }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::presets::compact_bfv;

    #[test]
    fn test_keygen() {
        let params = compact_bfv().unwrap();
        let sk = gen_secret_key(&params).unwrap();
        let pk = gen_public_key(&sk).unwrap();

        assert_eq!(sk.poly.ring_degree, 1024);
        assert_eq!(pk.pk0.ring_degree, 1024);
        assert_eq!(pk.pk1.ring_degree, 1024);
    }

    #[test]
    fn test_relin_keygen() {
        let params = compact_bfv().unwrap();
        let sk = gen_secret_key(&params).unwrap();
        let rlk = gen_relin_key(&sk).unwrap();

        assert!(!rlk.keys.is_empty());
    }

    #[test]
    fn test_automorphism() {
        // X → X^3 in Z_q[X]/(X^4+1)
        // X^0 → X^0, X^1 → X^3, X^2 → X^6 = -X^2, X^3 → X^9 = X^1 · X^8 = X · (X^4)^2 = X · 1 = X
        // Wait: X^9 mod (X^4+1): 9 mod 8 = 1, 9/8 = 1 odd → negate? Let me recompute.
        // new_exp = 3*3 = 9, 9 % (2*4) = 1, 1 < 4 so result[1] += coeff
        let p = CoeffPoly::from_coeffs(vec![1, 1, 1, 1], 17); // 1+X+X^2+X^3
        let r = apply_automorphism(&p, 3);
        // X^0 → X^0: result[0] += 1
        // X^1 → X^3: result[3] += 1
        // X^2 → X^6: 6 % 8 = 6, 6 >= 4, so result[2] -= 1 → 16 mod 17
        // X^3 → X^9: 9 % 8 = 1, 1 < 4, result[1] += 1
        assert_eq!(r.coeffs, vec![1, 1, 16, 1]);
    }
}
