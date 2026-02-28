use std::sync::Arc;
use concrete_ntt::prime64::Plan;

use crate::error::{ExactoError, Result};
use crate::ring::poly::CoeffPoly;

/// Polynomial in NTT (evaluation) representation over Z_q[X]/(X^n + 1).
///
/// Uses `concrete-ntt` for hardware-accelerated NTT (AVX2/AVX-512/NEON).
#[derive(Clone, Debug)]
pub struct NttPoly {
    pub evals: Vec<u64>,
    pub modulus: u64,
    pub plan: Arc<Plan>,
}

/// Cache of NTT plans keyed by (ring_degree, modulus).
/// In production you'd use a concurrent map; for now callers share plans via Arc.
pub fn make_plan(n: usize, modulus: u64) -> Result<Arc<Plan>> {
    if !n.is_power_of_two() || n < 2 {
        return Err(ExactoError::InvalidRingDegree(n));
    }
    // concrete-ntt requires modulus to be prime and ≡ 1 (mod 2n)
    let plan = Plan::try_new(n, modulus)
        .ok_or_else(|| ExactoError::InvalidParam(
            format!("cannot create NTT plan for n={n}, q={modulus} (need prime q ≡ 1 mod {}", 2 * n)
        ))?;
    Ok(Arc::new(plan))
}

impl NttPoly {
    /// Create a zero polynomial in NTT domain.
    pub fn zero(n: usize, modulus: u64, plan: Arc<Plan>) -> Self {
        Self {
            evals: vec![0u64; n],
            modulus,
            plan,
        }
    }

    /// Forward NTT: convert from coefficient to evaluation representation.
    pub fn from_coeff_poly(poly: &CoeffPoly, plan: Arc<Plan>) -> Result<Self> {
        if poly.modulus != plan.modulus() {
            return Err(ExactoError::ModulusMismatch);
        }
        let _n = poly.len();
        let mut evals = poly.coeffs.clone();
        // concrete-ntt operates in-place
        plan.fwd(&mut evals);
        Ok(Self {
            evals,
            modulus: poly.modulus,
            plan,
        })
    }

    /// Inverse NTT: convert back to coefficient representation.
    pub fn to_coeff_poly(&self) -> CoeffPoly {
        let mut coeffs = self.evals.clone();
        self.plan.inv(&mut coeffs);
        // concrete-ntt's inv already normalizes by 1/n
        self.plan.normalize(&mut coeffs);
        CoeffPoly {
            coeffs,
            modulus: self.modulus,
        }
    }

    /// Ring degree.
    pub fn len(&self) -> usize {
        self.evals.len()
    }

    /// Component-wise addition in NTT domain (= polynomial addition).
    pub fn add(&self, other: &Self) -> Result<Self> {
        if self.len() != other.len() || self.modulus != other.modulus {
            return Err(ExactoError::ModulusMismatch);
        }
        let mut evals = self.evals.clone();
        for (a, &b) in evals.iter_mut().zip(other.evals.iter()) {
            let sum = *a as u128 + b as u128;
            *a = if sum >= self.modulus as u128 {
                (sum - self.modulus as u128) as u64
            } else {
                sum as u64
            };
        }
        Ok(Self { evals, modulus: self.modulus, plan: self.plan.clone() })
    }

    /// Component-wise subtraction.
    pub fn sub(&self, other: &Self) -> Result<Self> {
        if self.len() != other.len() || self.modulus != other.modulus {
            return Err(ExactoError::ModulusMismatch);
        }
        let mut evals = self.evals.clone();
        for (a, &b) in evals.iter_mut().zip(other.evals.iter()) {
            if *a >= b {
                *a -= b;
            } else {
                *a = self.modulus - b + *a;
            }
        }
        Ok(Self { evals, modulus: self.modulus, plan: self.plan.clone() })
    }

    /// Component-wise negation.
    pub fn neg(&self) -> Self {
        let evals = self.evals.iter().map(|&a| {
            if a == 0 { 0 } else { self.modulus - a }
        }).collect();
        Self { evals, modulus: self.modulus, plan: self.plan.clone() }
    }

    /// Component-wise multiplication in NTT domain (= polynomial multiplication).
    ///
    /// This performs a plain pointwise multiply (without normalization).
    /// The normalization (1/n) is applied later during `to_coeff_poly` via `plan.normalize()`.
    pub fn mul(&self, other: &Self) -> Result<Self> {
        if self.len() != other.len() || self.modulus != other.modulus {
            return Err(ExactoError::ModulusMismatch);
        }
        let bk = crate::ring::modular::barrett_constant(self.modulus);
        let evals = self.evals.iter()
            .zip(other.evals.iter())
            .map(|(&a, &b)| crate::ring::modular::mod_mul(a, b, self.modulus, bk))
            .collect();
        Ok(Self { evals, modulus: self.modulus, plan: self.plan.clone() })
    }

    /// Multiply by scalar in NTT domain.
    pub fn scalar_mul(&self, scalar: u64) -> Self {
        let s = scalar % self.modulus;
        let bk = crate::ring::modular::barrett_constant(self.modulus);
        let evals = self.evals.iter().map(|&a| {
            crate::ring::modular::mod_mul(a, s, self.modulus, bk)
        }).collect();
        Self { evals, modulus: self.modulus, plan: self.plan.clone() }
    }

    /// Check if zero.
    pub fn is_zero(&self) -> bool {
        self.evals.iter().all(|&e| e == 0)
    }
}

impl PartialEq for NttPoly {
    fn eq(&self, other: &Self) -> bool {
        self.modulus == other.modulus && self.evals == other.evals
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // concrete-ntt requires polynomial_size >= 16, and modulus prime with p ≡ 1 (mod 2n)
    fn test_params() -> (usize, u64) {
        // n=16, q=65537 (Fermat prime, 2^16+1). 65537 ≡ 1 (mod 32). ✓
        (16, 65537)
    }

    fn make_vec(v: &[u64], n: usize) -> Vec<u64> {
        let mut r = vec![0u64; n];
        r[..v.len()].copy_from_slice(v);
        r
    }

    #[test]
    fn test_ntt_roundtrip() {
        let (n, q) = test_params();
        let plan = make_plan(n, q).unwrap();

        let original = CoeffPoly::from_coeffs(make_vec(&[1, 2, 3, 4, 5, 6, 7, 8], n), q);
        let ntt = NttPoly::from_coeff_poly(&original, plan).unwrap();
        let recovered = ntt.to_coeff_poly();
        assert_eq!(original.coeffs, recovered.coeffs);
    }

    #[test]
    fn test_ntt_mul_matches_naive() {
        let (n, q) = test_params();
        let plan = make_plan(n, q).unwrap();

        let a = CoeffPoly::from_coeffs(make_vec(&[1, 1], n), q);
        let b = CoeffPoly::from_coeffs(make_vec(&[1, 1], n), q);
        let naive_result = a.mul_naive(&b).unwrap();

        let a_ntt = NttPoly::from_coeff_poly(&a, plan.clone()).unwrap();
        let b_ntt = NttPoly::from_coeff_poly(&b, plan).unwrap();
        let c_ntt = a_ntt.mul(&b_ntt).unwrap();
        let ntt_result = c_ntt.to_coeff_poly();

        assert_eq!(naive_result.coeffs, ntt_result.coeffs);
    }

    #[test]
    fn test_ntt_add() {
        let (n, q) = test_params();
        let plan = make_plan(n, q).unwrap();

        let a = CoeffPoly::from_coeffs(make_vec(&[1, 2, 3], n), q);
        let b = CoeffPoly::from_coeffs(make_vec(&[4, 5, 6], n), q);

        let a_ntt = NttPoly::from_coeff_poly(&a, plan.clone()).unwrap();
        let b_ntt = NttPoly::from_coeff_poly(&b, plan).unwrap();
        let c_ntt = a_ntt.add(&b_ntt).unwrap();
        let c = c_ntt.to_coeff_poly();

        let expected = a.add(&b).unwrap();
        assert_eq!(c.coeffs, expected.coeffs);
    }
}
