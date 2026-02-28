use std::sync::Arc;
use concrete_ntt::prime64::Plan;

use crate::error::{ExactoError, Result};
use crate::ring::modular::{mod_mul, mod_add, mod_inv, barrett_constant};
use crate::ring::ntt::{NttPoly, make_plan};
use crate::ring::poly::CoeffPoly;

/// Polynomial in RNS (Residue Number System) representation.
///
/// Stores one NttPoly per RNS prime. The actual polynomial lives in Z_Q[X]/(X^n+1)
/// where Q = ∏ q_i.
#[derive(Clone, Debug)]
pub struct RnsPoly {
    pub components: Vec<NttPoly>,
    pub ring_degree: usize,
}

/// Stores precomputed data for an RNS basis.
#[derive(Clone, Debug)]
pub struct RnsBasis {
    pub moduli: Vec<u64>,
    pub plans: Vec<Arc<Plan>>,
    pub ring_degree: usize,
    /// Barrett constants for each modulus.
    pub barrett_ks: Vec<u64>,
    /// q_star_i = Q / q_i for each i.
    /// Stored as their residues mod each q_j (for CRT reconstruction).
    /// q_star_inv_i = (Q/q_i)^{-1} mod q_i.
    pub q_star_inv: Vec<u64>,
}

impl RnsBasis {
    /// Create a new RNS basis from a list of NTT-friendly primes.
    pub fn new(moduli: Vec<u64>, ring_degree: usize) -> Result<Self> {
        let plans: Vec<Arc<Plan>> = moduli.iter()
            .map(|&q| make_plan(ring_degree, q))
            .collect::<Result<Vec<_>>>()?;

        let barrett_ks: Vec<u64> = moduli.iter()
            .map(|&q| barrett_constant(q))
            .collect();

        // Compute q_star_inv_i = (Q/q_i)^{-1} mod q_i
        let q_star_inv: Vec<u64> = (0..moduli.len()).map(|i| {
            let mut prod = 1u64;
            let bk = barrett_ks[i];
            for (j, &qj) in moduli.iter().enumerate() {
                if i != j {
                    prod = mod_mul(prod, qj % moduli[i], moduli[i], bk);
                }
            }
            mod_inv(prod, moduli[i]).expect("RNS moduli must be coprime")
        }).collect();

        Ok(Self {
            moduli,
            plans,
            ring_degree,
            barrett_ks,
            q_star_inv,
        })
    }

    pub fn num_moduli(&self) -> usize {
        self.moduli.len()
    }
}

impl RnsPoly {
    /// Create a zero polynomial in RNS.
    pub fn zero(basis: &RnsBasis) -> Self {
        let components = basis.moduli.iter()
            .zip(basis.plans.iter())
            .map(|(&q, plan)| NttPoly::zero(basis.ring_degree, q, plan.clone()))
            .collect();
        Self {
            components,
            ring_degree: basis.ring_degree,
        }
    }

    /// Create an RnsPoly from a single CoeffPoly by reducing mod each RNS prime and NTT-ing.
    pub fn from_coeff_poly(poly: &CoeffPoly, basis: &RnsBasis) -> Result<Self> {
        if poly.len() != basis.ring_degree {
            return Err(ExactoError::DimensionMismatch {
                expected: basis.ring_degree,
                got: poly.len(),
            });
        }
        let components: Vec<NttPoly> = basis.moduli.iter()
            .zip(basis.plans.iter())
            .map(|(&q, plan)| {
                let reduced = CoeffPoly::from_coeffs(
                    poly.coeffs.iter().map(|&c| c % q).collect(),
                    q,
                );
                NttPoly::from_coeff_poly(&reduced, plan.clone())
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(Self {
            components,
            ring_degree: basis.ring_degree,
        })
    }

    /// Convert back to a CoeffPoly using CRT reconstruction.
    /// Result is in Z_Q where Q = product of all moduli.
    /// Returns coefficients in [0, Q).
    ///
    /// NOTE: For large Q (multi-limb), this returns approximate coefficients mod the largest
    /// representable value in u64. For full-precision CRT, use `to_coeff_poly_bigint`.
    /// This version works correctly when the true coefficients fit in u64.
    pub fn to_coeff_poly(&self, basis: &RnsBasis) -> CoeffPoly {
        let n = self.ring_degree;
        // First, inverse NTT each component
        let coeff_components: Vec<CoeffPoly> = self.components.iter()
            .map(|ntt| ntt.to_coeff_poly())
            .collect();

        // CRT: for each coefficient position j,
        // x_j = Σ_i (c_{i,j} * q_star_inv_i mod q_i) * (Q/q_i)  mod Q
        // When Q fits in u128, we can do this directly.
        // For now, compute in u128 and reduce mod first modulus (simple version).

        // Actually, let's just return the representation mod each prime.
        // For BFV usage, we typically keep things in RNS and only reconstruct when needed.
        // For a simple reconstruction assuming Q fits in u128:
        let num_mods = basis.moduli.len();
        if num_mods == 1 {
            return coeff_components[0].clone();
        }

        // Compute Q as u128 (works for up to 2 primes with 64-bit moduli)
        let big_q: u128 = basis.moduli.iter().map(|&q| q as u128).product();

        let mut coeffs = vec![0u64; n];
        for j in 0..n {
            let mut val = 0u128;
            for i in 0..num_mods {
                let c_ij = coeff_components[i].coeffs[j];
                let t = mod_mul(c_ij, basis.q_star_inv[i], basis.moduli[i], basis.barrett_ks[i]);
                // Q/q_i
                let q_star_i: u128 = big_q / basis.moduli[i] as u128;
                val = (val + t as u128 * q_star_i) % big_q;
            }
            coeffs[j] = val as u64;
        }

        CoeffPoly { coeffs, modulus: big_q as u64 }
    }

    /// Number of RNS components.
    pub fn num_components(&self) -> usize {
        self.components.len()
    }

    /// Component-wise addition in RNS.
    pub fn add(&self, other: &Self) -> Result<Self> {
        if self.components.len() != other.components.len() {
            return Err(ExactoError::DimensionMismatch {
                expected: self.components.len(),
                got: other.components.len(),
            });
        }
        let components = self.components.iter()
            .zip(other.components.iter())
            .map(|(a, b)| a.add(b))
            .collect::<Result<Vec<_>>>()?;
        Ok(Self { components, ring_degree: self.ring_degree })
    }

    /// Component-wise subtraction in RNS.
    pub fn sub(&self, other: &Self) -> Result<Self> {
        if self.components.len() != other.components.len() {
            return Err(ExactoError::DimensionMismatch {
                expected: self.components.len(),
                got: other.components.len(),
            });
        }
        let components = self.components.iter()
            .zip(other.components.iter())
            .map(|(a, b)| a.sub(b))
            .collect::<Result<Vec<_>>>()?;
        Ok(Self { components, ring_degree: self.ring_degree })
    }

    /// Negate.
    pub fn neg(&self) -> Self {
        let components = self.components.iter()
            .map(|a| a.neg())
            .collect();
        Self { components, ring_degree: self.ring_degree }
    }

    /// Component-wise multiplication in RNS (= polynomial multiplication).
    pub fn mul(&self, other: &Self) -> Result<Self> {
        if self.components.len() != other.components.len() {
            return Err(ExactoError::DimensionMismatch {
                expected: self.components.len(),
                got: other.components.len(),
            });
        }
        let components = self.components.iter()
            .zip(other.components.iter())
            .map(|(a, b)| a.mul(b))
            .collect::<Result<Vec<_>>>()?;
        Ok(Self { components, ring_degree: self.ring_degree })
    }

    /// Scalar multiply.
    pub fn scalar_mul(&self, scalar: u64) -> Self {
        let components = self.components.iter()
            .map(|c| c.scalar_mul(scalar))
            .collect();
        Self { components, ring_degree: self.ring_degree }
    }
}

/// Fast base extension from basis Q to basis P.
/// Given x mod q_i for each i, compute x mod p_j for each j,
/// using the approach from Bajard et al.
///
/// This is approximate: it computes x + α·Q where α ∈ {0, 1} is a correction term.
/// For BFV multiplication this approximation suffices.
pub fn fast_base_extend(
    poly: &RnsPoly,
    from: &RnsBasis,
    to: &RnsBasis,
) -> Result<RnsPoly> {
    let n = poly.ring_degree;

    // Step 1: compute t_i = c_i * q_star_inv_i mod q_i (in coeff domain)
    let coeff_components: Vec<CoeffPoly> = poly.components.iter()
        .map(|ntt| ntt.to_coeff_poly())
        .collect();

    let t_coeffs: Vec<Vec<u64>> = coeff_components.iter()
        .enumerate()
        .map(|(i, cp)| {
            let bk = from.barrett_ks[i];
            cp.coeffs.iter()
                .map(|&c| mod_mul(c, from.q_star_inv[i], from.moduli[i], bk))
                .collect()
        })
        .collect();

    // Step 2: for each target modulus p_j, compute
    //   x mod p_j ≈ Σ_i t_i * (Q/q_i) mod p_j
    let mut target_components = Vec::with_capacity(to.moduli.len());
    for (j, (&pj, plan)) in to.moduli.iter().zip(to.plans.iter()).enumerate() {
        let bk = to.barrett_ks[j];
        // Precompute (Q/q_i) mod p_j for each i
        let q_star_mod_pj: Vec<u64> = (0..from.moduli.len()).map(|i| {
            let mut prod = 1u64;
            for (k, &qk) in from.moduli.iter().enumerate() {
                if k != i {
                    prod = mod_mul(prod, qk % pj, pj, bk);
                }
            }
            prod
        }).collect();

        let mut coeffs = vec![0u64; n];
        for pos in 0..n {
            let mut val = 0u64;
            for i in 0..from.moduli.len() {
                let contrib = mod_mul(t_coeffs[i][pos], q_star_mod_pj[i], pj, bk);
                val = mod_add(val, contrib, pj);
            }
            coeffs[pos] = val;
        }

        let cp = CoeffPoly { coeffs, modulus: pj };
        target_components.push(NttPoly::from_coeff_poly(&cp, plan.clone())?);
    }

    Ok(RnsPoly {
        components: target_components,
        ring_degree: n,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    // concrete-ntt minimum polynomial_size is 16
    const TEST_N: usize = 16;
    const TEST_Q: u64 = 65537; // prime, ≡ 1 (mod 32)

    fn make_vec(v: &[u64]) -> Vec<u64> {
        let mut r = vec![0u64; TEST_N];
        r[..v.len()].copy_from_slice(v);
        r
    }

    #[test]
    fn test_rns_roundtrip_single() {
        let basis = RnsBasis::new(vec![TEST_Q], TEST_N).unwrap();

        let original = CoeffPoly::from_coeffs(make_vec(&[1, 2, 3, 4, 5, 6, 7, 8]), TEST_Q);
        let rns = RnsPoly::from_coeff_poly(&original, &basis).unwrap();
        let recovered = rns.to_coeff_poly(&basis);
        assert_eq!(original.coeffs, recovered.coeffs);
    }

    #[test]
    fn test_rns_add() {
        let basis = RnsBasis::new(vec![TEST_Q], TEST_N).unwrap();

        let a = CoeffPoly::from_coeffs(make_vec(&[1, 2, 3, 4, 5, 6, 7, 8]), TEST_Q);
        let b = CoeffPoly::from_coeffs(make_vec(&[10, 20, 30, 40, 50, 60, 70, 80]), TEST_Q);

        let ra = RnsPoly::from_coeff_poly(&a, &basis).unwrap();
        let rb = RnsPoly::from_coeff_poly(&b, &basis).unwrap();
        let rc = ra.add(&rb).unwrap();

        let result = rc.to_coeff_poly(&basis);
        let expected = a.add(&b).unwrap();
        assert_eq!(result.coeffs, expected.coeffs);
    }

    #[test]
    fn test_rns_mul() {
        let basis = RnsBasis::new(vec![TEST_Q], TEST_N).unwrap();

        let a = CoeffPoly::from_coeffs(make_vec(&[1, 1]), TEST_Q);
        let b = CoeffPoly::from_coeffs(make_vec(&[1, 1]), TEST_Q);

        let ra = RnsPoly::from_coeff_poly(&a, &basis).unwrap();
        let rb = RnsPoly::from_coeff_poly(&b, &basis).unwrap();
        let rc = ra.mul(&rb).unwrap();

        let result = rc.to_coeff_poly(&basis);
        let expected = a.mul_naive(&b).unwrap();
        assert_eq!(result.coeffs, expected.coeffs);
    }
}
