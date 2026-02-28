use crate::error::{ExactoError, Result};
use crate::ring::modular::{mod_add, mod_sub, mod_neg, mod_mul, barrett_constant};

/// Polynomial in coefficient representation over Z_q[X]/(X^n + 1).
#[derive(Clone, Debug)]
pub struct CoeffPoly {
    pub coeffs: Vec<u64>,
    pub modulus: u64,
}

impl CoeffPoly {
    /// Create a zero polynomial of degree < n in Z_q.
    pub fn zero(n: usize, modulus: u64) -> Self {
        Self {
            coeffs: vec![0u64; n],
            modulus,
        }
    }

    /// Create a polynomial from coefficients (reduced mod q).
    pub fn from_coeffs(coeffs: Vec<u64>, modulus: u64) -> Self {
        let mut p = Self { coeffs, modulus };
        p.reduce();
        p
    }

    /// Ring degree (number of coefficients).
    pub fn len(&self) -> usize {
        self.coeffs.len()
    }

    /// Reduce all coefficients mod q.
    pub fn reduce(&mut self) {
        for c in self.coeffs.iter_mut() {
            *c %= self.modulus;
        }
    }

    /// Add two polynomials in Z_q[X]/(X^n+1).
    pub fn add(&self, other: &Self) -> Result<Self> {
        if self.len() != other.len() {
            return Err(ExactoError::DimensionMismatch {
                expected: self.len(),
                got: other.len(),
            });
        }
        if self.modulus != other.modulus {
            return Err(ExactoError::ModulusMismatch);
        }
        let coeffs = self.coeffs.iter()
            .zip(other.coeffs.iter())
            .map(|(&a, &b)| mod_add(a, b, self.modulus))
            .collect();
        Ok(Self { coeffs, modulus: self.modulus })
    }

    /// Subtract two polynomials in Z_q[X]/(X^n+1).
    pub fn sub(&self, other: &Self) -> Result<Self> {
        if self.len() != other.len() {
            return Err(ExactoError::DimensionMismatch {
                expected: self.len(),
                got: other.len(),
            });
        }
        if self.modulus != other.modulus {
            return Err(ExactoError::ModulusMismatch);
        }
        let coeffs = self.coeffs.iter()
            .zip(other.coeffs.iter())
            .map(|(&a, &b)| mod_sub(a, b, self.modulus))
            .collect();
        Ok(Self { coeffs, modulus: self.modulus })
    }

    /// Negate polynomial.
    pub fn neg(&self) -> Self {
        let coeffs = self.coeffs.iter()
            .map(|&a| mod_neg(a, self.modulus))
            .collect();
        Self { coeffs, modulus: self.modulus }
    }

    /// Schoolbook multiply in Z_q[X]/(X^n+1).
    /// For testing/small polys only — use NTT for production.
    pub fn mul_naive(&self, other: &Self) -> Result<Self> {
        if self.len() != other.len() {
            return Err(ExactoError::DimensionMismatch {
                expected: self.len(),
                got: other.len(),
            });
        }
        if self.modulus != other.modulus {
            return Err(ExactoError::ModulusMismatch);
        }
        let n = self.len();
        let bk = barrett_constant(self.modulus);
        let mut result = vec![0u64; n];

        for i in 0..n {
            if self.coeffs[i] == 0 {
                continue;
            }
            for j in 0..n {
                if other.coeffs[j] == 0 {
                    continue;
                }
                let prod = mod_mul(self.coeffs[i], other.coeffs[j], self.modulus, bk);
                let idx = i + j;
                if idx < n {
                    result[idx] = mod_add(result[idx], prod, self.modulus);
                } else {
                    // X^n ≡ -1 in X^n+1, so wrap with negation
                    let idx = idx - n;
                    result[idx] = mod_sub(result[idx], prod, self.modulus);
                }
            }
        }

        Ok(Self { coeffs: result, modulus: self.modulus })
    }

    /// Multiply by a scalar.
    pub fn scalar_mul(&self, scalar: u64) -> Self {
        let bk = barrett_constant(self.modulus);
        let s = scalar % self.modulus;
        let coeffs = self.coeffs.iter()
            .map(|&c| mod_mul(c, s, self.modulus, bk))
            .collect();
        Self { coeffs, modulus: self.modulus }
    }

    /// Check if all coefficients are zero.
    pub fn is_zero(&self) -> bool {
        self.coeffs.iter().all(|&c| c == 0)
    }

    /// Get centered representation: map [0, q) -> [-(q-1)/2, (q-1)/2]
    pub fn centered_coeffs(&self) -> Vec<i64> {
        let half = self.modulus / 2;
        self.coeffs.iter().map(|&c| {
            if c > half {
                c as i64 - self.modulus as i64
            } else {
                c as i64
            }
        }).collect()
    }
}

impl PartialEq for CoeffPoly {
    fn eq(&self, other: &Self) -> bool {
        self.modulus == other.modulus && self.coeffs == other.coeffs
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zero() {
        let p = CoeffPoly::zero(4, 17);
        assert!(p.is_zero());
        assert_eq!(p.len(), 4);
    }

    #[test]
    fn test_add_sub() {
        let a = CoeffPoly::from_coeffs(vec![1, 2, 3, 4], 17);
        let b = CoeffPoly::from_coeffs(vec![5, 6, 7, 8], 17);
        let c = a.add(&b).unwrap();
        assert_eq!(c.coeffs, vec![6, 8, 10, 12]);
        let d = c.sub(&b).unwrap();
        assert_eq!(d.coeffs, a.coeffs);
    }

    #[test]
    fn test_neg() {
        let a = CoeffPoly::from_coeffs(vec![1, 0, 3, 16], 17);
        let neg_a = a.neg();
        let zero = a.add(&neg_a).unwrap();
        assert!(zero.is_zero());
    }

    #[test]
    fn test_mul_naive() {
        // In Z_17[X]/(X^4+1):
        // (1 + X) * (1 + X) = 1 + 2X + X^2
        let a = CoeffPoly::from_coeffs(vec![1, 1, 0, 0], 17);
        let c = a.mul_naive(&a).unwrap();
        assert_eq!(c.coeffs, vec![1, 2, 1, 0]);
    }

    #[test]
    fn test_mul_naive_wraparound() {
        // In Z_17[X]/(X^4+1):
        // X^3 * X^3 = X^6 = X^4 * X^2 = -X^2 (mod X^4+1)
        let a = CoeffPoly::from_coeffs(vec![0, 0, 0, 1], 17); // X^3
        let c = a.mul_naive(&a).unwrap();
        // Should be -X^2 = [0, 0, 17-1, 0] = [0, 0, 16, 0]
        assert_eq!(c.coeffs, vec![0, 0, 16, 0]);
    }

    #[test]
    fn test_scalar_mul() {
        let a = CoeffPoly::from_coeffs(vec![1, 2, 3, 4], 17);
        let c = a.scalar_mul(3);
        assert_eq!(c.coeffs, vec![3, 6, 9, 12]);
    }

    #[test]
    fn test_centered() {
        // q=17, half=8. Values > 8 map to negative.
        // 0→0, 1→1, 16→-1, 9→-8
        let a = CoeffPoly::from_coeffs(vec![0, 1, 16, 9], 17);
        let c = a.centered_coeffs();
        assert_eq!(c, vec![0, 1, -1, -8]);
    }
}
