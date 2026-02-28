use num_bigint::BigInt;
use num_traits::{Zero, One};

use crate::error::{ExactoError, Result};

/// Lattice basis for L = {f ∈ Z^d : Σ f_i · b^i ≡ 0 mod p}.
///
/// Adding any L-vector to a dBFV digit representation preserves the plaintext,
/// since L-vectors evaluate to 0 mod p at B=b.
#[derive(Clone, Debug)]
pub struct LatticeBasis {
    /// d×d basis matrix. Row i is the i-th basis vector.
    pub basis: Vec<Vec<BigInt>>,
    pub base: u64,
    pub d: usize,
    pub p: u64,
}

/// Precomputed small representatives of B^j mod L for j ≥ d.
#[derive(Clone, Debug)]
pub struct SmallReps {
    pub reps: Vec<Vec<i64>>,
    pub base: u64,
    pub d: usize,
}

/// Precomputed data for efficient lattice reduction.
#[derive(Clone, Debug)]
pub struct LatticeReducer {
    pub basis: LatticeBasis,
    pub small_reps: SmallReps,
    /// Gram-Schmidt orthogonalized basis vectors (as f64).
    pub gs_basis: Vec<Vec<f64>>,
    /// Gram-Schmidt squared norms ||b*_i||^2.
    pub gs_norms_sq: Vec<f64>,
    /// Gram-Schmidt coefficients μ_{i,j}.
    pub gs_mu: Vec<Vec<f64>>,
}

impl LatticeBasis {
    /// Construct the lattice basis.
    ///
    /// Basis rows (from paper §4.4):
    ///   r_i = b·e_i - e_{i+1}  for i=0,...,d-2
    ///   r_{d-1} = p·e_0
    ///
    /// These satisfy:
    ///   <r_i, (1,b,...,b^{d-1})> = b·b^i - b^{i+1} = 0 for i < d-1
    ///   <r_{d-1}, (1,b,...,b^{d-1})> = p ≡ 0 mod p
    ///
    /// Determinant = p (the lattice has index p in Z^d).
    pub fn new(base: u64, d: usize, p: u64) -> Result<Self> {
        if d < 1 {
            return Err(ExactoError::InvalidParam("d must be >= 1".into()));
        }
        if base < 2 {
            return Err(ExactoError::InvalidParam("base must be >= 2".into()));
        }

        let b = BigInt::from(base);
        let p_big = BigInt::from(p);
        let mut basis = Vec::with_capacity(d);

        for i in 0..d.saturating_sub(1) {
            let mut row = vec![BigInt::zero(); d];
            row[i] = b.clone();
            row[i + 1] = -BigInt::one();
            basis.push(row);
        }

        let mut last_row = vec![BigInt::zero(); d];
        last_row[0] = p_big;
        basis.push(last_row);

        Ok(Self { basis, base, d, p })
    }
}

impl SmallReps {
    /// Compute small representatives of B^j for j = d, ..., 2d-2.
    ///
    /// B^j evaluated at B=b gives b^j. We need (a_0,...,a_{d-1}) with small entries
    /// such that Σ a_i · b^i ≡ b^j mod p. This is the digit decomposition of b^j mod p,
    /// then reduced via the lattice.
    pub fn compute(base: u64, d: usize, p: u64, reducer: &LatticeReducer) -> Self {
        let mut reps = Vec::new();

        for j in d..=(2 * d - 2) {
            let val = mod_pow_u128(base as u128, j as u128, p as u128) as u64;
            let digits: Vec<i64> = crate::dbfv::decomposition::digit_decompose(val, base, d)
                .into_iter()
                .map(|x| x as i64)
                .collect();
            // Reduce these digits via the lattice to get small representatives
            let reduced = reduce_digits(&digits, reducer);
            reps.push(reduced);
        }

        Self { reps, base, d }
    }

    /// Simple version without reducer (uses raw digit decomposition).
    /// p=0 means p=2^64 (full u64 range).
    pub fn compute_simple(base: u64, d: usize, p: u64) -> Self {
        let mut reps = Vec::new();

        for j in d..=(2 * d - 2) {
            let val = if p == 0 {
                // p = 2^64: compute base^j mod 2^64 via wrapping
                (base as u64).wrapping_pow(j as u32)
            } else {
                mod_pow_u128(base as u128, j as u128, p as u128) as u64
            };
            let digits: Vec<i64> = crate::dbfv::decomposition::digit_decompose(val, base, d)
                .into_iter()
                .map(|x| x as i64)
                .collect();
            reps.push(digits);
        }

        Self { reps, base, d }
    }
}

impl LatticeReducer {
    /// Build a LatticeReducer from parameters.
    pub fn new(base: u64, d: usize, p: u64) -> Result<Self> {
        let basis = LatticeBasis::new(base, d, p)?;
        let (gs_basis, gs_mu, gs_norms_sq) = gram_schmidt_bigint(&basis);
        let small_reps = SmallReps::compute_simple(base, d, p);

        Ok(Self {
            basis,
            small_reps,
            gs_basis,
            gs_norms_sq,
            gs_mu,
        })
    }
}

/// Compute Gram-Schmidt orthogonalization of the lattice basis.
/// Uses BigInt internally for precision, converts to f64 for the final result.
fn gram_schmidt_bigint(basis: &LatticeBasis) -> (Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<f64>) {
    let d = basis.d;
    let mut b_star: Vec<Vec<f64>> = Vec::with_capacity(d);
    let mut mu: Vec<Vec<f64>> = vec![vec![0.0; d]; d];
    let mut norms_sq: Vec<f64> = Vec::with_capacity(d);

    for i in 0..d {
        let mut bi: Vec<f64> = basis.basis[i].iter()
            .map(|x| bigint_to_f64(x))
            .collect();

        for j in 0..i {
            let dot_ij = dot_f64(&bi, &b_star[j]);
            mu[i][j] = if norms_sq[j].abs() > 1e-15 { dot_ij / norms_sq[j] } else { 0.0 };

            for k in 0..d {
                bi[k] -= mu[i][j] * b_star[j][k];
            }
        }

        let norm_sq = dot_f64(&bi, &bi);
        norms_sq.push(norm_sq);
        b_star.push(bi);
    }

    (b_star, mu, norms_sq)
}

/// Hybrid lattice reduction on a vector of digit values.
///
/// Given v = (v_0, ..., v_{d-1}) with possibly large entries, find v' = v - w
/// where w ∈ L such that all |v'_i| ≈ b.
///
/// Uses Babai's nearest-plane algorithm with the precomputed Gram-Schmidt basis.
///
/// Paper §4.6.2: This is a PUBLIC operation — it operates on digit representations
/// directly, not on encrypted data.
pub fn reduce_digits(digits: &[i64], reducer: &LatticeReducer) -> Vec<i64> {
    let d = reducer.basis.d;
    assert_eq!(digits.len(), d);

    // Babai's nearest-plane algorithm:
    // Process basis vectors from last to first.
    // For each basis vector b_i, compute the coefficient c_i = round(<t, b*_i> / ||b*_i||^2)
    // and subtract c_i * b_i from the target.

    let mut t: Vec<f64> = digits.iter().map(|&x| x as f64).collect();

    // Compute integer coefficients of the closest lattice vector
    let mut coeffs = vec![0i64; d];

    for i in (0..d).rev() {
        let dot_tb = dot_f64(&t, &reducer.gs_basis[i]);
        let c = if reducer.gs_norms_sq[i].abs() > 1e-15 {
            (dot_tb / reducer.gs_norms_sq[i]).round() as i64
        } else {
            0
        };
        coeffs[i] = c;

        // Subtract c * basis[i] from target (using exact integer basis values)
        for k in 0..d {
            let bik = bigint_to_f64(&reducer.basis.basis[i][k]);
            t[k] -= c as f64 * bik;
        }
    }

    // The reduced vector is t (= original - lattice_vector)
    // Convert back to i64 with rounding
    t.iter().map(|&x| x.round() as i64).collect()
}

/// Reduce a vector of digit values, returning the result in [0, modulus) form.
/// This is the main entry point for dBFV lattice reduction.
pub fn reduce_digits_unsigned(digits: &[i64], reducer: &LatticeReducer, bfv_plain_mod: u64) -> Vec<u64> {
    let reduced = reduce_digits(digits, reducer);
    reduced.iter().map(|&x| {
        ((x as i128 % bfv_plain_mod as i128 + bfv_plain_mod as i128) % bfv_plain_mod as i128) as u64
    }).collect()
}

fn dot_f64(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn bigint_to_f64(x: &BigInt) -> f64 {
    // For our use case, values fit comfortably in f64
    x.to_string().parse::<f64>().unwrap_or(0.0)
}

fn mod_pow_u128(mut base: u128, mut exp: u128, modulus: u128) -> u128 {
    let mut result = 1u128;
    base %= modulus;
    while exp > 0 {
        if exp & 1 == 1 {
            result = result * base % modulus;
        }
        exp >>= 1;
        if exp > 0 {
            base = base * base % modulus;
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lattice_basis() {
        let lb = LatticeBasis::new(16, 2, 256).unwrap();
        assert_eq!(lb.basis.len(), 2);
        assert_eq!(lb.basis[0], vec![BigInt::from(16), BigInt::from(-1)]);
        assert_eq!(lb.basis[1], vec![BigInt::from(256), BigInt::from(0)]);
    }

    #[test]
    fn test_lattice_vectors_in_kernel() {
        let b = 16u64;
        let p = 256u64;
        let lb = LatticeBasis::new(b, 2, p).unwrap();

        for row in &lb.basis {
            let eval: i128 = row.iter().enumerate()
                .map(|(j, a)| {
                    let a_i64: i64 = a.to_string().parse().unwrap();
                    a_i64 as i128 * (b as i128).pow(j as u32)
                })
                .sum();
            assert_eq!(eval.rem_euclid(p as i128), 0,
                "basis vector {:?} evaluates to {} mod {}", row, eval % p as i128, p);
        }
    }

    #[test]
    fn test_reduce_digits_identity() {
        // Digits already in range should be mostly unchanged
        let reducer = LatticeReducer::new(16, 2, 256).unwrap();
        let digits = vec![3i64, 7];
        let reduced = reduce_digits(&digits, &reducer);
        assert_eq!(reduced, vec![3, 7]);
    }

    #[test]
    fn test_reduce_digits_overflow() {
        // After convolution: 3*7=21 at position 0. This needs reduction.
        // [21, 0] should reduce to [5, 1] since 5 + 1*16 = 21 ≡ 21 mod 256.
        let reducer = LatticeReducer::new(16, 2, 256).unwrap();
        let digits = vec![21i64, 0];
        let reduced = reduce_digits(&digits, &reducer);

        // Verify: reduced evaluates to same value mod p
        let original_eval: i64 = 21 + 0 * 16;
        let reduced_eval: i64 = reduced.iter().enumerate()
            .map(|(j, &a)| a * 16i64.pow(j as u32))
            .sum();
        assert_eq!(
            original_eval.rem_euclid(256),
            reduced_eval.rem_euclid(256),
            "reduction changed the value: original={original_eval}, reduced={reduced_eval}"
        );

        // Verify: all digits are small (in range roughly [-b/2, b/2])
        for &d in &reduced {
            assert!(d.abs() <= 16,
                "digit {} exceeds base 16", d);
        }
    }

    #[test]
    fn test_reduce_digits_subtraction() {
        // After sub: 50-20. Digits [2,3] - [4,1] = [-2, 2] (as signed integers)
        // But BFV gives [14, 2] (in Z_16). We need to handle this.
        // Actually, reduce_digits works on signed integers.
        // Input: [14, 2] → evaluate = 14 + 32 = 46.
        // Need to find [v0, v1] with small entries such that v0 + 16*v1 ≡ 46 mod 256.
        // 46 = -2 + 3*16 = [-2, 3], or 46 = 14 + 2*16.
        // The nearest-plane should give [-2, 3] or similar small representation.
        let reducer = LatticeReducer::new(16, 2, 256).unwrap();
        let digits = vec![14i64, 2];
        let reduced = reduce_digits(&digits, &reducer);

        let original_eval: i64 = digits.iter().enumerate()
            .map(|(j, &a)| a * 16i64.pow(j as u32))
            .sum();
        let reduced_eval: i64 = reduced.iter().enumerate()
            .map(|(j, &a)| a * 16i64.pow(j as u32))
            .sum();
        assert_eq!(
            original_eval.rem_euclid(256),
            reduced_eval.rem_euclid(256),
        );
    }

    #[test]
    fn test_reduce_preserves_plaintext() {
        // General test: for various vectors, reduction preserves evaluation mod p
        let reducer = LatticeReducer::new(16, 2, 256).unwrap();

        for v0 in [-30i64, -5, 0, 5, 21, 50, 100] {
            for v1 in [-10i64, 0, 5, 15, 30] {
                let digits = vec![v0, v1];
                let reduced = reduce_digits(&digits, &reducer);

                let original_eval: i64 = v0 + v1 * 16;
                let reduced_eval: i64 = reduced[0] + reduced[1] * 16;

                assert_eq!(
                    original_eval.rem_euclid(256),
                    reduced_eval.rem_euclid(256),
                    "reduction changed value for [{v0}, {v1}]: original_eval={original_eval}, reduced_eval={reduced_eval}"
                );
            }
        }
    }
}
