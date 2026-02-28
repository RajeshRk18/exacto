use rand::Rng;
use crate::ring::poly::CoeffPoly;

/// Sample a polynomial with coefficients from a discrete Gaussian distribution
/// using the FACCT (Fast Constant-time Approximate CDT) method.
///
/// Parameters:
/// - `n`: ring degree
/// - `modulus`: coefficient modulus
/// - `sigma`: standard deviation (typically 3.2 for FHE)
/// - `rng`: random number generator
///
/// The FACCT method provides constant-time sampling to avoid timing side channels.
/// For σ = 3.2, the tail bound is effectively zero beyond ±6σ ≈ ±19.
pub fn sample_gaussian_poly<R: Rng>(n: usize, modulus: u64, sigma: f64, rng: &mut R) -> CoeffPoly {
    let coeffs: Vec<u64> = (0..n)
        .map(|_| {
            let sample = sample_discrete_gaussian(sigma, rng);
            if sample >= 0 {
                sample as u64 % modulus
            } else {
                (modulus as i64 + (sample % modulus as i64)) as u64 % modulus
            }
        })
        .collect();
    CoeffPoly { coeffs, modulus }
}

/// Sample a single value from the discrete Gaussian distribution over Z
/// with standard deviation σ, centered at 0.
///
/// Uses a constant-time CDT (cumulative distribution table) approach for small σ.
/// For σ = 3.2, we only need to consider values in [-19, 19].
///
/// The scan is branchless: every CDF entry is visited and a conditional select
/// (using bitwise ops on integer masks) determines the result, avoiding
/// data-dependent branches that could leak timing information.
fn sample_discrete_gaussian<R: Rng>(sigma: f64, rng: &mut R) -> i64 {
    let tail = (6.0 * sigma).ceil() as i64;

    // CDT: precompute cumulative probabilities (unnormalized)
    let table_size = (2 * tail + 1) as usize;
    let mut cdf = Vec::with_capacity(table_size);
    let mut cumulative = 0.0f64;
    let two_sigma_sq = 2.0 * sigma * sigma;

    for x in -tail..=tail {
        let prob = (-((x * x) as f64) / two_sigma_sq).exp();
        cumulative += prob;
        cdf.push(cumulative);
    }

    let total = cumulative;
    let u: f64 = rng.random::<f64>() * total;

    // Branchless constant-time scan: iterate in reverse, always selecting
    // the lowest index i where u < cdf[i]. Uses integer masks to avoid
    // data-dependent branches.
    let mut result = tail; // default: maximum value
    for i in (0..table_size).rev() {
        // Convert comparison to an all-ones or all-zeros mask
        let cmp = u < cdf[i];
        let mask = (cmp as i64).wrapping_neg(); // 0 → 0, 1 → -1 (all bits set)
        let candidate = -tail + i as i64;
        // Conditional select: result = cmp ? candidate : result
        result = (candidate & mask) | (result & !mask);
    }

    result
}

/// Sample from a rounded Gaussian using Box-Muller (NOT constant time).
/// Use only for testing or non-security-critical paths.
pub fn sample_gaussian_box_muller<R: Rng>(sigma: f64, rng: &mut R) -> i64 {
    let u1: f64 = rng.random::<f64>();
    let u2: f64 = rng.random::<f64>();
    let z = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
    (z * sigma).round() as i64
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_gaussian_distribution() {
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let sigma = 3.2;
        let n = 10000;

        let samples: Vec<i64> = (0..n)
            .map(|_| sample_discrete_gaussian(sigma, &mut rng))
            .collect();

        // Check mean is close to 0
        let mean: f64 = samples.iter().map(|&x| x as f64).sum::<f64>() / n as f64;
        assert!(mean.abs() < 0.5, "mean = {mean}");

        // Check variance is close to σ²
        let var: f64 = samples.iter().map(|&x| (x as f64 - mean).powi(2)).sum::<f64>() / n as f64;
        let expected_var = sigma * sigma;
        assert!((var - expected_var).abs() < 2.0, "var = {var}, expected ≈ {expected_var}");

        // Check all samples are within [-6σ, 6σ]
        let tail = (6.0 * sigma).ceil() as i64;
        for &s in &samples {
            assert!(s.abs() <= tail, "sample {s} exceeds tail bound {tail}");
        }
    }

    #[test]
    fn test_gaussian_poly() {
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let poly = sample_gaussian_poly(1024, 65537, 3.2, &mut rng);
        assert_eq!(poly.len(), 1024);
        assert_eq!(poly.modulus, 65537);
        // All coefficients should be in [0, q)
        for &c in &poly.coeffs {
            assert!(c < 65537);
        }
    }
}
