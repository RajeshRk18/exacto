/// LWE security estimation.
///
/// Estimates the security level (in bits) for given LWE parameters using
/// the lattice estimator heuristic from Albrecht et al.
///
/// This is a simplified estimator. For production use, cross-validate with
/// the lattice-estimator sage tool.

/// Estimate security level for RLWE parameters.
///
/// Uses the core-SVP methodology:
/// security ≈ n · log2(q/σ) · (1/δ) where δ ≈ 1.005 for BKZ-300
pub fn estimate_security_bits(ring_degree: usize, log2_q: f64, sigma: f64) -> f64 {
    // Hermite factor for BKZ-300 (conservative estimate)
    let log2_delta = 0.00722; // log2(1.005)

    // Root Hermite factor model:
    // BKZ dimension β needed to achieve Hermite factor δ = δ_0^{2n}
    // Security ≈ 0.292 * β (Core-SVP model)

    let n = ring_degree as f64;
    let log2_sigma = sigma.log2();

    // For RLWE with ring dimension n and modulus q:
    // The attacker needs to find a vector of norm ≈ σ·√n in a lattice of volume q^n
    // Using the GSA model: δ^{2nβ} · (q^n)^{1/β} ≈ σ·√n

    // Simplified Albrecht model for dimension-n RLWE:
    // security ≈ n · log2(q / (σ·√n)) / (2·log2(δ))
    // This is then converted to bit security via the core-SVP model

    let log2_norm = log2_sigma + 0.5 * n.log2();
    let log2_advantage = log2_q - log2_norm / n;

    // BKZ block size estimate
    let beta = log2_advantage / (2.0 * log2_delta);

    // Core-SVP bit security
    0.292 * beta * n / ring_degree as f64
}

/// Check if parameters meet a minimum security level.
pub fn check_security(ring_degree: usize, log2_q: f64, sigma: f64, min_bits: f64) -> bool {
    estimate_security_bits(ring_degree, log2_q, sigma) >= min_bits
}

/// Find minimum ring degree for given security level and modulus size.
pub fn min_ring_degree(log2_q: f64, sigma: f64, target_bits: f64) -> usize {
    let mut n = 1024;
    while n <= 65536 {
        if estimate_security_bits(n, log2_q, sigma) >= target_bits {
            return n;
        }
        n *= 2;
    }
    n
}
