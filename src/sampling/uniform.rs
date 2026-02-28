use rand::Rng;
use crate::ring::poly::CoeffPoly;

/// Sample a polynomial with uniformly random coefficients in [0, modulus).
pub fn sample_uniform_poly<R: Rng>(n: usize, modulus: u64, rng: &mut R) -> CoeffPoly {
    // Rejection sampling to avoid bias
    let mask = if modulus.is_power_of_two() {
        modulus - 1
    } else {
        (1u64 << (64 - modulus.leading_zeros())) - 1
    };

    let coeffs: Vec<u64> = (0..n)
        .map(|_| {
            loop {
                let val = rng.random::<u64>() & mask;
                if val < modulus {
                    break val;
                }
            }
        })
        .collect();
    CoeffPoly { coeffs, modulus }
}

/// Sample a polynomial with ternary coefficients {-1, 0, 1}.
/// Each coefficient is independently: -1 with prob 1/3, 0 with prob 1/3, 1 with prob 1/3.
/// Stored as {q-1, 0, 1} mod q.
pub fn sample_ternary_poly<R: Rng>(n: usize, modulus: u64, rng: &mut R) -> CoeffPoly {
    let coeffs: Vec<u64> = (0..n)
        .map(|_| {
            // Use rejection sampling on 2 bits for uniform {0,1,2}
            let val = loop {
                let r = rng.random::<u8>() & 0x03;
                if r < 3 { break r; }
            };
            match val {
                0 => modulus - 1, // -1 mod q
                1 => 0,
                2 => 1,
                _ => unreachable!(),
            }
        })
        .collect();
    CoeffPoly { coeffs, modulus }
}

/// Sample a polynomial with binary coefficients {0, 1}.
pub fn sample_binary_poly<R: Rng>(n: usize, modulus: u64, rng: &mut R) -> CoeffPoly {
    let coeffs: Vec<u64> = (0..n)
        .map(|_| rng.random::<u64>() & 1)
        .collect();
    CoeffPoly { coeffs, modulus }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_uniform() {
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let poly = sample_uniform_poly(1024, 65537, &mut rng);
        assert_eq!(poly.len(), 1024);
        for &c in &poly.coeffs {
            assert!(c < 65537);
        }
    }

    #[test]
    fn test_ternary() {
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let q = 65537u64;
        let poly = sample_ternary_poly(1024, q, &mut rng);
        for &c in &poly.coeffs {
            assert!(c == 0 || c == 1 || c == q - 1, "invalid ternary coeff: {c}");
        }
        // Check roughly 1/3 each
        let zeros = poly.coeffs.iter().filter(|&&c| c == 0).count();
        let ones = poly.coeffs.iter().filter(|&&c| c == 1).count();
        let neg_ones = poly.coeffs.iter().filter(|&&c| c == q - 1).count();
        assert!(zeros > 200 && zeros < 500, "zeros = {zeros}");
        assert!(ones > 200 && ones < 500, "ones = {ones}");
        assert!(neg_ones > 200 && neg_ones < 500, "neg_ones = {neg_ones}");
    }

    #[test]
    fn test_binary() {
        let mut rng = ChaCha20Rng::seed_from_u64(42);
        let poly = sample_binary_poly(1024, 65537, &mut rng);
        for &c in &poly.coeffs {
            assert!(c == 0 || c == 1);
        }
    }
}
