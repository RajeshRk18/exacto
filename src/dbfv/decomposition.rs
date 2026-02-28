
/// Decompose a value in Z_p into d digits in base b.
///
/// Returns digits [d_0, d_1, ..., d_{d-1}] such that
/// value = Σ d_i · b^i (mod p).
///
/// Each digit d_i ∈ [0, b).
pub fn digit_decompose(value: u64, base: u64, num_digits: usize) -> Vec<u64> {
    let mut digits = Vec::with_capacity(num_digits);
    let mut remaining = value;
    for _ in 0..num_digits {
        digits.push(remaining % base);
        remaining /= base;
    }
    digits
}

/// Recompose digits in base b to a value in Z_p.
///
/// value = Σ d_i · b^i (mod p).
/// modulus=0 means p=2^64 (full u64 range, reduction via truncation).
pub fn digit_recompose(digits: &[u64], base: u64, modulus: u64) -> u64 {
    let mut result = 0u128;
    let mut power = 1u128;
    for &d in digits {
        result += d as u128 * power;
        power *= base as u128;
    }
    if modulus == 0 {
        result as u64 // truncation = mod 2^64
    } else {
        (result % modulus as u128) as u64
    }
}

/// Recompose digits using centered (signed) representation.
///
/// Each digit is in Z_t (BFV plaintext modulus). Values > t/2 are interpreted
/// as negative (digit - t). This is needed because homomorphic operations
/// (especially subtraction) can produce "negative" digit values that wrap
/// around mod t.
///
/// value = Σ centered(d_i) · b^i (mod p)
/// modulus=0 means p=2^64 (full u64 range, reduction via truncation).
pub fn digit_recompose_signed(
    digits: &[u64],
    base: u64,
    modulus: u64,
    bfv_plain_mod: u64,
) -> u64 {
    let half_t = bfv_plain_mod / 2;
    let mut result = 0i128;
    let mut power = 1i128;
    for &d in digits {
        let centered = if d > half_t {
            d as i128 - bfv_plain_mod as i128
        } else {
            d as i128
        };
        result += centered * power;
        power *= base as i128;
    }
    if modulus == 0 {
        result as u64 // truncation = mod 2^64
    } else {
        ((result % modulus as i128 + modulus as i128) % modulus as i128) as u64
    }
}

/// Decompose a polynomial coefficient-wise: each coefficient of the Z_p polynomial
/// is decomposed into d digits in base b.
///
/// Returns d polynomials, where the i-th polynomial contains the i-th digit of
/// each coefficient.
pub fn poly_digit_decompose(
    coeffs: &[u64],
    base: u64,
    num_digits: usize,
) -> Vec<Vec<u64>> {
    let mut digit_polys = vec![Vec::with_capacity(coeffs.len()); num_digits];

    for &c in coeffs {
        let digits = digit_decompose(c, base, num_digits);
        for (i, &d) in digits.iter().enumerate() {
            digit_polys[i].push(d);
        }
    }

    digit_polys
}

/// Recompose d digit polynomials back into a single polynomial.
pub fn poly_digit_recompose(
    digit_polys: &[Vec<u64>],
    base: u64,
    modulus: u64,
) -> Vec<u64> {
    if digit_polys.is_empty() {
        return Vec::new();
    }
    let n = digit_polys[0].len();
    (0..n).map(|i| {
        let digits: Vec<u64> = digit_polys.iter().map(|dp| dp[i]).collect();
        digit_recompose(&digits, base, modulus)
    }).collect()
}

/// Recompose digit polynomials using centered (signed) digit interpretation.
///
/// Each digit entry is interpreted in Z_t where values > t/2 are treated as
/// negative. This mirrors [`digit_recompose_signed`] coefficient-wise.
pub fn poly_digit_recompose_signed(
    digit_polys: &[Vec<u64>],
    base: u64,
    modulus: u64,
    bfv_plain_mod: u64,
) -> Vec<u64> {
    if digit_polys.is_empty() {
        return Vec::new();
    }

    let n = digit_polys.iter().map(|dp| dp.len()).min().unwrap_or(0);
    (0..n).map(|i| {
        let digits: Vec<u64> = digit_polys.iter().map(|dp| dp[i]).collect();
        digit_recompose_signed(&digits, base, modulus, bfv_plain_mod)
    }).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decompose_recompose() {
        let base = 16u64;
        let num_digits = 4;
        let p = 65536u64; // 16^4

        for value in [0, 1, 42, 255, 12345, 65535] {
            let digits = digit_decompose(value, base, num_digits);
            let recovered = digit_recompose(&digits, base, p);
            assert_eq!(recovered, value, "failed for value={value}");
        }
    }

    #[test]
    fn test_decompose_digits() {
        // 42 = 10 + 2*16 in base 16
        let digits = digit_decompose(42, 16, 2);
        assert_eq!(digits, vec![10, 2]);

        // 255 = 15 + 15*16 in base 16
        let digits = digit_decompose(255, 16, 2);
        assert_eq!(digits, vec![15, 15]);
    }

    #[test]
    fn test_poly_decompose_recompose() {
        let base = 16u64;
        let d = 2;
        let p = 256u64;

        let coeffs = vec![42, 100, 255, 0];
        let digit_polys = poly_digit_decompose(&coeffs, base, d);
        let recovered = poly_digit_recompose(&digit_polys, base, p);
        assert_eq!(recovered, coeffs);
    }

    #[test]
    fn test_poly_decompose_recompose_signed() {
        // Coeff0: [-2, 2] in base 16 -> 30 mod 256
        // Coeff1: [-1, 0] in base 16 -> 255 mod 256
        // Digits encoded in Z_t with t=929: -2 -> 927, -1 -> 928.
        let digit_polys = vec![
            vec![927, 928],
            vec![2, 0],
        ];
        let recovered = poly_digit_recompose_signed(&digit_polys, 16, 256, 929);
        assert_eq!(recovered, vec![30, 255]);
    }
}
