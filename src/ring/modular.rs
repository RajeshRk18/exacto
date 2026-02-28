/// Barrett reduction: compute a mod m using precomputed Barrett constant.
///
/// Barrett constant k = floor(2^64 / m). For moduli ≤ 2^32, the single-word
/// Barrett trick is exact for inputs a < m^2. For larger moduli, we fall back
/// to u128 division (still single-instruction on modern x86-64/aarch64).
#[inline(always)]
pub fn barrett_reduce(a: u128, m: u64, barrett_k: u64) -> u64 {
    if m > (1u64 << 32) {
        // For moduli > 2^32, products can exceed 2^64 and the Barrett
        // approximation error grows beyond 1. Use exact u128 division.
        (a % m as u128) as u64
    } else {
        // q_hat = (a * k) >> 64
        let q_hat = ((a * barrett_k as u128) >> 64) as u64;
        let r = (a as u64).wrapping_sub(q_hat.wrapping_mul(m));
        // At most one conditional subtraction needed
        if r >= m { r.wrapping_sub(m) } else { r }
    }
}

/// Compute Barrett constant for modulus m: floor(2^64 / m)
#[inline]
pub fn barrett_constant(m: u64) -> u64 {
    assert!(m > 1, "modulus must be > 1");
    // 2^64 / m — we compute via 128-bit division
    ((1u128 << 64) / m as u128) as u64
}

/// Montgomery form: aR mod m where R = 2^64
/// Montgomery reduction: given T < m*R, compute T*R^{-1} mod m
///
/// Requires m odd and m_inv = -m^{-1} mod 2^64.
#[inline(always)]
pub fn montgomery_reduce(t: u128, m: u64, m_inv_neg: u64) -> u64 {
    let t_lo = t as u64;
    let k = t_lo.wrapping_mul(m_inv_neg);
    let km = k as u128 * m as u128;
    let r = ((t.wrapping_add(km)) >> 64) as u64;
    if r >= m { r.wrapping_sub(m) } else { r }
}

/// Compute -m^{-1} mod 2^64 for Montgomery reduction
pub fn montgomery_inv_neg(m: u64) -> u64 {
    assert!(m & 1 == 1, "Montgomery requires odd modulus");
    // Newton's method to find m^{-1} mod 2^64
    let mut inv = m; // m * m ≡ m^2, start with m as initial guess
    // Each iteration doubles the number of correct bits
    for _ in 0..6 {
        inv = inv.wrapping_mul(2u64.wrapping_sub(m.wrapping_mul(inv)));
    }
    // inv is now m^{-1} mod 2^64, negate it
    inv.wrapping_neg()
}

/// Modular addition: (a + b) mod m, assumes a, b < m
#[inline(always)]
pub fn mod_add(a: u64, b: u64, m: u64) -> u64 {
    let sum = a as u128 + b as u128;
    let r = sum as u64;
    if sum >= m as u128 { r.wrapping_sub(m) } else { r }
}

/// Modular subtraction: (a - b) mod m, assumes a, b < m
#[inline(always)]
pub fn mod_sub(a: u64, b: u64, m: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        m - b + a
    }
}

/// Modular negation: (-a) mod m, assumes a < m
#[inline(always)]
pub fn mod_neg(a: u64, m: u64) -> u64 {
    if a == 0 { 0 } else { m - a }
}

/// Modular multiplication: (a * b) mod m using Barrett reduction
#[inline(always)]
pub fn mod_mul(a: u64, b: u64, m: u64, barrett_k: u64) -> u64 {
    let product = a as u128 * b as u128;
    barrett_reduce(product, m, barrett_k)
}

/// Modular exponentiation: a^exp mod m
pub fn mod_pow(mut base: u64, mut exp: u64, m: u64) -> u64 {
    let bk = barrett_constant(m);
    let mut result = 1u64;
    base %= m;
    while exp > 0 {
        if exp & 1 == 1 {
            result = mod_mul(result, base, m, bk);
        }
        exp >>= 1;
        base = mod_mul(base, base, m, bk);
    }
    result
}

/// Find modular inverse: a^{-1} mod m using extended Euclidean algorithm
pub fn mod_inv(a: u64, m: u64) -> Option<u64> {
    let (mut old_r, mut r) = (a as i128, m as i128);
    let (mut old_s, mut s) = (1i128, 0i128);

    while r != 0 {
        let q = old_r / r;
        let tmp = r;
        r = old_r - q * r;
        old_r = tmp;
        let tmp = s;
        s = old_s - q * s;
        old_s = tmp;
    }

    if old_r != 1 {
        return None; // not coprime
    }

    Some(((old_s % m as i128 + m as i128) % m as i128) as u64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_barrett_reduce() {
        let m = 65537u64;
        let bk = barrett_constant(m);
        assert_eq!(barrett_reduce(0, m, bk), 0);
        assert_eq!(barrett_reduce(1, m, bk), 1);
        assert_eq!(barrett_reduce(m as u128, m, bk), 0);
        assert_eq!(barrett_reduce(m as u128 + 1, m, bk), 1);
        assert_eq!(barrett_reduce(123456789u128, m, bk), (123456789u128 % m as u128) as u64);
    }

    #[test]
    fn test_mod_mul() {
        let m = 65537u64;
        let bk = barrett_constant(m);
        assert_eq!(mod_mul(1234, 5678, m, bk), ((1234u128 * 5678) % m as u128) as u64);
        assert_eq!(mod_mul(0, 5678, m, bk), 0);
        assert_eq!(mod_mul(1, 5678, m, bk), 5678);
    }

    #[test]
    fn test_mod_add_sub() {
        let m = 65537u64;
        assert_eq!(mod_add(100, 200, m), 300);
        assert_eq!(mod_add(m - 1, 2, m), 1);
        assert_eq!(mod_sub(200, 100, m), 100);
        assert_eq!(mod_sub(100, 200, m), m - 100);
    }

    #[test]
    fn test_mod_neg() {
        let m = 65537u64;
        assert_eq!(mod_neg(0, m), 0);
        assert_eq!(mod_neg(1, m), m - 1);
        assert_eq!(mod_add(100, mod_neg(100, m), m), 0);
    }

    #[test]
    fn test_mod_pow() {
        let m = 65537u64;
        assert_eq!(mod_pow(2, 10, m), 1024);
        assert_eq!(mod_pow(2, 16, m), 65536 % m);
        assert_eq!(mod_pow(3, 0, m), 1);
    }

    #[test]
    fn test_mod_inv() {
        let m = 65537u64;
        let bk = barrett_constant(m);
        let a = 12345u64;
        let inv = mod_inv(a, m).unwrap();
        assert_eq!(mod_mul(a, inv, m, bk), 1);
    }

    #[test]
    fn test_montgomery() {
        let m = 65537u64; // odd modulus
        let m_inv_neg = montgomery_inv_neg(m);
        // Verify: m * (-m^{-1}) ≡ -1 (mod 2^64), i.e., m * m_inv_neg ≡ 2^64 - 1... no:
        // m_inv_neg = -m^{-1} mod 2^64, so m * m_inv_neg ≡ -1 mod 2^64
        // which means m * m_inv_neg + 1 ≡ 0 mod 2^64
        assert_eq!(m.wrapping_mul(m_inv_neg).wrapping_add(1), 0);

        // Test Montgomery reduce: T = a*R where R = 2^64
        // montgomery_reduce(a*R, m, m_inv_neg) should give a mod m
        let a = 12345u64;
        // a * R mod ... we need T < m*R
        // Let's just test the round-trip: to_mont(a) = a*R mod m, from_mont(t) = t*R^{-1} mod m
        let _a_r = (a as u128 * (1u128 << 64)) % m as u128; // a*R mod m
        // For a simple test, let's verify reduce of (a * R^2 mod m) gives a*R mod m
        // Actually, let's just check that montgomery_reduce recovers from a known product
        let t = a as u128 * m as u128; // t = a*m, so t*R^{-1} mod m = 0
        // Wait, t = a*m ≡ 0 mod m, so montgomery_reduce should give 0
        let r = montgomery_reduce(t, m, m_inv_neg);
        assert_eq!(r, 0);
    }
}
