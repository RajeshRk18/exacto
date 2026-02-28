use crate::error::{ExactoError, Result};
use crate::ring::poly::CoeffPoly;
use crate::ring::rns::RnsPoly;
use crate::bfv::{BfvCiphertext, keygen::RelinKey};

/// Gadget decomposition of a polynomial using balanced digits.
///
/// For coefficient c (interpreted in centered form), compute digits d_0, d_1, ..., d_{k-1} such that
/// c = Σ d_i · base^i (mod q), where each d_i is centered in [-base/2, base/2).
/// Digits are returned modulo q (negative digits are represented as q - |d_i|).
pub fn gadget_decompose(
    poly: &CoeffPoly,
    base: u64,
    num_digits: usize,
) -> Vec<CoeffPoly> {
    let q = poly.modulus;
    let base_i = base as i128;
    let half_base = base_i / 2;
    let q_i = q as i128;
    let half_q = q / 2;

    let mut digits_coeffs = vec![vec![0u64; poly.len()]; num_digits];

    for (pos, &c_mod_q) in poly.coeffs.iter().enumerate() {
        // Center coefficient to signed representative.
        let mut remaining = if c_mod_q > half_q {
            c_mod_q as i128 - q_i
        } else {
            c_mod_q as i128
        };

        for d in 0..num_digits {
            let mut rem = remaining % base_i;
            if rem < -half_base {
                rem += base_i;
            } else if rem >= half_base {
                rem -= base_i;
            }

            let rem_mod_q = ((rem % q_i) + q_i) % q_i;
            digits_coeffs[d][pos] = rem_mod_q as u64;
            remaining = (remaining - rem) / base_i;
        }
    }

    let mut digits = Vec::with_capacity(num_digits);
    for coeffs in digits_coeffs {
        digits.push(CoeffPoly { coeffs, modulus: q });
    }

    digits
}

/// Relinearize a degree-2 ciphertext (c0, c1, c2) to degree-1 (c0', c1').
///
/// Uses the relinearization key to eliminate c2:
/// c0' = c0 + Σ_i decompose_i(c2) · rlk0_i
/// c1' = c1 + Σ_i decompose_i(c2) · rlk1_i
pub fn relinearize(
    ct: &BfvCiphertext,
    rlk: &RelinKey,
) -> Result<BfvCiphertext> {
    if ct.c.len() < 3 {
        return Ok(ct.clone()); // Already degree-1
    }
    if ct.c.len() > 3 {
        return Err(ExactoError::InvalidParam(
            "relinearization only supports degree-2 ciphertexts".into()
        ));
    }

    let params = &ct.params;
    let basis = &params.ct_basis;

    // Extract c2 in coefficient form for decomposition
    let c2_coeffs = ct.c[2].to_coeff_poly(basis);

    // Gadget decompose c2
    let digits = gadget_decompose(&c2_coeffs, params.gadget_base, params.gadget_digits);

    // Compute: c0' = c0 + Σ digit_i * rlk0_i
    //          c1' = c1 + Σ digit_i * rlk1_i
    let mut c0_new = ct.c[0].clone();
    let mut c1_new = ct.c[1].clone();

    for (i, digit) in digits.iter().enumerate() {
        if i >= rlk.keys.len() {
            break;
        }
        let digit_rns = RnsPoly::from_coeff_poly(digit, basis)?;
        let prod0 = digit_rns.mul(&rlk.keys[i].0)?;
        let prod1 = digit_rns.mul(&rlk.keys[i].1)?;
        c0_new = c0_new.add(&prod0)?;
        c1_new = c1_new.add(&prod1)?;
    }

    Ok(BfvCiphertext {
        c: vec![c0_new, c1_new],
        params: params.clone(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ring::modular::{barrett_constant, mod_mul, mod_add};

    #[test]
    fn test_gadget_decompose() {
        // Balanced decomposition of 42 in base 16, 2 digits:
        // 42 = (-6) + 3*16, so digits are [-6, 3] modulo q.
        let poly = CoeffPoly::from_coeffs(vec![42], 65537);
        let digits = gadget_decompose(&poly, 16, 2);
        assert_eq!(digits.len(), 2);
        assert_eq!(digits[0].coeffs[0], 65531); // -6 mod q
        assert_eq!(digits[1].coeffs[0], 3);
    }

    #[test]
    fn test_gadget_decompose_reconstruct() {
        let q = 65537u64;
        let base = 16u64;
        let num_digits = 4; // base^4 = 65536 ≈ q

        let poly = CoeffPoly::from_coeffs(vec![12345, 54321, 100, 0], q);
        let digits = gadget_decompose(&poly, base, num_digits);

        // Reconstruct: Σ digit_i * base^i should give back original
        let bk = barrett_constant(q);
        for pos in 0..poly.len() {
            let mut reconstructed = 0u64;
            let mut power = 1u64;
            for d in 0..num_digits {
                let contrib = mod_mul(digits[d].coeffs[pos], power, q, bk);
                reconstructed = mod_add(reconstructed, contrib, q);
                power = mod_mul(power, base, q, bk);
            }
            assert_eq!(reconstructed, poly.coeffs[pos],
                "mismatch at position {pos}: got {reconstructed}, expected {}", poly.coeffs[pos]);
        }
    }

    #[test]
    fn test_gadget_decompose_balanced_negative() {
        // q-1 corresponds to -1 in centered form.
        let q = 65537u64;
        let poly = CoeffPoly::from_coeffs(vec![q - 1], q);
        let digits = gadget_decompose(&poly, 16, 4);

        // Least-significant balanced digit should be -1 mod q.
        assert_eq!(digits[0].coeffs[0], q - 1);
    }
}
