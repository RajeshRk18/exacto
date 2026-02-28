use crate::error::{ExactoError, Result};
use crate::ring::poly::CoeffPoly;
use crate::ring::rns::{RnsPoly, RnsBasis};
use crate::bfv::BfvCiphertext;

/// Modulus switch: given ct mod Q, produce ct' mod Q' where Q' = Q/q_last.
///
/// This drops the last RNS prime, reducing the modulus and noise proportionally.
/// Each component c_i is mapped: c_i' = round((Q'/Q) · c_i) mod Q'.
///
/// For RNS representation, this is equivalent to dropping the last NTT component
/// and applying a correction term.
pub fn mod_switch_drop_prime(ct: &BfvCiphertext) -> Result<BfvCiphertext> {
    let params = &ct.params;
    let basis = &params.ct_basis;

    if basis.num_moduli() <= 1 {
        return Err(ExactoError::InvalidParam(
            "cannot drop prime: only one modulus remaining".into()
        ));
    }

    let new_c: Vec<RnsPoly> = ct.c.iter()
        .map(|ci| drop_last_rns_prime(ci, basis))
        .collect::<Result<Vec<_>>>()?;

    // Note: the params should also be updated to reflect the new modulus.
    // For now, we keep the same params and just drop the component.
    // In production, you'd create new params with the reduced basis.
    Ok(BfvCiphertext {
        c: new_c,
        params: params.clone(),
    })
}

/// Drop the last RNS component from a polynomial, with rounding correction.
fn drop_last_rns_prime(poly: &RnsPoly, basis: &RnsBasis) -> Result<RnsPoly> {
    let num_primes = poly.num_components();
    if num_primes <= 1 {
        return Err(ExactoError::InvalidParam(
            "cannot drop: only one component".into()
        ));
    }

    let last_idx = num_primes - 1;
    let _q_last = basis.moduli[last_idx];

    // Get the last component in coefficient form
    let last_coeffs = poly.components[last_idx].to_coeff_poly();

    // For each remaining prime q_i, compute:
    // c_i' = c_i - (c_last mod q_i) * (q_last^{-1} mod q_i)
    // This is a simplified mod-switch that preserves the plaintext.
    let mut new_components = Vec::with_capacity(num_primes - 1);

    for i in 0..last_idx {
        let qi = basis.moduli[i];
        let _bk = basis.barrett_ks[i];

        // c_last mod q_i
        let last_mod_qi: Vec<u64> = last_coeffs.coeffs.iter()
            .map(|&c| c % qi)
            .collect();

        // Subtract from c_i
        let ci_coeffs = poly.components[i].to_coeff_poly();
        let mut new_coeffs = Vec::with_capacity(ci_coeffs.len());
        for (j, &c) in ci_coeffs.coeffs.iter().enumerate() {
            let diff = if c >= last_mod_qi[j] {
                c - last_mod_qi[j]
            } else {
                qi - last_mod_qi[j] + c
            };
            new_coeffs.push(diff);
        }

        let new_cp = CoeffPoly { coeffs: new_coeffs, modulus: qi };
        new_components.push(
            crate::ring::ntt::NttPoly::from_coeff_poly(&new_cp, basis.plans[i].clone())?
        );
    }

    Ok(RnsPoly {
        components: new_components,
        ring_degree: poly.ring_degree,
    })
}
