use crate::error::{ExactoError, Result};
use num_bigint::{BigInt, BigUint};
use num_traits::{One, Signed, ToPrimitive, Zero};
use crate::ring::modular::{mod_mul, mod_inv, barrett_constant};
use crate::ring::ntt::NttPoly;
use crate::ring::poly::CoeffPoly;
use crate::ring::rns::{RnsBasis, RnsPoly};
use crate::bfv::{BfvCiphertext, keygen::{RelinKey, GaloisKey}};
use crate::bfv::encrypt::scale_plaintext;
use crate::bfv::keyswitch::{relinearize, gadget_decompose};

/// Homomorphic addition: ct_out = ct1 + ct2.
/// Component-wise addition of ciphertext polynomials.
pub fn bfv_add(ct1: &BfvCiphertext, ct2: &BfvCiphertext) -> Result<BfvCiphertext> {
    let max_len = ct1.c.len().max(ct2.c.len());
    let mut c = Vec::with_capacity(max_len);

    for i in 0..max_len {
        match (ct1.c.get(i), ct2.c.get(i)) {
            (Some(a), Some(b)) => c.push(a.add(b)?),
            (Some(a), None) => c.push(a.clone()),
            (None, Some(b)) => c.push(b.clone()),
            (None, None) => unreachable!(),
        }
    }

    Ok(BfvCiphertext {
        c,
        params: ct1.params.clone(),
    })
}

/// Homomorphic subtraction: ct_out = ct1 - ct2.
pub fn bfv_sub(ct1: &BfvCiphertext, ct2: &BfvCiphertext) -> Result<BfvCiphertext> {
    let max_len = ct1.c.len().max(ct2.c.len());
    let mut c = Vec::with_capacity(max_len);

    for i in 0..max_len {
        match (ct1.c.get(i), ct2.c.get(i)) {
            (Some(a), Some(b)) => c.push(a.sub(b)?),
            (Some(a), None) => c.push(a.clone()),
            (None, Some(b)) => c.push(b.neg()),
            (None, None) => unreachable!(),
        }
    }

    Ok(BfvCiphertext {
        c,
        params: ct1.params.clone(),
    })
}

/// Negate a ciphertext.
pub fn bfv_neg(ct: &BfvCiphertext) -> BfvCiphertext {
    let c = ct.c.iter().map(|ci| ci.neg()).collect();
    BfvCiphertext {
        c,
        params: ct.params.clone(),
    }
}

/// Homomorphic multiplication followed by relinearization.
///
/// Uses the simple tensor product approach:
/// For ct1 = (c0, c1), ct2 = (d0, d1):
///   - c_out_0 = round(p/q · c0·d0)
///   - c_out_1 = round(p/q · (c0·d1 + c1·d0))
///   - c_out_2 = round(p/q · c1·d1)
/// Then relinearize c_out_2 using the relin key.
///
/// For the RNS variant (HPS), we compute the tensor product in the Q basis,
/// then scale by p/q approximately using the RNS representation.
pub fn bfv_mul_and_relin(
    ct1: &BfvCiphertext,
    ct2: &BfvCiphertext,
    rlk: &RelinKey,
) -> Result<BfvCiphertext> {
    // First multiply to get a degree-2 ciphertext
    let ct_mul = bfv_mul_no_relin(ct1, ct2)?;
    // Then relinearize
    relinearize(&ct_mul, rlk)
}

/// Homomorphic multiplication WITHOUT relinearization.
/// Returns a degree-2 ciphertext (3 components).
///
/// Uses HPS RNS multiplication when an auxiliary basis is available (O(n log n)),
/// falls back to schoolbook multiplication (O(n^2)) otherwise.
pub fn bfv_mul_no_relin(
    ct1: &BfvCiphertext,
    ct2: &BfvCiphertext,
) -> Result<BfvCiphertext> {
    if ct1.c.len() != 2 || ct2.c.len() != 2 {
        return Err(ExactoError::InvalidParam(
            "multiplication requires degree-1 ciphertexts".into()
        ));
    }

    let params = &ct1.params;
    if params.ct_basis.num_moduli() > 1 {
        return bfv_mul_generic_rns(ct1, ct2);
    }
    if params.aux_basis.is_some() {
        bfv_mul_hps(ct1, ct2)
    } else {
        bfv_mul_schoolbook(ct1, ct2)
    }
}

/// Multi-prime BFV multiplication via exact CRT reconstruction and BigInt scaling.
///
/// This path is used when ciphertext modulus Q has multiple RNS primes.
fn bfv_mul_generic_rns(
    ct1: &BfvCiphertext,
    ct2: &BfvCiphertext,
) -> Result<BfvCiphertext> {
    let params = &ct1.params;
    let basis = &params.ct_basis;
    let n = params.ring_degree;
    let p = params.plain_modulus;

    let q_big_u = modulus_product_biguint(&basis.moduli);
    let q_big = BigInt::from(q_big_u.clone());
    let half_q = &q_big >> 1;

    let c0 = reconstruct_centered_bigint(&ct1.c[0], basis, &q_big_u)?;
    let c1 = reconstruct_centered_bigint(&ct1.c[1], basis, &q_big_u)?;
    let d0 = reconstruct_centered_bigint(&ct2.c[0], basis, &q_big_u)?;
    let d1 = reconstruct_centered_bigint(&ct2.c[1], basis, &q_big_u)?;

    let t0 = poly_mul_bigint(&c0, &d0, n);
    let t1 = poly_add_bigint(&poly_mul_bigint(&c0, &d1, n), &poly_mul_bigint(&c1, &d0, n));
    let t2 = poly_mul_bigint(&c1, &d1, n);

    let r0 = scale_tensor_component_bigint(&t0, p, &q_big, &half_q);
    let r1 = scale_tensor_component_bigint(&t1, p, &q_big, &half_q);
    let r2 = scale_tensor_component_bigint(&t2, p, &q_big, &half_q);

    let r0_rns = centered_bigint_to_rns(&r0, basis)?;
    let r1_rns = centered_bigint_to_rns(&r1, basis)?;
    let r2_rns = centered_bigint_to_rns(&r2, basis)?;

    Ok(BfvCiphertext {
        c: vec![r0_rns, r1_rns, r2_rns],
        params: params.clone(),
    })
}

/// HPS RNS multiplication: O(n log n) using auxiliary basis P.
///
/// Algorithm (Halevi-Polyakov-Shoup):
/// 1. Base-extend ct1, ct2 components from Q to P (centered representation).
/// 2. Compute tensor products in both Q and P bases (NTT multiply, O(n log n)).
/// 3. For each tensor coefficient, use CRT to recover the true value's scale:
///    Given a = t mod q, b = t mod P, compute m = (b-a)·q^{-1} mod P (centered).
///    Then round(p·t/q) = round(p·a_centered/q) + p·m.
fn bfv_mul_hps(
    ct1: &BfvCiphertext,
    ct2: &BfvCiphertext,
) -> Result<BfvCiphertext> {
    let params = &ct1.params;
    let ct_basis = &params.ct_basis;
    let aux_basis = params.aux_basis.as_ref()
        .ok_or_else(|| ExactoError::InvalidParam("HPS requires aux_basis".into()))?;
    let p = params.plain_modulus;
    let q = ct_basis.moduli[0];
    let big_p = aux_basis.moduli[0];

    // With a single aux prime, HPS centering requires P > n*Q/2.
    if aux_basis.moduli.len() == 1 {
        let min_required = (params.ring_degree as u128 * q as u128) / 2;
        if big_p as u128 <= min_required {
            return Err(ExactoError::InvalidParam(format!(
                "single aux prime too small for HPS centering: P={} <= n*Q/2={}",
                big_p, min_required
            )));
        }
    }

    // Step 1: Base-extend ct components from Q to P (centered representation)
    let c0_p = base_extend_centered(&ct1.c[0], q, aux_basis)?;
    let c1_p = base_extend_centered(&ct1.c[1], q, aux_basis)?;
    let d0_p = base_extend_centered(&ct2.c[0], q, aux_basis)?;
    let d1_p = base_extend_centered(&ct2.c[1], q, aux_basis)?;

    // Step 2: Compute tensor products in Q basis (NTT multiply)
    let t0_q = ct1.c[0].mul(&ct2.c[0])?;
    let t1_q_a = ct1.c[0].mul(&ct2.c[1])?;
    let t1_q_b = ct1.c[1].mul(&ct2.c[0])?;
    let t1_q = t1_q_a.add(&t1_q_b)?;
    let t2_q = ct1.c[1].mul(&ct2.c[1])?;

    // Step 3: Compute tensor products in P basis (NTT multiply)
    let t0_p = c0_p.mul(&d0_p)?;
    let t1_p_a = c0_p.mul(&d1_p)?;
    let t1_p_b = c1_p.mul(&d0_p)?;
    let t1_p = t1_p_a.add(&t1_p_b)?;
    let t2_p = c1_p.mul(&d1_p)?;

    // Step 4: Scale each tensor component using HPS formula
    let r0 = hps_scale(&t0_q, &t0_p, p, q, big_p, ct_basis, aux_basis)?;
    let r1 = hps_scale(&t1_q, &t1_p, p, q, big_p, ct_basis, aux_basis)?;
    let r2 = hps_scale(&t2_q, &t2_p, p, q, big_p, ct_basis, aux_basis)?;

    Ok(BfvCiphertext {
        c: vec![r0, r1, r2],
        params: params.clone(),
    })
}

/// Base-extend a polynomial from Q (single modulus q) to P (one or more primes)
/// using centered representation.
///
/// For each coefficient c ∈ [0, q):
/// - If c ≤ q/2: the "true" value is c, so c mod p_j = c % p_j.
/// - If c > q/2: the "true" value is c - q (negative), so (c-q) mod p_j = p_j - (q-c) % p_j.
fn base_extend_centered(
    poly: &RnsPoly,
    q: u64,
    basis_p: &RnsBasis,
) -> Result<RnsPoly> {
    let n = poly.ring_degree;
    let coeffs_q = poly.components[0].to_coeff_poly();
    let half_q = q / 2;

    let mut p_components = Vec::with_capacity(basis_p.moduli.len());
    for (j, (&pj, plan)) in basis_p.moduli.iter().zip(basis_p.plans.iter()).enumerate() {
        let _ = j;
        let mut coeffs_p = vec![0u64; n];
        for i in 0..n {
            let c = coeffs_q.coeffs[i];
            if c > half_q {
                // Centered value: c - q (negative). Compute (c - q) mod p_j = p_j - ((q - c) % p_j)
                let neg_val = q - c; // positive
                let rem = neg_val % pj;
                coeffs_p[i] = if rem == 0 { 0 } else { pj - rem };
            } else {
                coeffs_p[i] = c % pj;
            }
        }

        let cp = CoeffPoly { coeffs: coeffs_p, modulus: pj };
        p_components.push(NttPoly::from_coeff_poly(&cp, plan.clone())?);
    }

    Ok(RnsPoly { components: p_components, ring_degree: n })
}

/// HPS scaling: given tensor product mod Q and mod P, compute round(p·t/q) mod q.
///
/// Supports single-prime Q and one or two aux primes P.
/// For each coefficient i:
///   a = t_i mod q  (from Q-basis)
///   b_j = t_i mod p_j  (from P-basis)
///   Recover m = (t_i - a) / q via CRT on P residues
///   result = round(p · center(a,q) / q) + p · m  (mod q)
fn hps_scale(
    t_q: &RnsPoly,
    t_p: &RnsPoly,
    p: u64,      // BFV plaintext modulus
    q: u64,      // single ciphertext prime
    _big_p: u64, // unused (kept for signature compat)
    ct_basis: &RnsBasis,
    aux_basis: &RnsBasis,
) -> Result<RnsPoly> {
    let n = t_q.ring_degree;
    let num_aux = aux_basis.moduli.len();

    // Get tensor product mod q (single component INTT)
    let a_poly = t_q.components[0].to_coeff_poly();

    // Get tensor product mod each aux prime (INTT each component)
    let b_polys: Vec<CoeffPoly> = t_p.components.iter()
        .map(|c| c.to_coeff_poly())
        .collect();

    // Precompute q^{-1} mod p_j for each aux prime
    let mut q_inv_pj = Vec::with_capacity(num_aux);
    let mut bk_pj = Vec::with_capacity(num_aux);
    for j in 0..num_aux {
        let pj = aux_basis.moduli[j];
        q_inv_pj.push(
            mod_inv(q % pj, pj)
                .ok_or_else(|| ExactoError::InvalidParam("q not invertible mod p_j".into()))?
        );
        bk_pj.push(barrett_constant(pj));
    }

    let half_q = q / 2;
    let p_128 = p as i128;
    let q_128 = q as i128;

    let mut result = vec![0u64; n];

    if num_aux == 1 {
        // Single aux prime: direct formula (original algorithm)
        let big_p = aux_basis.moduli[0];
        let bkp = bk_pj[0];
        let half_p = big_p / 2;

        for i in 0..n {
            let a = a_poly.coeffs[i];
            let b = b_polys[0].coeffs[i];
            let a_centered: i128 = if a > half_q { a as i128 - q_128 } else { a as i128 };

            // Extend a_centered to P basis (must match centering used in round term)
            let a_ext = if a > half_q {
                let neg = q - a;
                let rem = neg % big_p;
                if rem == 0 { 0 } else { big_p - rem }
            } else {
                a % big_p
            };
            let diff = if b >= a_ext { b - a_ext } else { big_p - a_ext + b };
            let m_raw = mod_mul(diff, q_inv_pj[0], big_p, bkp);

            let m_centered: i128 = if m_raw > half_p {
                m_raw as i128 - big_p as i128
            } else {
                m_raw as i128
            };

            let pa = p_128 * a_centered;
            let round_pa_q = if pa >= 0 {
                (pa + q_128 / 2) / q_128
            } else {
                -((-pa + q_128 / 2) / q_128)
            };

            let scaled = round_pa_q + p_128 * m_centered;
            result[i] = ((scaled % q_128 + q_128) % q_128) as u64;
        }
    } else if num_aux == 2 {
        // Two aux primes: reconstruct m in i128 via CRT
        let p0 = aux_basis.moduli[0];
        let p1 = aux_basis.moduli[1];
        let bk0 = bk_pj[0];
        let bk1 = bk_pj[1];

        // CRT: m = m_0 * p1 * (p1^{-1} mod p0) + m_1 * p0 * (p0^{-1} mod p1)  (mod P)
        // P = p0 * p1
        let p1_inv_p0 = mod_inv(p1 % p0, p0)
            .ok_or_else(|| ExactoError::InvalidParam("p1 not invertible mod p0".into()))?;
        let p0_inv_p1 = mod_inv(p0 % p1, p1)
            .ok_or_else(|| ExactoError::InvalidParam("p0 not invertible mod p1".into()))?;
        let big_p_128 = p0 as i128 * p1 as i128;
        let half_big_p = big_p_128 / 2;

        for i in 0..n {
            let a = a_poly.coeffs[i];
            let a_centered: i128 = if a > half_q { a as i128 - q_128 } else { a as i128 };

            // m mod p0 and m mod p1 — use centered a for consistency
            let b0 = b_polys[0].coeffs[i];
            let a_ext0 = if a > half_q {
                let neg = q - a;
                let rem = neg % p0;
                if rem == 0 { 0 } else { p0 - rem }
            } else {
                a % p0
            };
            let diff0 = if b0 >= a_ext0 { b0 - a_ext0 } else { p0 - a_ext0 + b0 };
            let m0 = mod_mul(diff0, q_inv_pj[0], p0, bk0);

            let b1 = b_polys[1].coeffs[i];
            let a_ext1 = if a > half_q {
                let neg = q - a;
                let rem = neg % p1;
                if rem == 0 { 0 } else { p1 - rem }
            } else {
                a % p1
            };
            let diff1 = if b1 >= a_ext1 { b1 - a_ext1 } else { p1 - a_ext1 + b1 };
            let m1 = mod_mul(diff1, q_inv_pj[1], p1, bk1);

            // CRT reconstruct m in [0, P) using i128
            // t0 = m0 * p1_inv_p0 mod p0
            // t1 = m1 * p0_inv_p1 mod p1
            // m_CRT = t0 * p1 + t1 * p0  (mod P)
            let t0 = mod_mul(m0, p1_inv_p0, p0, bk0) as i128;
            let t1 = mod_mul(m1, p0_inv_p1, p1, bk1) as i128;
            let crt_sum = t0 * p1 as i128 + t1 * p0 as i128;
            let m_crt = crt_sum % big_p_128; // in [0, P)

            // Center: m ∈ [-P/2, P/2), then reduce mod q to avoid i128 overflow
            // (p * m_centered can exceed i128 when P ≈ 2^110)
            let m_centered = if m_crt > half_big_p { m_crt - big_p_128 } else { m_crt };
            let m_mod_q = ((m_centered % q_128) + q_128) % q_128;

            let pa = p_128 * a_centered;
            let round_pa_q = if pa >= 0 {
                (pa + q_128 / 2) / q_128
            } else {
                -((-pa + q_128 / 2) / q_128)
            };

            // Compute result = round(p*a/q) + p*m (mod q)
            let round_mod_q = ((round_pa_q % q_128 + q_128) % q_128) as u64;
            let pm_mod_q = mod_mul(p, m_mod_q as u64, q, barrett_constant(q));
            result[i] = {
                let sum = round_mod_q as u128 + pm_mod_q as u128;
                (sum % q as u128) as u64
            };
        }
    } else {
        return Err(ExactoError::InvalidParam(
            format!("HPS scaling supports 1 or 2 aux primes, got {}", num_aux)
        ));
    }

    let result_poly = CoeffPoly { coeffs: result, modulus: q };
    RnsPoly::from_coeff_poly(&result_poly, ct_basis)
}

/// O(n^2) schoolbook fallback for BFV multiplication (used when no aux_basis).
fn bfv_mul_schoolbook(
    ct1: &BfvCiphertext,
    ct2: &BfvCiphertext,
) -> Result<BfvCiphertext> {
    let params = &ct1.params;
    let basis = &params.ct_basis;
    let p = params.plain_modulus;
    let q = basis.moduli[0];
    let n = params.ring_degree;

    if schoolbook_overflow_risk(p, q, n) {
        return Err(ExactoError::NotImplemented(
            "schoolbook BFV multiplication can overflow i128 for these parameters; use HPS auxiliary basis"
                .into(),
        ));
    }

    let c0 = to_centered(&ct1.c[0].to_coeff_poly(basis), q);
    let c1 = to_centered(&ct1.c[1].to_coeff_poly(basis), q);
    let d0 = to_centered(&ct2.c[0].to_coeff_poly(basis), q);
    let d1 = to_centered(&ct2.c[1].to_coeff_poly(basis), q);

    let t0 = poly_mul_i128(&c0, &d0, n);
    let t1 = poly_add_i128(&poly_mul_i128(&c0, &d1, n), &poly_mul_i128(&c1, &d0, n));
    let t2 = poly_mul_i128(&c1, &d1, n);

    let r0 = scale_tensor_component(&t0, p, q);
    let r1 = scale_tensor_component(&t1, p, q);
    let r2 = scale_tensor_component(&t2, p, q);

    let r0_rns = RnsPoly::from_coeff_poly(&CoeffPoly::from_coeffs(r0, q), basis)?;
    let r1_rns = RnsPoly::from_coeff_poly(&CoeffPoly::from_coeffs(r1, q), basis)?;
    let r2_rns = RnsPoly::from_coeff_poly(&CoeffPoly::from_coeffs(r2, q), basis)?;

    Ok(BfvCiphertext {
        c: vec![r0_rns, r1_rns, r2_rns],
        params: params.clone(),
    })
}

/// Conservative i128-overflow check for schoolbook multiplication + scaling.
fn schoolbook_overflow_risk(p: u64, q: u64, n: usize) -> bool {
    let max_coeff = (q / 2) as u128;
    let max_tensor = (n as u128)
        .saturating_mul(max_coeff)
        .saturating_mul(max_coeff);
    let max_scaled = max_tensor.saturating_mul(p as u128);
    max_tensor > i128::MAX as u128 || max_scaled > i128::MAX as u128
}

/// Multiply ciphertext by a plaintext polynomial.
/// ct_out = ct * pt (component-wise multiplication by Δ·pt).
pub fn bfv_plain_mul(
    ct: &BfvCiphertext,
    plaintext: &CoeffPoly,
) -> Result<BfvCiphertext> {
    let params = &ct.params;
    let basis = &params.ct_basis;

    // Convert plaintext to RNS-NTT (without Δ scaling — raw plaintext multiply)
    let pt_rns = RnsPoly::from_coeff_poly(plaintext, basis)?;

    let c: Vec<RnsPoly> = ct.c.iter()
        .map(|ci| ci.mul(&pt_rns))
        .collect::<Result<Vec<_>>>()?;

    Ok(BfvCiphertext {
        c,
        params: params.clone(),
    })
}

/// Add a plaintext to a ciphertext: ct + Δ·m.
pub fn bfv_plain_add(
    ct: &BfvCiphertext,
    plaintext: &CoeffPoly,
) -> Result<BfvCiphertext> {
    let params = &ct.params;
    let delta_m = scale_plaintext(plaintext, params)?;

    let mut c = ct.c.clone();
    c[0] = c[0].add(&delta_m)?;

    Ok(BfvCiphertext {
        c,
        params: params.clone(),
    })
}

/// Apply a Galois automorphism σ_k: X → X^k to a BFV ciphertext.
///
/// Given ct = (c0, c1) encrypting μ under secret key s:
///   1. Compute c0' = σ_k(c0), c1_temp = σ_k(c1)
///   2. Key-switch c1_temp from σ_k(s) to s using the Galois key
///
/// The result encrypts σ_k(μ) under the original key s.
pub fn bfv_apply_automorphism(
    ct: &BfvCiphertext,
    gk: &GaloisKey,
) -> Result<BfvCiphertext> {
    if ct.c.len() != 2 {
        return Err(ExactoError::InvalidParam(
            "automorphism requires degree-1 ciphertext".into()
        ));
    }

    let params = &ct.params;
    let basis = &params.ct_basis;
    let k = gk.element;

    // Apply automorphism to c0 and c1 (convert to coeffs, permute, convert back)
    let c0_coeffs = ct.c[0].to_coeff_poly(basis);
    let c1_coeffs = ct.c[1].to_coeff_poly(basis);

    let c0_auto = crate::bfv::keygen::apply_automorphism(&c0_coeffs, k);
    let c1_auto = crate::bfv::keygen::apply_automorphism(&c1_coeffs, k);

    let c0_new_init = RnsPoly::from_coeff_poly(&c0_auto, basis)?;

    // Key-switch c1_auto from σ_k(s) to s:
    // c0_final = c0_auto + Σ digit_i * gk[i].0
    // c1_final = Σ digit_i * gk[i].1
    let digits = gadget_decompose(&c1_auto, params.gadget_base, params.gadget_digits);

    let mut c0_final = c0_new_init;
    let mut c1_final: Option<RnsPoly> = None;

    for (i, digit) in digits.iter().enumerate() {
        if i >= gk.keys.len() {
            break;
        }
        let digit_rns = RnsPoly::from_coeff_poly(digit, basis)?;
        let prod0 = digit_rns.mul(&gk.keys[i].0)?;
        let prod1 = digit_rns.mul(&gk.keys[i].1)?;
        c0_final = c0_final.add(&prod0)?;
        c1_final = Some(match c1_final {
            Some(acc) => acc.add(&prod1)?,
            None => prod1,
        });
    }

    Ok(BfvCiphertext {
        c: vec![c0_final, c1_final.unwrap()],
        params: params.clone(),
    })
}

/// Compute the homomorphic trace: Tr(ct) = Σ_{i ∈ G} σ_i(ct).
///
/// For the power-of-2 cyclotomic ring R = Z[X]/(X^n+1), the Galois group is
/// {σ_k : k odd, 1 ≤ k < 2n}. The trace sums over a subgroup, extracting
/// specific coefficients of the plaintext polynomial.
///
/// This is used in CoeffsToSlots to isolate individual polynomial coefficients.
///
/// `galois_elements` specifies the elements k for which to apply σ_k and sum.
/// `galois_keys` maps element → GaloisKey.
pub fn bfv_trace(
    ct: &BfvCiphertext,
    galois_elements: &[usize],
    galois_keys: &std::collections::HashMap<usize, GaloisKey>,
) -> Result<BfvCiphertext> {
    let mut result = ct.clone();
    for &k in galois_elements {
        let gk = galois_keys.get(&k).ok_or_else(|| {
            ExactoError::InvalidParam(format!("missing Galois key for element {k}"))
        })?;
        let rotated = bfv_apply_automorphism(&result, gk)?;
        result = bfv_add(&result, &rotated)?;
    }
    Ok(result)
}

/// Inner product of a vector of ciphertexts with plaintext coefficients.
///
/// Computes Σ_i pt_i * ct_i. Used in homomorphic linear transforms
/// (e.g., CoeffsToSlots diagonal method).
pub fn bfv_inner_product(
    cts: &[BfvCiphertext],
    pts: &[CoeffPoly],
) -> Result<BfvCiphertext> {
    if cts.is_empty() || cts.len() != pts.len() {
        return Err(ExactoError::InvalidParam("mismatched ct/pt lengths".into()));
    }
    let mut acc = bfv_plain_mul(&cts[0], &pts[0])?;
    for i in 1..cts.len() {
        let term = bfv_plain_mul(&cts[i], &pts[i])?;
        acc = bfv_add(&acc, &term)?;
    }
    Ok(acc)
}

/// Multiply a BFV ciphertext by the monomial X^j in R = Z[X]/(X^n+1).
///
/// This shifts the plaintext polynomial coefficients: if ct encrypts m(X),
/// the result encrypts X^j · m(X) mod (X^n+1). Since X^n = -1, this is
/// a cyclic rotation with negation of wrapped coefficients.
pub fn bfv_monomial_mul(ct: &BfvCiphertext, j: usize) -> Result<BfvCiphertext> {
    let params = &ct.params;
    let basis = &params.ct_basis;
    let n = params.ring_degree;
    let j = j % (2 * n); // reduce mod 2n

    if j == 0 {
        return Ok(ct.clone());
    }

    let c: Vec<RnsPoly> = ct.c.iter()
        .map(|ci| {
            let coeffs = ci.to_coeff_poly(basis);
            let rotated = monomial_mul_poly(&coeffs, j, n);
            RnsPoly::from_coeff_poly(&rotated, basis)
        })
        .collect::<Result<Vec<_>>>()?;

    Ok(BfvCiphertext { c, params: params.clone() })
}

/// Multiply a coefficient polynomial by X^j mod (X^n+1).
fn monomial_mul_poly(poly: &CoeffPoly, j: usize, n: usize) -> CoeffPoly {
    let q = poly.modulus;
    let mut result = vec![0u64; n];

    for (i, &c) in poly.coeffs.iter().enumerate() {
        if c == 0 { continue; }
        let new_idx = (i + j) % (2 * n);
        if new_idx < n {
            result[new_idx] = crate::ring::modular::mod_add(result[new_idx], c, q);
        } else {
            // X^{n+k} = -X^k
            let k = new_idx - n;
            result[k] = crate::ring::modular::mod_sub(result[k], c, q);
        }
    }

    CoeffPoly { coeffs: result, modulus: q }
}

// --- Internal helpers ---

/// Convert a polynomial from [0, q) representation to centered [-q/2, q/2) as i128.
fn to_centered(poly: &CoeffPoly, q: u64) -> Vec<i128> {
    let half_q = q / 2;
    poly.coeffs.iter().map(|&c| {
        if c > half_q {
            c as i128 - q as i128
        } else {
            c as i128
        }
    }).collect()
}

/// Polynomial multiply in Z[X]/(X^n+1) using i128 arithmetic (exact, no modular reduction).
fn poly_mul_i128(a: &[i128], b: &[i128], n: usize) -> Vec<i128> {
    let mut result = vec![0i128; n];
    for i in 0..n {
        if a[i] == 0 { continue; }
        for j in 0..n {
            if b[j] == 0 { continue; }
            let prod = a[i] * b[j];
            let idx = i + j;
            if idx < n {
                result[idx] += prod;
            } else {
                // X^n ≡ -1 in X^n+1
                result[idx - n] -= prod;
            }
        }
    }
    result
}

/// Add two i128 polynomial vectors.
fn poly_add_i128(a: &[i128], b: &[i128]) -> Vec<i128> {
    a.iter().zip(b.iter()).map(|(&x, &y)| x + y).collect()
}

/// Scale a tensor product component: for each coefficient x (an integer over Z),
/// compute round(p * x / q) mod q.
fn scale_tensor_component(coeffs: &[i128], p: u64, q: u64) -> Vec<u64> {
    let p = p as i128;
    let q_i = q as i128;
    coeffs.iter().map(|&x| {
        // round(p * x / q) with correct rounding for negative values
        let numerator = p * x;
        let rounded = if numerator >= 0 {
            (numerator + q_i / 2) / q_i
        } else {
            -(-numerator + q_i / 2) / q_i
        };
        // Reduce mod q to [0, q)
        ((rounded % q_i + q_i) % q_i) as u64
    }).collect()
}

fn modulus_product_biguint(moduli: &[u64]) -> BigUint {
    let mut q_big = BigUint::one();
    for &q in moduli {
        q_big *= BigUint::from(q);
    }
    q_big
}

fn reconstruct_centered_bigint(
    poly: &RnsPoly,
    basis: &RnsBasis,
    q_big_u: &BigUint,
) -> Result<Vec<BigInt>> {
    let n = poly.ring_degree;
    let coeff_components: Vec<CoeffPoly> = poly.components.iter()
        .map(|c| c.to_coeff_poly())
        .collect();

    let mut crt_terms = Vec::with_capacity(basis.moduli.len());
    for &qi in &basis.moduli {
        let qi_big = BigUint::from(qi);
        let q_star = q_big_u / &qi_big;
        let q_star_mod_qi = (&q_star % &qi_big)
            .to_u64()
            .ok_or_else(|| ExactoError::InvalidParam(
                "failed to reduce CRT factor modulo q_i".into()
            ))?;
        let inv = mod_inv(q_star_mod_qi, qi)
            .ok_or_else(|| ExactoError::InvalidParam("non-coprime ciphertext moduli".into()))?;
        crt_terms.push(q_star * BigUint::from(inv));
    }

    let q_big = BigInt::from(q_big_u.clone());
    let half_q = &q_big >> 1;
    let mut out = vec![BigInt::zero(); n];

    for coeff_idx in 0..n {
        let mut x = BigUint::zero();
        for (i, coeffs_i) in coeff_components.iter().enumerate() {
            x += &crt_terms[i] * BigUint::from(coeffs_i.coeffs[coeff_idx]);
        }
        x %= q_big_u;

        let mut centered = BigInt::from(x);
        if centered > half_q {
            centered -= &q_big;
        }
        out[coeff_idx] = centered;
    }

    Ok(out)
}

fn centered_bigint_to_rns(
    coeffs: &[BigInt],
    basis: &RnsBasis,
) -> Result<RnsPoly> {
    let n = basis.ring_degree;
    let mut components = Vec::with_capacity(basis.num_moduli());

    for (idx, (&qi, plan)) in basis.moduli.iter().zip(basis.plans.iter()).enumerate() {
        let _ = idx;
        let qi_big = BigInt::from(qi);
        let mut coeffs_qi = vec![0u64; n];
        for (j, c) in coeffs.iter().enumerate() {
            let mut r = c % &qi_big;
            if r.is_negative() {
                r += &qi_big;
            }
            coeffs_qi[j] = r.to_u64().ok_or_else(|| {
                ExactoError::InvalidParam("coefficient conversion to u64 failed".into())
            })?;
        }
        let cp = CoeffPoly { coeffs: coeffs_qi, modulus: qi };
        components.push(NttPoly::from_coeff_poly(&cp, plan.clone())?);
    }

    Ok(RnsPoly {
        components,
        ring_degree: basis.ring_degree,
    })
}

fn poly_mul_bigint(a: &[BigInt], b: &[BigInt], n: usize) -> Vec<BigInt> {
    let mut result = vec![BigInt::zero(); n];
    for i in 0..n {
        if a[i].is_zero() { continue; }
        for j in 0..n {
            if b[j].is_zero() { continue; }
            let prod = &a[i] * &b[j];
            let idx = i + j;
            if idx < n {
                result[idx] += &prod;
            } else {
                result[idx - n] -= &prod;
            }
        }
    }
    result
}

fn poly_add_bigint(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
}

fn scale_tensor_component_bigint(
    coeffs: &[BigInt],
    p: u64,
    q_big: &BigInt,
    half_q: &BigInt,
) -> Vec<BigInt> {
    let p_big = BigInt::from(p);
    coeffs.iter().map(|x| {
        let num = &p_big * x;
        if num.is_negative() {
            -(((-num) + half_q) / q_big)
        } else {
            (num + half_q) / q_big
        }
    }).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bfv::keygen::*;
    use crate::bfv::encrypt::*;
    use crate::bfv::encoding::*;
    use crate::params::BfvParamsBuilder;
    use crate::params::presets::compact_bfv;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_homomorphic_add() {
        let params = compact_bfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let sk = gen_secret_key_with_rng(&params, &mut rng).unwrap();

        let pt1 = encode_scalar(10, &params).unwrap();
        let pt2 = encode_scalar(20, &params).unwrap();

        let ct1 = encrypt_sk_with_rng(&pt1, &sk, &params, &mut rng).unwrap();
        let ct2 = encrypt_sk_with_rng(&pt2, &sk, &params, &mut rng).unwrap();

        let ct_sum = bfv_add(&ct1, &ct2).unwrap();
        let decrypted = decrypt(&ct_sum, &sk).unwrap();

        assert_eq!(decode_scalar(&decrypted), 30);
    }

    #[test]
    fn test_homomorphic_sub() {
        let params = compact_bfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let sk = gen_secret_key_with_rng(&params, &mut rng).unwrap();

        let pt1 = encode_scalar(50, &params).unwrap();
        let pt2 = encode_scalar(20, &params).unwrap();

        let ct1 = encrypt_sk_with_rng(&pt1, &sk, &params, &mut rng).unwrap();
        let ct2 = encrypt_sk_with_rng(&pt2, &sk, &params, &mut rng).unwrap();

        let ct_diff = bfv_sub(&ct1, &ct2).unwrap();
        let decrypted = decrypt(&ct_diff, &sk).unwrap();

        assert_eq!(decode_scalar(&decrypted), 30);
    }

    #[test]
    fn test_homomorphic_mul() {
        let params = compact_bfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let sk = gen_secret_key_with_rng(&params, &mut rng).unwrap();
        let rlk = gen_relin_key_with_rng(&sk, &mut rng).unwrap();

        let pt1 = encode_scalar(3, &params).unwrap();
        let pt2 = encode_scalar(7, &params).unwrap();

        let ct1 = encrypt_sk_with_rng(&pt1, &sk, &params, &mut rng).unwrap();
        let ct2 = encrypt_sk_with_rng(&pt2, &sk, &params, &mut rng).unwrap();

        let ct_prod = bfv_mul_and_relin(&ct1, &ct2, &rlk).unwrap();
        let decrypted = decrypt(&ct_prod, &sk).unwrap();

        assert_eq!(decode_scalar(&decrypted), 21);
    }

    #[test]
    fn test_homomorphic_mul_multi_prime_q() {
        let params = BfvParamsBuilder::new()
            .ring_degree(16)
            .plain_modulus(257)
            .ct_moduli(vec![65537, 1099509805057]) // both ≡ 1 mod 32
            .sigma(3.2)
            .gadget_base(8)
            .build()
            .unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(1234);

        let sk = gen_secret_key_with_rng(&params, &mut rng).unwrap();
        let rlk = gen_relin_key_with_rng(&sk, &mut rng).unwrap();

        for (a, b, expected) in [(3u64, 7u64, 21u64), (10, 20, 200), (0, 5, 0)] {
            let pt1 = encode_scalar(a, &params).unwrap();
            let pt2 = encode_scalar(b, &params).unwrap();
            let ct1 = encrypt_sk_with_rng(&pt1, &sk, &params, &mut rng).unwrap();
            let ct2 = encrypt_sk_with_rng(&pt2, &sk, &params, &mut rng).unwrap();

            let ct_prod = bfv_mul_and_relin(&ct1, &ct2, &rlk).unwrap();
            let decrypted = decrypt(&ct_prod, &sk).unwrap();
            assert_eq!(decode_scalar(&decrypted), expected);
        }
    }

    #[test]
    fn test_apply_automorphism() {
        let params = compact_bfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let sk = gen_secret_key_with_rng(&params, &mut rng).unwrap();

        // Encrypt a polynomial with known coefficients
        // Use scalar encoding: m = 10 (constant polynomial)
        let pt = encode_scalar(10, &params).unwrap();
        let ct = encrypt_sk_with_rng(&pt, &sk, &params, &mut rng).unwrap();

        // Apply σ_3: X → X^3. For a constant polynomial (all coefficients in slot 0),
        // the automorphism should preserve the value.
        let gk = crate::bfv::keygen::gen_galois_key_with_rng(&sk, 3, &mut rng).unwrap();
        let ct_auto = bfv_apply_automorphism(&ct, &gk).unwrap();

        let decrypted = decrypt(&ct_auto, &sk).unwrap();
        // The automorphism permutes the plaintext polynomial.
        // For a constant (scalar-encoded) polynomial, σ_k(m) = m.
        assert_eq!(decode_scalar(&decrypted), 10,
            "automorphism should preserve scalar value");
    }

    #[test]
    fn test_apply_automorphism_poly() {
        let params = compact_bfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let sk = gen_secret_key_with_rng(&params, &mut rng).unwrap();

        // Encrypt m(X) = 1 + 2X (polynomial with two nonzero coefficients)
        let mut pt_coeffs = vec![0u64; params.ring_degree];
        pt_coeffs[0] = 1;
        pt_coeffs[1] = 2;
        let pt = CoeffPoly { coeffs: pt_coeffs, modulus: params.plain_modulus };
        let ct = encrypt_sk_with_rng(&pt, &sk, &params, &mut rng).unwrap();

        // Apply σ_3: X → X^3. Result should be 1 + 2X^3.
        let gk = crate::bfv::keygen::gen_galois_key_with_rng(&sk, 3, &mut rng).unwrap();
        let ct_auto = bfv_apply_automorphism(&ct, &gk).unwrap();

        let decrypted = decrypt(&ct_auto, &sk).unwrap();
        assert_eq!(decrypted.coeffs[0], 1, "constant term should be 1");
        assert_eq!(decrypted.coeffs[1], 0, "X^1 coeff should be 0");
        assert_eq!(decrypted.coeffs[2], 0, "X^2 coeff should be 0");
        assert_eq!(decrypted.coeffs[3], 2, "X^3 coeff should be 2");
    }

    #[test]
    fn test_plain_add() {
        let params = compact_bfv().unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(42);

        let sk = gen_secret_key_with_rng(&params, &mut rng).unwrap();

        let pt1 = encode_scalar(10, &params).unwrap();
        let pt2 = encode_scalar(5, &params).unwrap();

        let ct = encrypt_sk_with_rng(&pt1, &sk, &params, &mut rng).unwrap();
        let ct_sum = bfv_plain_add(&ct, &pt2).unwrap();
        let decrypted = decrypt(&ct_sum, &sk).unwrap();

        assert_eq!(decode_scalar(&decrypted), 15);
    }
}
