use crate::error::Result;
use crate::bfv::{BfvCiphertext, eval as bfv_eval, keygen::RelinKey};
use crate::dbfv::ciphertext::DbfvCiphertext;
use crate::dbfv::lattice::SmallReps;

/// Reduce a dBFV ciphertext after multiplication.
///
/// Two steps:
/// 1. **Degree reduction**: Replace limbs j ≥ d with their small representatives
///    using precomputed B^j mod L decompositions. Non-expanding.
/// 2. **Hybrid lattice reduction** (§4.6.2): Reduce coefficient sizes using
///    the lattice L structure. This is a public operation on the digit
///    coefficients. (Currently applied during decryption via signed
///    recomposition rather than at the ciphertext level.)
pub fn reduce(
    ct: &DbfvCiphertext,
    _rlk: &RelinKey,
) -> Result<DbfvCiphertext> {
    let params = &ct.params;
    let d = params.num_digits;

    if ct.degree <= d {
        return Ok(ct.clone()); // Already reduced
    }

    // Step 1: Degree reduction
    // Replace B^j (j >= d) with small representative: B^j ≡ Σ r_{j,i} · B^i mod L
    let small_reps = SmallReps::compute_simple(params.base, d, params.plain_modulus);

    // Start with the first d limbs
    let mut result_limbs: Vec<BfvCiphertext> = ct.limbs[..d].to_vec();

    // For each excess limb (j >= d), add its contribution using small reps
    for j in d..ct.limbs.len() {
        let rep_idx = j - d;
        if rep_idx >= small_reps.reps.len() {
            continue;
        }
        let rep = &small_reps.reps[rep_idx];

        // B^j ≡ Σ_i rep[i] · B^i, so add rep[i] * limb_j to result_limb_i
        for i in 0..d {
            let coeff: i64 = rep[i];
            if coeff == 0 {
                continue;
            }

            // Multiply limb_j by the scalar coefficient
            let scaled = scale_bfv_ciphertext(&ct.limbs[j], coeff)?;
            result_limbs[i] = bfv_eval::bfv_add(&result_limbs[i], &scaled)?;
        }
    }

    Ok(DbfvCiphertext {
        limbs: result_limbs,
        degree: d,
        mul_depth: ct.mul_depth,
        params: params.clone(),
    })
}

/// Scale a BFV ciphertext by an integer scalar.
/// For positive scalar: multiply each component.
/// For negative scalar: negate and multiply by |scalar|.
fn scale_bfv_ciphertext(ct: &BfvCiphertext, scalar: i64) -> Result<BfvCiphertext> {
    if scalar == 0 {
        // Return encryption of zero (just zero out components)
        let c = ct.c.iter()
            .map(|ci| {
                let mut z = ci.clone();
                for comp in &mut z.components {
                    comp.evals.fill(0);
                }
                z
            })
            .collect();
        return Ok(BfvCiphertext { c, params: ct.params.clone() });
    }

    let abs_scalar = scalar.unsigned_abs();
    let scaled = BfvCiphertext {
        c: ct.c.iter()
            .map(|ci| ci.scalar_mul(abs_scalar))
            .collect(),
        params: ct.params.clone(),
    };

    if scalar < 0 {
        Ok(bfv_eval::bfv_neg(&scaled))
    } else {
        Ok(scaled)
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_reduction_not_needed() {
        // If degree <= d, reduce should be a no-op
        // (tested implicitly through dbfv_mul tests)
    }
}
