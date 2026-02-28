use crate::error::{ExactoError, Result};
use crate::ring::poly::CoeffPoly;
use crate::params::BfvParams;

/// Encode a scalar value into a plaintext polynomial.
/// The scalar m is placed in all n coefficients (constant polynomial m).
pub fn encode_scalar(m: u64, params: &BfvParams) -> Result<CoeffPoly> {
    if m >= params.plain_modulus {
        return Err(ExactoError::InvalidParam(
            format!("plaintext {} >= plain_modulus {}", m, params.plain_modulus)
        ));
    }
    let mut coeffs = vec![0u64; params.ring_degree];
    coeffs[0] = m;
    Ok(CoeffPoly {
        coeffs,
        modulus: params.plain_modulus,
    })
}

/// Decode a scalar from a plaintext polynomial (read coefficient 0).
pub fn decode_scalar(poly: &CoeffPoly) -> u64 {
    poly.coeffs[0]
}

/// Encode a vector of values using SIMD/CRT packing.
///
/// For plaintext modulus p prime and p ≡ 1 (mod 2n), the ring Z_p[X]/(X^n+1)
/// splits into n slots via CRT. Each slot can hold an independent Z_p value.
///
/// For now, this is a simple coefficient packing (not true SIMD).
/// True SIMD requires computing roots of X^n+1 mod p.
pub fn encode_simd(values: &[u64], params: &BfvParams) -> Result<CoeffPoly> {
    if values.len() > params.ring_degree {
        return Err(ExactoError::DimensionMismatch {
            expected: params.ring_degree,
            got: values.len(),
        });
    }
    for &v in values {
        if v >= params.plain_modulus {
            return Err(ExactoError::InvalidParam(
                format!("plaintext {} >= plain_modulus {}", v, params.plain_modulus)
            ));
        }
    }
    let mut coeffs = vec![0u64; params.ring_degree];
    coeffs[..values.len()].copy_from_slice(values);
    Ok(CoeffPoly {
        coeffs,
        modulus: params.plain_modulus,
    })
}

/// Decode SIMD-packed values.
pub fn decode_simd(poly: &CoeffPoly, num_slots: usize) -> Vec<u64> {
    poly.coeffs[..num_slots].to_vec()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::presets::compact_bfv;

    #[test]
    fn test_encode_decode_scalar() {
        let params = compact_bfv().unwrap();
        let pt = encode_scalar(42, &params).unwrap();
        assert_eq!(decode_scalar(&pt), 42);
    }

    #[test]
    fn test_encode_simd() {
        let params = compact_bfv().unwrap();
        let values = vec![1, 2, 3, 4, 5];
        let pt = encode_simd(&values, &params).unwrap();
        let decoded = decode_simd(&pt, 5);
        assert_eq!(decoded, values);
    }

    #[test]
    fn test_encode_out_of_range() {
        let params = compact_bfv().unwrap();
        assert!(encode_scalar(300, &params).is_err()); // p=257
    }
}
