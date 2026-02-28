pub mod security;
pub mod presets;

use std::sync::Arc;
use num_bigint::BigUint;
use num_traits::One;
use crate::error::{ExactoError, Result};
use crate::ring::rns::RnsBasis;

/// Parameters for the standard BFV scheme.
#[derive(Clone, Debug)]
pub struct BfvParams {
    /// Ring degree n (must be power of 2).
    pub ring_degree: usize,
    /// Plaintext modulus p.
    pub plain_modulus: u64,
    /// RNS basis for ciphertext modulus Q = ∏ q_i.
    pub ct_basis: Arc<RnsBasis>,
    /// Auxiliary RNS basis P for HPS multiplication (fast base extension).
    pub aux_basis: Option<Arc<RnsBasis>>,
    /// Gaussian noise standard deviation.
    pub sigma: f64,
    /// Gadget decomposition base (for key switching).
    pub gadget_base: u64,
    /// Number of gadget digits.
    pub gadget_digits: usize,
}

/// Builder for BfvParams.
pub struct BfvParamsBuilder {
    ring_degree: usize,
    plain_modulus: u64,
    ct_moduli: Vec<u64>,
    aux_moduli: Vec<u64>,
    sigma: f64,
    gadget_base: u64,
}

impl BfvParamsBuilder {
    pub fn new() -> Self {
        Self {
            ring_degree: 4096,
            plain_modulus: 65537,
            ct_moduli: Vec::new(),
            aux_moduli: Vec::new(),
            sigma: 3.2,
            gadget_base: 0, // auto-compute
        }
    }

    pub fn ring_degree(mut self, n: usize) -> Self {
        self.ring_degree = n;
        self
    }

    pub fn plain_modulus(mut self, p: u64) -> Self {
        self.plain_modulus = p;
        self
    }

    pub fn ct_moduli(mut self, moduli: Vec<u64>) -> Self {
        self.ct_moduli = moduli;
        self
    }

    pub fn aux_moduli(mut self, moduli: Vec<u64>) -> Self {
        self.aux_moduli = moduli;
        self
    }

    pub fn sigma(mut self, sigma: f64) -> Self {
        self.sigma = sigma;
        self
    }

    pub fn gadget_base(mut self, base: u64) -> Self {
        self.gadget_base = base;
        self
    }

    pub fn build(self) -> Result<Arc<BfvParams>> {
        if !self.ring_degree.is_power_of_two() || self.ring_degree < 2 {
            return Err(ExactoError::InvalidRingDegree(self.ring_degree));
        }
        if self.ct_moduli.is_empty() {
            return Err(ExactoError::InvalidParam("must specify at least one ciphertext modulus".into()));
        }
        if self.plain_modulus < 2 {
            return Err(ExactoError::InvalidParam("plaintext modulus must be >= 2".into()));
        }

        let ct_basis = Arc::new(RnsBasis::new(self.ct_moduli.clone(), self.ring_degree)?);

        let aux_basis = if self.aux_moduli.is_empty() {
            None
        } else {
            Some(Arc::new(RnsBasis::new(self.aux_moduli.clone(), self.ring_degree)?))
        };

        // Auto-compute gadget parameters if not set.
        // Digits are based on full ciphertext modulus Q = ∏ q_i.
        let (gadget_base, gadget_digits) = if self.gadget_base == 0 {
            // Use base=2^16 as specified in the plan (paper §4).
            // Relin noise ∝ d_gadget * base * σ * √n, which must be << Δ = q/p.
            // Requires sufficiently large q for correct decryption.
            let base = 1u64 << 16;
            let digits = compute_gadget_digits(&self.ct_moduli, base);
            (base, digits.max(1))
        } else {
            let digits = compute_gadget_digits(&self.ct_moduli, self.gadget_base);
            (self.gadget_base, digits.max(1))
        };

        Ok(Arc::new(BfvParams {
            ring_degree: self.ring_degree,
            plain_modulus: self.plain_modulus,
            ct_basis,
            aux_basis,
            sigma: self.sigma,
            gadget_base,
            gadget_digits,
        }))
    }
}

fn compute_gadget_digits(ct_moduli: &[u64], base: u64) -> usize {
    let mut q_big = BigUint::one();
    for &q in ct_moduli {
        q_big *= BigUint::from(q);
    }

    let base_big = BigUint::from(base);
    let mut pow = BigUint::one();
    let mut digits = 0usize;
    while pow < q_big {
        pow *= &base_big;
        digits += 1;
    }
    digits.max(1)
}

/// Parameters for the dBFV scheme.
#[derive(Clone, Debug)]
pub struct DbfvParams {
    /// Underlying BFV parameters (each digit encrypted as a BFV ciphertext).
    pub bfv_params: Arc<BfvParams>,
    /// Base b ≈ p^{1/d} for digit decomposition.
    pub base: u64,
    /// Number of digits d.
    pub num_digits: usize,
    /// Original plaintext modulus p.
    pub plain_modulus: u64,
}

impl DbfvParams {
    /// Create dBFV parameters.
    /// - `bfv_params`: BFV params with plain_modulus = base (b)
    /// - `base`: digit base b ≈ p^{1/d}
    /// - `num_digits`: d such that b^d ≥ p
    /// - `plain_modulus`: the original large plaintext modulus p.
    ///   Use 0 to represent p = 2^64 (full u64 range).
    pub fn new(
        bfv_params: Arc<BfvParams>,
        base: u64,
        num_digits: usize,
        plain_modulus: u64,
    ) -> Result<Arc<Self>> {
        if base < 2 {
            return Err(ExactoError::InvalidParam("base must be >= 2".into()));
        }
        if num_digits < 1 {
            return Err(ExactoError::InvalidParam("num_digits must be >= 1".into()));
        }
        // Verify b^d >= p (plain_modulus=0 means p=2^64)
        let mut bd = 1u128;
        for _ in 0..num_digits {
            bd = bd.saturating_mul(base as u128);
        }
        let p128 = if plain_modulus == 0 { 1u128 << 64 } else { plain_modulus as u128 };
        if bd < p128 {
            return Err(ExactoError::InvalidParam(
                format!("base^digits = {} < plain_modulus = {}", bd, p128)
            ));
        }

        Ok(Arc::new(Self {
            bfv_params,
            base,
            num_digits,
            plain_modulus,
        }))
    }
}
