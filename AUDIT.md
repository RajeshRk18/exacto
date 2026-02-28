# Consistency Audit: Exacto vs. Paper

Audit of the implementation against Peikert-Zarchy-Zyskind (Feb 2026),
"High-Precision Exact FHE Made Simple, General, and Fast."

## Features that Match the Paper

| Paper Section | Feature | Code Location | Notes |
|---|---|---|---|
| §3.2 | BFV secret keys & ciphertexts | `bfv/keygen.rs`, `bfv/mod.rs` | k=1 (Ring-LWE) specialization |
| §3.3 | BFV encrypt/decrypt | `bfv/encrypt.rs` | Standard Delta-scaling (equivalent to paper's Q_p) |
| §3.4.1 | BFV addition | `bfv/eval.rs` `bfv_add` | Component-wise, handles different degrees |
| §3.4.1 | BFV multiplication | `bfv/eval.rs` `bfv_mul_hps`, `bfv_mul_schoolbook` | HPS RNS (O(n log n)) + schoolbook fallback |
| §3.4.2 | Key switching / relinearization | `bfv/keyswitch.rs` | Gadget decomposition + relin key |
| §3.4.4 | BFV automorphisms | `bfv/eval.rs` `bfv_apply_automorphism` | sigma_k + key-switch from sigma_k(s) to s |
| §4.1 | Digit decomposition (Def 4.1) | `dbfv/decomposition.rs` | Scalar Z_p only, not full R/p |
| §4.2 | dBFV ciphertexts (Def 4.3) | `dbfv/ciphertext.rs` | Vec of BFV limbs, one per digit |
| §4.2.1 | dBFV decryption | `dbfv/decrypt.rs` | Rounds coefficients BEFORE evaluating at b (correct) |
| §4.3 | dBFV encryption (Method 1) | `dbfv/encrypt.rs` | Via ordinary BFV per-digit encryption |
| §4.4.1 | dBFV addition | `dbfv/eval.rs` `dbfv_add` | Limb-wise BFV addition |
| §4.4.2 | dBFV multiplication | `dbfv/eval.rs` `dbfv_mul` | d^2 BFV muls via convolution, parallelized with rayon |
| §4.6.1 | Degree reduction (p=b^d) | `dbfv/reduction.rs` | Truncation: b^j mod p = 0 for j >= d |
| §4.6.1 | Lattice basis + Babai nearest-plane | `dbfv/lattice.rs` | Operates on plaintext digit vectors |
| §5.1 | dBFV bootstrapping | `bootstrap/bfv_host.rs` `dbfv_bootstrap` | Bootstrap each BFV limb independently (parallel) |
| §5.2 | BFV bootstrapping | `bootstrap/bfv_host.rs` `bfv_bootstrap` | Modswitch + rounding poly + CoeffsToSlots |
| §5.2 | CoeffsToSlots / SlotsToCoeffs | `bootstrap/coeffs_to_slots.rs` | Trace-based coefficient extraction (naive O(n^2)) |

## Deviations (Correct but Different)

### 1. Standard Delta-scaling vs paper's scale-free Q_p

The paper eliminates explicit ciphertext modulus q and scale factor Delta, working over
K/p with fractional (rational) entries. The code uses the standard concrete BFV formulation
with `Delta = floor(q/t)`. These are mathematically equivalent: `c_paper = c_code / q`.

### 2. Scalar-only dBFV

Plaintexts are single Z_p values, not full ring elements R/p. Each dBFV ciphertext
carries one scalar decomposed across d BFV limbs. No SIMD/packed plaintexts at the
dBFV level.

### 3. Single ciphertext modulus

RNS uses a single prime Q for the ciphertext modulus. Multi-prime CRT for Q is not
supported. (Multi-prime is only used for auxiliary primes P in HPS multiplication.)

### 4. Unsigned gadget digits

Key switching uses unsigned gadget digits in [0, base) rather than balanced
[-base/2, base/2). This produces slightly more noise but is correct.

### 5. Lattice reduction deferred to decryption

Instead of ciphertext-level lattice reduction between multiplications, `digit_recompose_signed`
handles non-canonical digit representations at decrypt time. The constraint
`t > 2*d*(b-1)^2` ensures one multiplication depth is safe.

### 6. Naive trace in CoeffsToSlots

O(n) automorphisms per coefficient extraction (naive summation over Galois group) instead
of the paper's O(log n) iterative doubling down the tower of two-power cyclotomics.
Results in O(n^2) total bootstrap complexity instead of O(n log n).

### 7. Lattice basis differs from paper

The code's lattice basis last row is `[p, 0, ..., 0]` vs the paper's `b*B^{d-1} - r_d(B)`
(Equation 4.4). Both span the same lattice but produce different Gram-Schmidt
orthogonalization quality. No practical impact since ciphertext-level lattice reduction
is not applied.

## Missing Features

### Not Implemented

| Paper Section | Feature | Description |
|---|---|---|
| §4.3 Method 2 | Variant LWE encryption | Directly constructs degree < d dBFV ciphertext via variant LWE. Error rate 1/p^{1/d} vs 1/p. |
| §4.4.3 | dBFV automorphisms | Apply tau to coefficients of c(B,S). Building blocks exist (per-limb BFV auto) but no dBFV wrapper. |
| §4.4.4 | Division by divisors of p (phi_b) | Divide ciphertext's constant-in-B term by b, shift remaining degrees down. Needed for bootstrapping optimizations. |
| §4.5 | Homomorphic change of base | Convert dBFV ciphertext from base b to base b'. Useful for comparisons and base trading. |
| §4.6.2 | General hybrid lattice reduction | Paper's specialized O(d) algorithm: nearest-plane on rightmost column of R + round-off via dual matrix D*. Only standard Babai nearest-plane is implemented. |
| §4.6 | Ciphertext-level lattice reduction | Lattice reduction not applied between homomorphic multiplications (deferred to decryption). |
| §5.3 | GSW-style bootstrapping | FHEW/TFHE blind rotation approach using dGSW as host. |
| App. A | Decomposed GSW (dGSW) | Entire dGSW scheme (Def A.1, encryption, homomorphic operations) not implemented. |
| §1.3.5 | Homomorphic comparisons/sign | Iterative clipping strategy using PBS for sign extraction on decomposed plaintexts. |

### Performance Optimizations Not Implemented

| Feature | Current | Paper's approach |
|---|---|---|
| CoeffsToSlots trace | O(n) automorphisms per coeff | O(log n) via iterative doubling down cyclotomic tower |
| Gadget decomposition | Unsigned [0, base) | Balanced [-base/2, base/2) for less noise |
| dBFV encryption | BFV per-digit (error rate ~1/p) | Variant LWE (error rate ~1/p^{1/d}) |

## Correctness Concerns (Updated)

### 1. Multiplicative depth safety

Chained dBFV multiplication still requires ciphertext-level lattice reduction from §4.6.2
for full support. The implementation now **enforces** this constraint explicitly by tracking
`mul_depth` and rejecting unsupported chained multiplications with a deterministic error
instead of risking silent decryption failure.

### 2. Multi-modulus plaintext scaling

Fixed. BFV plaintext scaling now computes `Δ = floor(Q/p)` with `Q = prod(q_i)` and uses
`Δ mod q_i` in each RNS component. This replaces the prior per-modulus `floor(q_i/p)`
approximation.

### 3. HPS overflow handling and unsupported parameter guards

The existing two-aux-prime overflow handling remains correct (`m_centered` reduced modulo `q`
before multiplying by `p`). Additional guards now reject unsupported regimes explicitly:
- Single-aux HPS configs where `P <= n*Q/2` (unstable centering).
- Schoolbook multiplication configurations that can overflow `i128`.

## Test Coverage

81 tests pass, 0 ignored:
- Formerly ignored high-scale cases are now active tests that assert expected explicit errors.

## Conclusion

The core dBFV scheme (§3-4 core, §5.1-5.2) is **faithfully implemented and
algorithmically correct**. All implemented algorithms match the paper's mathematical
definitions. The missing features (dGSW, dBFV automorphisms, phi_b division, change
of base, variant LWE encryption, general hybrid reduction) are extensions and
optimizations. The foundational scheme now has explicit safety checks for unsupported
depth/parameter regimes and correct multi-modulus BFV plaintext scaling.
