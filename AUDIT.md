# Exacto Consistency Audit (Current)

Audit baseline: Peikert-Zarchy-Zyskind (2026), “High-Precision Exact FHE Made Simple, General, and Fast.”

## Audit Verdict

- Open blocking findings: **none**
- Test status: `cargo test -q` passes (`99` total)
- Reproduction harness: `cargo run --release --bin paper_repro` generates
  `reports/paper_reproduction.md`

## Implemented and Verified

- BFV core: keygen, encrypt/decrypt, add/sub/neg, mul+relin, plain mul/add, automorphisms.
- dBFV core: digit decomposition, encrypt/decrypt (scalar + ring plaintext), add/sub/neg, mul.
- dBFV advanced ops: automorphism wrapper, division by base, homomorphic base change.
- dBFV safety path for deep chains: `dbfv_mul_then_bootstrap` and
  `dbfv_mul_chain_then_bootstrap`.
- Bootstrapping host: BFV bootstrap + dBFV limb bootstrap flow.
- CoeffsToSlots: optimized relative-trace path for larger power-of-two rings with
  correctness-preserving fallback for small rings.
- Multi-prime ciphertext modulus support in BFV multiplication path.

## Scope Notes (Not Blocking)

The following are explicit scope boundaries for this codebase version and are not tracked
as defects:

- dGSW/FHEW-style bootstrapping host.
- Specialized “variant LWE” encryption path from paper Method 2.
- Full paper-wide feature parity outside the BFV/dBFV + bootstrap host scope.

These are roadmap expansions, not regressions against supported functionality.
