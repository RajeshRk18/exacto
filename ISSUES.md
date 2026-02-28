# Exacto Engineering Issue Log

Status: **No open blocking issues**

This file tracks resolved engineering issues that affected correctness, safety, or
benchmark reproducibility.

## Closed Issues

1. BFV plaintext scaling on multi-modulus Q
- Fixed by computing `Delta = floor(Q/p)` over full `Q = prod(q_i)` and reducing per limb.

2. Unsupported/unsafe HPS parameter regimes
- Fixed by explicit guards for single-aux instability and schoolbook overflow risk.

3. Balanced gadget decomposition and keygen gadget progression
- Fixed by centered gadget digits and overflow-safe gadget progression in key generation.

4. dBFV multiplicative-depth safety
- Fixed by depth tracking and deterministic rejection of unsupported chained `dbfv_mul`.
- Safe chaining path added: `dbfv_mul_then_bootstrap` and
  `dbfv_mul_chain_then_bootstrap`.

5. Ring-plaintext API gap in dBFV
- Fixed by `dbfv_encrypt_poly`, `dbfv_encrypt_poly_sk`, and `dbfv_decrypt_poly`.
- Scalar APIs retained as compatibility wrappers.

6. Bootstrap host consistency for dBFV
- Fixed by rebuilding dBFV params on refresh and resetting depth after bootstrap.

7. CoeffsToSlots performance/correctness envelope
- Fixed by optimized relative-trace path for larger power-of-two rings with a
  correctness fallback for small rings.

8. Benchmark reproducibility gap vs paper-style reporting
- Fixed by dedicated reproduction binary:
  `cargo run --release --bin paper_repro`
- Report artifact:
  `reports/paper_reproduction.md`

## Roadmap (Non-Issue Enhancements)

- dGSW/FHEW host path.
- Additional hardware/back-end tuning and larger profile sweeps.
- Optional expanded paper-feature coverage beyond current supported scope.
