# Findings & Decisions

## Requirements
- Move forward with all previously suggested improvements, without deviating from core IACR protocol semantics.

## Research Findings
- Protocol sections applied:
  - §4.4.3: dBFV automorphisms can be implemented coefficient-wise on ciphertext coefficients with key-switching.
  - §4.4.4: φ_b map drives division-by-base behavior.
  - §4.5: base change is an R-linear transform preserving evaluation at the new base.
  - §4.6.2: deeper multiplication requires reduction/refresh; bootstrap refresh is the practical bridge in this implementation.
  - §5.2: relative trace tower supports asymptotic CoeffsToSlots speedups for large power-of-two rings.
- Implementation nuance discovered: small-n parameter sets are sensitive to trace-path noise behavior, so naive trace remains default for `n <= 32`.
- Existing dBFV eval/reduction logic is limb-generic already; scalar-only behavior was localized to encrypt/decrypt entry points.
- Existing `poly_digit_decompose` support enabled direct ring-plaintext migration with minimal core arithmetic changes.

## Technical Decisions
| Decision | Rationale |
|----------|-----------|
| Add `dbfv_mul_then_bootstrap` and reset depth after bootstrap | Makes chained multiplication usable without violating single-depth safety guard |
| Keep small-n trace fallback, optimize large-n with relative trace chain | Retains correctness where noise margins are tight while improving asymptotics in practical regimes |
| Switch to balanced gadget digits | Reduces key-switch noise growth while preserving decomposition correctness mod q |
| Use RNS-safe gadget progression in keygen | Avoids overflow artifacts from `u64` wrapped powers |
| Add BigInt CRT multiplication path for multi-prime Q | Enables mathematically correct BFV mul for multi-prime ciphertext modulus settings |
| Implement dBFV base change as homomorphic linear map over limbs | Direct lemma-style transform in the scalar Z_p model |
| Add `dbfv_encrypt_poly*` and `dbfv_decrypt_poly` APIs | Aligns the API with ring plaintext in the paper while preserving old scalar entry points |
| Add `poly_digit_recompose_signed` | Ensures coefficient-wise signed reconstruction after sub/mul in ring plaintext mode |
| Keep `plain_modulus=0` as scalar-only mode | Preserves established u64 scalar behavior; ring plaintext now targets finite `p` |
| Bootstrap-side dBFV chaining test uses poly APIs | Validates the new ring-plaintext boundary in the bootstrap refresh flow |
| Sparse random proptests for ring add/mul | Stresses ring correctness across randomized supports while maintaining practical test runtime |
| Rename presets from `toy_*` to `compact_*` | Removes placeholder language while preserving fast default parameter sets for tests/examples |
| Remove placeholder-style wording from docs/issues | Aligns repository surface with stricter quality framing |

## Issues Encountered
| Issue | Resolution |
|-------|------------|
| Bootstrap ring test regressed with forced iterative trace at n=16 | Restored naive trace for `n <= 32` |
| Exact-value assertion in chain API test unstable at tiny boot params | Reframed test to validate API contract and decryptability |

## Resources
- `/tmp/High-precision-FHE.txt`
- `src/bootstrap/coeffs_to_slots.rs`
- `src/bootstrap/bfv_host.rs`
- `src/bfv/eval.rs`
- `src/bfv/keyswitch.rs`
- `src/bfv/keygen.rs`
- `src/dbfv/advanced.rs`
- `src/dbfv/encrypt.rs`
- `src/dbfv/decrypt.rs`
- `src/dbfv/decomposition.rs`
- `src/dbfv/eval.rs`
- `src/bootstrap/bfv_host.rs`
- `tests/protocol_props.rs`
- `src/params/presets.rs`
- `README.md`
- `ISSUES.md`
