# Paper Reproduction Report

Command: `cargo run --release --bin paper_repro`
Scope: p=2^64 scalar dBFV mode, n=4096, CPU-only local run.

| profile | enc (ms) | add (us) | mul (ms) | growth mean | growth max | unsafe depth* |
|---|---:|---:|---:|---:|---:|---:|
| d=4, b=2^16 | 3.482 | 17.996 | 8.881 | 3701399.717 | 3870182.385 | 0 |
| d=8, b=2^8 | 5.324 | 35.828 | 31.395 | 3606163025.620 | 4195950717.538 | 1 |
| d=16, b=2^4 | 10.637 | 70.777 | 160.679 | 62062512.905 | 68028834.357 | 1 |

* `unsafe depth` is measured by bypassing the public `mul_depth` guard for experimental comparison only.
  Use `dbfv_mul_then_bootstrap` for supported chaining in application code.
