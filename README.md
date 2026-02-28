# Exacto

High-precision fully homomorphic encryption based on Decomposed BFV (dBFV) from
[Peikert-Zarchy-Zyskind (2026)](https://eprint.iacr.org/2026/XXX).

Standard BFV encrypts plaintexts in Z_p and suffers noise growth proportional to p,
making large-modulus arithmetic expensive. dBFV represents each Z_p plaintext as d
digits in base b = p^{1/d} and encrypts each digit as a separate BFV ciphertext.
Noise growth drops to d * b = d * p^{1/d} -- exponentially better.

**Translation: you can multiply two encrypted u64s without the ciphertext exploding.**

## Usage

### Multiply two encrypted numbers

```rust
use exacto::params::presets::compact_dbfv;
use exacto::dbfv::keygen::dbfv_keygen;
use exacto::dbfv::{dbfv_encrypt, dbfv_decrypt, dbfv_mul};

fn main() {
    // p=256, base=16, d=2 digits -- compact preset for fast iteration
    let params = compact_dbfv().unwrap();

    // Generate secret key, public key, and relinearization key
    let (sk, pk, rlk) = dbfv_keygen(&params).unwrap();

    // Encrypt two values (plaintexts are u64 mod p=256)
    let ct_a = dbfv_encrypt(13, &pk, &params).unwrap();
    let ct_b = dbfv_encrypt(7, &pk, &params).unwrap();

    // Multiply homomorphically -- nobody sees 13 or 7
    let ct_prod = dbfv_mul(&ct_a, &ct_b, &rlk).unwrap();

    // Decrypt
    let result = dbfv_decrypt(&ct_prod, &sk).unwrap();
    assert_eq!(result, 91); // 13 * 7 = 91
}
```

### Multiply a chain safely (auto-refresh between steps)

Use `dbfv_mul_chain_then_bootstrap(cts, rlk, bsk)` to evaluate long multiplicative
chains while keeping ciphertexts in the supported depth envelope. The function applies
`dbfv_mul_then_bootstrap` at each step.

### Add and subtract

```rust
use exacto::params::presets::compact_dbfv;
use exacto::dbfv::keygen::dbfv_keygen;
use exacto::dbfv::{dbfv_encrypt, dbfv_decrypt, dbfv_add, dbfv_sub};

fn main() {
    let params = compact_dbfv().unwrap();
    let (sk, pk, _rlk) = dbfv_keygen(&params).unwrap();

    let ct_a = dbfv_encrypt(200, &pk, &params).unwrap();
    let ct_b = dbfv_encrypt(55, &pk, &params).unwrap();

    let sum = dbfv_add(&ct_a, &ct_b).unwrap();
    assert_eq!(dbfv_decrypt(&sum, &sk).unwrap(), 255); // 200 + 55

    let diff = dbfv_sub(&ct_a, &ct_b).unwrap();
    assert_eq!(dbfv_decrypt(&diff, &sk).unwrap(), 145); // 200 - 55
}
```

### Use plain BFV (single ciphertext, smaller plaintext space)

```rust
use exacto::params::presets::compact_bfv;
use exacto::bfv::keygen::{gen_secret_key, gen_public_key, gen_relin_key};
use exacto::bfv::encoding::{encode_scalar, decode_scalar};
use exacto::bfv::encrypt::{encrypt_pk, decrypt};
use exacto::bfv::eval::{bfv_add, bfv_mul_and_relin};

fn main() {
    // n=1024, p=257, 40-bit q -- compact preset for fast iteration
    let params = compact_bfv().unwrap();

    let sk = gen_secret_key(&params).unwrap();
    let pk = gen_public_key(&sk).unwrap();
    let rlk = gen_relin_key(&sk).unwrap();

    let pt_a = encode_scalar(10, &params).unwrap();
    let pt_b = encode_scalar(20, &params).unwrap();

    let ct_a = encrypt_pk(&pt_a, &pk, &params).unwrap();
    let ct_b = encrypt_pk(&pt_b, &pk, &params).unwrap();

    // Homomorphic add
    let ct_sum = bfv_add(&ct_a, &ct_b).unwrap();
    assert_eq!(decode_scalar(&decrypt(&ct_sum, &sk).unwrap()), 30);

    // Homomorphic multiply + relinearize
    let ct_prod = bfv_mul_and_relin(&ct_a, &ct_b, &rlk).unwrap();
    assert_eq!(decode_scalar(&decrypt(&ct_prod, &sk).unwrap()), 200);
}
```

### Custom parameters

```rust
use exacto::params::{BfvParamsBuilder, DbfvParams};

// Custom BFV parameters
let bfv = BfvParamsBuilder::new()
    .ring_degree(4096)
    .plain_modulus(65537)
    .ct_moduli(vec![576460752308273153]) // 60-bit NTT prime
    .sigma(3.2)
    .gadget_base(65536)                 // 2^16
    .build()
    .unwrap();

// Wrap as dBFV: p=65536, base=256, d=2
let dbfv = DbfvParams::new(bfv, 256, 2, 65536).unwrap();
```

## Architecture

```
src/
  ring/         Barrett arithmetic, CoeffPoly, NTT (concrete-ntt), RNS
  sampling/     Discrete Gaussian (constant-time CDT), uniform, ternary
  bfv/          Standard BFV: keygen, encrypt, decrypt, eval (HPS RNS mul)
  dbfv/         dBFV: digit decomposition, lattice reduction, encrypt/decrypt/eval
  bootstrap/    BFV bootstrapping: Lagrange interpolation, Paterson-Stockmeyer eval
  params/       Parameter sets, builders, presets
```

### Why dBFV is better

| Scheme | Plaintext mod | Noise growth per mul | Mul of two u64s |
|--------|--------------|---------------------|-----------------|
| BFV    | p            | ~ p                 | Painful         |
| dBFV   | p = b^d      | ~ d * b = d * p^{1/d} | Feasible     |

For p = 2^64, d = 8, b = 2^8 = 256: noise factor drops from 2^64 to 8 * 256 = 2048.
That's a factor of ~10^16 improvement.

## Benchmarks

Run with `cargo bench`. Compact preset (n=1024, 40-bit q):

| Operation         | Time    |
|-------------------|---------|
| BFV keygen (sk)   | ~80 us  |
| BFV encrypt (sk)  | ~150 us |
| BFV add           | ~4 us   |
| BFV mul + relin   | ~390 us |
| dBFV encrypt      | ~630 us |
| dBFV add          | ~4 us   |
| dBFV mul          | ~830 us |

Paper-style reproduction (u64 profiles, latency + growth + depth probe):

```bash
cargo run --release --bin paper_repro
cat reports/paper_reproduction.md
```

## Testing

```
cargo test          # full unit/integration/doctest suite
cargo bench         # criterion benchmarks
```

## Status

This is a research implementation, not production-ready. See [ISSUES.md](ISSUES.md)
for known deviations from the paper and remaining work.

### What works

- Full BFV: keygen, encrypt (pk/sk), decrypt, add, sub, mul (HPS RNS), relin,
  plain mul/add, encoding (scalar + SIMD), modulus switching, gadget decomposition
- Full dBFV: digit decomposition, encrypt (pk/sk), decrypt (signed recomposition),
  add, sub, mul (d^2 parallel BFV muls + degree reduction), lattice reduction (Babai)
- Bootstrap: key generation, modulus switching, rounding polynomial (Lagrange),
  homomorphic polynomial evaluation (Paterson-Stockmeyer), trivial/scalar bootstrap

### Future work

- Ring bootstrapping via CoeffsToSlots (Chen-Han 2018)
- Ciphertext-level lattice reduction
- Security parameter estimation
- Larger preset parameter sets matching paper benchmarks
