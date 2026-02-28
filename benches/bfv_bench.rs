use criterion::{criterion_group, criterion_main, Criterion, black_box};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

use exacto::params::presets::{compact_bfv, compact_dbfv};
use exacto::bfv::keygen::{gen_secret_key_with_rng, gen_public_key_with_rng, gen_relin_key_with_rng, gen_galois_key_with_rng};
use exacto::bfv::encrypt::{encrypt_sk_with_rng, encrypt_pk_with_rng, decrypt};
use exacto::bfv::encoding::encode_scalar;
use exacto::bfv::eval::{bfv_add, bfv_mul_and_relin, bfv_plain_mul};
use exacto::ring::poly::CoeffPoly;

fn bfv_keygen(c: &mut Criterion) {
    let params = compact_bfv().unwrap();
    let mut rng = ChaCha20Rng::seed_from_u64(0);

    c.bench_function("bfv_keygen_secret", |b| {
        b.iter(|| gen_secret_key_with_rng(black_box(&params), &mut rng))
    });

    let sk = gen_secret_key_with_rng(&params, &mut rng).unwrap();
    c.bench_function("bfv_keygen_public", |b| {
        b.iter(|| gen_public_key_with_rng(black_box(&sk), &mut rng))
    });

    c.bench_function("bfv_keygen_relin", |b| {
        b.iter(|| gen_relin_key_with_rng(black_box(&sk), &mut rng))
    });
}

fn bfv_encrypt_decrypt(c: &mut Criterion) {
    let params = compact_bfv().unwrap();
    let mut rng = ChaCha20Rng::seed_from_u64(1);
    let sk = gen_secret_key_with_rng(&params, &mut rng).unwrap();
    let pk = gen_public_key_with_rng(&sk, &mut rng).unwrap();
    let pt = encode_scalar(42, &params).unwrap();

    c.bench_function("bfv_encrypt_sk", |b| {
        b.iter(|| encrypt_sk_with_rng(black_box(&pt), &sk, &params, &mut rng))
    });

    c.bench_function("bfv_encrypt_pk", |b| {
        b.iter(|| encrypt_pk_with_rng(black_box(&pt), &pk, &params, &mut rng))
    });

    let ct = encrypt_sk_with_rng(&pt, &sk, &params, &mut rng).unwrap();
    c.bench_function("bfv_decrypt", |b| {
        b.iter(|| decrypt(black_box(&ct), &sk))
    });
}

fn bfv_eval(c: &mut Criterion) {
    let params = compact_bfv().unwrap();
    let mut rng = ChaCha20Rng::seed_from_u64(2);
    let sk = gen_secret_key_with_rng(&params, &mut rng).unwrap();
    let rlk = gen_relin_key_with_rng(&sk, &mut rng).unwrap();
    let pt1 = encode_scalar(10, &params).unwrap();
    let pt2 = encode_scalar(20, &params).unwrap();
    let ct1 = encrypt_sk_with_rng(&pt1, &sk, &params, &mut rng).unwrap();
    let ct2 = encrypt_sk_with_rng(&pt2, &sk, &params, &mut rng).unwrap();

    c.bench_function("bfv_add", |b| {
        b.iter(|| bfv_add(black_box(&ct1), black_box(&ct2)))
    });

    c.bench_function("bfv_mul_and_relin", |b| {
        b.iter(|| bfv_mul_and_relin(black_box(&ct1), black_box(&ct2), &rlk))
    });

    let plain = CoeffPoly {
        coeffs: {
            let mut v = vec![0u64; params.ring_degree];
            v[0] = 3;
            v
        },
        modulus: params.plain_modulus,
    };
    c.bench_function("bfv_plain_mul", |b| {
        b.iter(|| bfv_plain_mul(black_box(&ct1), black_box(&plain)))
    });
}

fn dbfv_operations(c: &mut Criterion) {
    let dbfv_params = compact_dbfv().unwrap();
    let bfv_params = &dbfv_params.bfv_params;
    let mut rng = ChaCha20Rng::seed_from_u64(3);
    let sk = gen_secret_key_with_rng(bfv_params, &mut rng).unwrap();
    let rlk = gen_relin_key_with_rng(&sk, &mut rng).unwrap();

    use exacto::dbfv::encrypt::dbfv_encrypt_sk;
    use exacto::dbfv::decrypt::dbfv_decrypt;
    use exacto::dbfv::advanced::{dbfv_apply_automorphism, dbfv_change_base, dbfv_div_by_base};
    use exacto::dbfv::eval::{dbfv_add, dbfv_mul};

    let ct1 = dbfv_encrypt_sk(100, &sk, &dbfv_params).unwrap();
    let ct2 = dbfv_encrypt_sk(50, &sk, &dbfv_params).unwrap();
    let ct_div = dbfv_encrypt_sk(128, &sk, &dbfv_params).unwrap();
    let gk = gen_galois_key_with_rng(&sk, 3, &mut rng).unwrap();

    c.bench_function("dbfv_encrypt", |b| {
        b.iter(|| dbfv_encrypt_sk(black_box(42), &sk, &dbfv_params))
    });

    c.bench_function("dbfv_decrypt", |b| {
        b.iter(|| dbfv_decrypt(black_box(&ct1), &sk))
    });

    c.bench_function("dbfv_add", |b| {
        b.iter(|| dbfv_add(black_box(&ct1), black_box(&ct2)))
    });

    c.bench_function("dbfv_mul", |b| {
        b.iter(|| dbfv_mul(black_box(&ct1), black_box(&ct2), &rlk))
    });

    c.bench_function("dbfv_apply_automorphism", |b| {
        b.iter(|| dbfv_apply_automorphism(black_box(&ct1), &gk))
    });

    c.bench_function("dbfv_div_by_base", |b| {
        b.iter(|| dbfv_div_by_base(black_box(&ct_div)))
    });

    c.bench_function("dbfv_change_base_16_to_4", |b| {
        b.iter(|| dbfv_change_base(black_box(&ct1), 4, 4))
    });
}

criterion_group!(benches, bfv_keygen, bfv_encrypt_decrypt, bfv_eval, dbfv_operations);
criterion_main!(benches);
