use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;
use std::sync::Arc;
use std::time::{Duration, Instant};

use exacto::bfv::encrypt::decrypt as bfv_decrypt;
use exacto::bfv::keygen::{gen_relin_key_with_rng, gen_secret_key_with_rng, RelinKey, SecretKey};
use exacto::dbfv::decrypt::dbfv_decrypt;
use exacto::dbfv::encrypt::dbfv_encrypt_sk_with_rng;
use exacto::dbfv::eval::{dbfv_add, dbfv_mul};
use exacto::error::Result as ExactoResult;
use exacto::error::ExactoError;
use exacto::params::{BfvParamsBuilder, DbfvParams};
use exacto::params::presets::u64_dbfv;
use rand::Rng;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

type AppResult<T> = std::result::Result<T, Box<dyn std::error::Error>>;

#[derive(Clone)]
struct Profile {
    name: &'static str,
    base: u64,
    digits: usize,
    bfv_plain_mod: u64,
    gadget_base: u64,
}

struct ProfileResult {
    profile: Profile,
    mul_mean_ms: f64,
    add_mean_us: f64,
    enc_mean_ms: f64,
    growth_factor_mean: f64,
    growth_factor_max: f64,
    unsafe_supported_depth: usize,
}

fn main() -> AppResult<()> {
    let profiles = vec![
        Profile {
            name: "d=4, b=2^16",
            base: 1u64 << 16,
            digits: 4,
            // > 2*d*(b-1)^2 = 34,358,689,800
            bfv_plain_mod: 34_359_738_367,
            gadget_base: 256,
        },
        Profile {
            name: "d=8, b=2^8",
            base: 1u64 << 8,
            digits: 8,
            bfv_plain_mod: 1_040_407,
            gadget_base: 256,
        },
        Profile {
            name: "d=16, b=2^4",
            base: 1u64 << 4,
            digits: 16,
            // > 2*d*(b-1)^2 = 7,200
            bfv_plain_mod: 12_289,
            gadget_base: 16,
        },
    ];

    let mut results = Vec::new();
    for profile in profiles {
        let params = if profile.digits == 8 && profile.base == 256 && profile.bfv_plain_mod == 1_040_407 {
            u64_dbfv()?
        } else {
            make_profile_params(&profile)?
        };
        let measured = run_profile(profile.clone(), &params)?;
        results.push(measured);
    }

    print_table(&results);
    write_report(&results)?;
    Ok(())
}

fn make_profile_params(profile: &Profile) -> ExactoResult<Arc<DbfvParams>> {
    let bfv = BfvParamsBuilder::new()
        .ring_degree(4096)
        .plain_modulus(profile.bfv_plain_mod)
        .ct_moduli(vec![1_152_921_504_606_830_593])
        .aux_moduli(vec![18_014_398_509_998_081, 36_028_797_018_972_161])
        .gadget_base(profile.gadget_base)
        .sigma(3.2)
        .build()?;

    DbfvParams::new(bfv, profile.base, profile.digits, 0)
}

fn run_profile(profile: Profile, params: &Arc<DbfvParams>) -> ExactoResult<ProfileResult> {
    let mut rng = ChaCha20Rng::seed_from_u64(1337 + profile.digits as u64);
    let sk = gen_secret_key_with_rng(&params.bfv_params, &mut rng)?;
    let rlk = gen_relin_key_with_rng(&sk, &mut rng)?;

    let enc_mean_ms = time_encrypt(params, &sk, &mut rng)?;
    let add_mean_us = time_add(params, &sk, &mut rng)?;
    let mul_mean_ms = time_mul(params, &sk, &rlk, &mut rng)?;
    let (growth_factor_mean, growth_factor_max) = measure_growth_factor(params, &sk, &rlk, &mut rng)?;
    let unsafe_supported_depth = measure_unsafe_supported_depth(params, &sk, &rlk, &mut rng)?;

    Ok(ProfileResult {
        profile,
        mul_mean_ms,
        add_mean_us,
        enc_mean_ms,
        growth_factor_mean,
        growth_factor_max,
        unsafe_supported_depth,
    })
}

fn time_encrypt(params: &Arc<DbfvParams>, sk: &SecretKey, rng: &mut ChaCha20Rng) -> ExactoResult<f64> {
    let trials = 20usize;
    let mut total = Duration::ZERO;
    for _ in 0..trials {
        let m = rng.random::<u64>();
        let start = Instant::now();
        let _ = dbfv_encrypt_sk_with_rng(m, sk, params, rng)?;
        total += start.elapsed();
    }
    Ok((total.as_secs_f64() * 1_000.0) / trials as f64)
}

fn time_add(params: &Arc<DbfvParams>, sk: &SecretKey, rng: &mut ChaCha20Rng) -> ExactoResult<f64> {
    let trials = 50usize;
    let ct_a = dbfv_encrypt_sk_with_rng(123_456_789, sk, params, rng)?;
    let ct_b = dbfv_encrypt_sk_with_rng(987_654_321, sk, params, rng)?;
    let mut total = Duration::ZERO;
    for _ in 0..trials {
        let start = Instant::now();
        let _ = dbfv_add(&ct_a, &ct_b)?;
        total += start.elapsed();
    }
    Ok((total.as_secs_f64() * 1_000_000.0) / trials as f64)
}

fn time_mul(
    params: &Arc<DbfvParams>,
    sk: &SecretKey,
    rlk: &RelinKey,
    rng: &mut ChaCha20Rng,
) -> ExactoResult<f64> {
    let trials = 20usize;
    let ct_a = dbfv_encrypt_sk_with_rng(13_579, sk, params, rng)?;
    let ct_b = dbfv_encrypt_sk_with_rng(24_680, sk, params, rng)?;
    let mut total = Duration::ZERO;
    for _ in 0..trials {
        let mut lhs = ct_a.clone();
        let mut rhs = ct_b.clone();
        lhs.mul_depth = 0;
        rhs.mul_depth = 0;
        let start = Instant::now();
        let _ = dbfv_mul(&lhs, &rhs, rlk)?;
        total += start.elapsed();
    }
    Ok((total.as_secs_f64() * 1_000.0) / trials as f64)
}

fn measure_growth_factor(
    params: &Arc<DbfvParams>,
    sk: &SecretKey,
    rlk: &RelinKey,
    rng: &mut ChaCha20Rng,
) -> ExactoResult<(f64, f64)> {
    let trials = 8usize;
    let mut growths = Vec::with_capacity(trials);

    for _ in 0..trials {
        let a = rng.random::<u64>();
        let b = rng.random::<u64>();
        let ct_a = dbfv_encrypt_sk_with_rng(a, sk, params, rng)?;
        let ct_b = dbfv_encrypt_sk_with_rng(b, sk, params, rng)?;
        let noise_a = max_limb_noise(&ct_a, sk)?;
        let noise_b = max_limb_noise(&ct_b, sk)?;
        let noise_in = noise_a.max(noise_b).max(1);

        if let Ok(ct_mul) = dbfv_mul(&ct_a, &ct_b, rlk) {
            if let Ok(noise_out) = max_limb_noise(&ct_mul, sk) {
                growths.push(noise_out as f64 / noise_in as f64);
            }
        }
    }

    if growths.is_empty() {
        return Err(ExactoError::DecryptionError);
    }

    let mean = growths.iter().sum::<f64>() / growths.len() as f64;
    let max = growths
        .iter()
        .copied()
        .fold(0.0_f64, |acc, x| if x > acc { x } else { acc });
    Ok((mean, max))
}

fn measure_unsafe_supported_depth(
    params: &Arc<DbfvParams>,
    sk: &SecretKey,
    rlk: &RelinKey,
    rng: &mut ChaCha20Rng,
) -> ExactoResult<usize> {
    let factor = 3u64;
    let mut ct_acc = dbfv_encrypt_sk_with_rng(5, sk, params, rng)?;
    let ct_factor = dbfv_encrypt_sk_with_rng(factor, sk, params, rng)?;
    let mut expected = 5u64;
    let mut depth = 0usize;

    for d in 1..=32 {
        // Benchmark-only mode: bypass API depth guard to probe empirical failure depth.
        let mut lhs = ct_acc.clone();
        let mut rhs = ct_factor.clone();
        lhs.mul_depth = 0;
        rhs.mul_depth = 0;

        let ct_next = match dbfv_mul(&lhs, &rhs, rlk) {
            Ok(ct) => ct,
            Err(_) => break,
        };

        expected = expected.wrapping_mul(factor);
        let got = dbfv_decrypt(&ct_next, sk)?;
        if got != expected {
            break;
        }
        depth = d;
        ct_acc = ct_next;
    }
    Ok(depth)
}

fn max_limb_noise(ct: &exacto::dbfv::DbfvCiphertext, sk: &SecretKey) -> ExactoResult<u64> {
    let mut max_noise = 0u64;
    for limb in &ct.limbs {
        let limb_noise = bfv_noise_inf(limb, sk)?;
        if limb_noise > max_noise {
            max_noise = limb_noise;
        }
    }
    Ok(max_noise)
}

fn bfv_noise_inf(ct: &exacto::bfv::BfvCiphertext, sk: &SecretKey) -> ExactoResult<u64> {
    let params = &ct.params;
    let q = params.ct_basis.moduli[0];
    let t = params.plain_modulus;
    let delta = q / t;

    let pt = bfv_decrypt(ct, sk)?;
    let mut phase = ct.c[0].clone();
    let mut s_power = sk.poly.clone();
    for i in 1..ct.c.len() {
        let c_i_s = ct.c[i].mul(&s_power)?;
        phase = phase.add(&c_i_s)?;
        if i < ct.c.len() - 1 {
            s_power = s_power.mul(&sk.poly)?;
        }
    }

    let phase_coeffs = phase.to_coeff_poly(&params.ct_basis);
    let mut max_abs = 0u64;
    for i in 0..params.ring_degree {
        let target = ((pt.coeffs[i] as u128 * delta as u128) % q as u128) as u64;
        let raw = if phase_coeffs.coeffs[i] >= target {
            phase_coeffs.coeffs[i] - target
        } else {
            phase_coeffs.coeffs[i] + q - target
        };
        let centered = raw.min(q - raw);
        if centered > max_abs {
            max_abs = centered;
        }
    }
    Ok(max_abs)
}

fn print_table(results: &[ProfileResult]) {
    println!("Paper-style Reproduction (p = 2^64 scalar mode, n = 4096)");
    println!(
        "{:<16} {:>12} {:>12} {:>12} {:>16} {:>16} {:>12}",
        "profile",
        "enc_ms",
        "add_us",
        "mul_ms",
        "growth_mean",
        "growth_max",
        "depth*"
    );
    for r in results {
        println!(
            "{:<16} {:>12.3} {:>12.3} {:>12.3} {:>16.3} {:>16.3} {:>12}",
            r.profile.name,
            r.enc_mean_ms,
            r.add_mean_us,
            r.mul_mean_ms,
            r.growth_factor_mean,
            r.growth_factor_max,
            r.unsafe_supported_depth
        );
    }
    println!("* depth is an UNSAFE benchmark-only estimate (guard bypass), not API-supported depth.");
}

fn write_report(results: &[ProfileResult]) -> AppResult<()> {
    let mut dir = PathBuf::from("reports");
    fs::create_dir_all(&dir)?;
    dir.push("paper_reproduction.md");
    let mut f = File::create(&dir)?;

    writeln!(f, "# Paper Reproduction Report")?;
    writeln!(f)?;
    writeln!(f, "Command: `cargo run --release --bin paper_repro`")?;
    writeln!(f, "Scope: p=2^64 scalar dBFV mode, n=4096, CPU-only local run.")?;
    writeln!(f)?;
    writeln!(f, "| profile | enc (ms) | add (us) | mul (ms) | growth mean | growth max | unsafe depth* |")?;
    writeln!(f, "|---|---:|---:|---:|---:|---:|---:|")?;
    for r in results {
        writeln!(
            f,
            "| {} | {:.3} | {:.3} | {:.3} | {:.3} | {:.3} | {} |",
            r.profile.name,
            r.enc_mean_ms,
            r.add_mean_us,
            r.mul_mean_ms,
            r.growth_factor_mean,
            r.growth_factor_max,
            r.unsafe_supported_depth
        )?;
    }
    writeln!(f)?;
    writeln!(
        f,
        "* `unsafe depth` is measured by bypassing the public `mul_depth` guard for experimental comparison only."
    )?;
    writeln!(
        f,
        "  Use `dbfv_mul_then_bootstrap` for supported chaining in application code."
    )?;
    Ok(())
}
