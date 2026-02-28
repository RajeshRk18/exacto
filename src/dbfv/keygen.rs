use std::sync::Arc;
use crate::error::Result;
use crate::params::DbfvParams;
use crate::bfv::keygen::{
    self, SecretKey, PublicKey, RelinKey, GaloisKey,
};

/// Generate keys for dBFV (wraps BFV keygen with dBFV parameters).
pub fn dbfv_keygen(params: &Arc<DbfvParams>) -> Result<(SecretKey, PublicKey, RelinKey)> {
    let sk = keygen::gen_secret_key(&params.bfv_params)?;
    let pk = keygen::gen_public_key(&sk)?;
    let rlk = keygen::gen_relin_key(&sk)?;
    Ok((sk, pk, rlk))
}

/// Generate a full key set including Galois keys for rotations.
pub fn dbfv_keygen_full(
    params: &Arc<DbfvParams>,
    galois_elements: &[usize],
) -> Result<(SecretKey, PublicKey, RelinKey, Vec<GaloisKey>)> {
    let sk = keygen::gen_secret_key(&params.bfv_params)?;
    let pk = keygen::gen_public_key(&sk)?;
    let rlk = keygen::gen_relin_key(&sk)?;

    let gks: Vec<GaloisKey> = galois_elements.iter()
        .map(|&elem| keygen::gen_galois_key(&sk, elem))
        .collect::<Result<Vec<_>>>()?;

    Ok((sk, pk, rlk, gks))
}
