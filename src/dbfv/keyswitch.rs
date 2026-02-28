use crate::error::Result;
use crate::bfv::keygen::RelinKey;
use crate::bfv::keyswitch::relinearize;
use crate::dbfv::ciphertext::DbfvCiphertext;

/// Per-limb key switching for dBFV ciphertexts.
///
/// Applies BFV relinearization to each limb independently.
pub fn dbfv_relinearize(
    ct: &DbfvCiphertext,
    rlk: &RelinKey,
) -> Result<DbfvCiphertext> {
    let limbs = ct.limbs.iter()
        .map(|limb| relinearize(limb, rlk))
        .collect::<Result<Vec<_>>>()?;

    Ok(DbfvCiphertext {
        limbs,
        degree: ct.degree,
        mul_depth: ct.mul_depth,
        params: ct.params.clone(),
    })
}
