//! A procedural blend to its record and back (#2388): [`to_record`] writes the definition
//! and the net as plain data, [`from_record`] makes the blend from it exactly, and
//! [`rebuild`] also re-certifies the net against the definition
//! ([`ProceduralBlend::recertify`]), refusing a net that is off it. This is the record ↔
//! kernel conversion `docs/layers.md` gives `solidcore-ops`: it copies numbers and computes
//! nothing, and moves there when that crate exists.
//!
//! Not upstream. The fork carries a definition through its database only
//! (`SmProceduralBlend::WriteToDB`, [`ProceduralBlend::read_from_db`]).
//!
//! Deviations: none (no upstream counterpart). Recursion: nested supports are walked at
//! most [`BLEND_MAX_DEPTH`] deep, with an explicit depth refused past it.

use super::record::ProceduralBlendRecord;
use super::ProceduralBlend;
use crate::error::{Result, SolidError};
use crate::kernel::Kernel;

const OP_TO: &str = "surface::procedural_blend::convert::to_record";

/// The record of a procedural blend (not written yet: #2388).
///
/// # Errors
///
/// Always `Unsupported`.
pub fn to_record(_k: &Kernel, _pb: &ProceduralBlend) -> Result<ProceduralBlendRecord> {
    Err(SolidError::Unsupported {
        operation: OP_TO,
        reason: "not written yet",
    })
}

/// The blend a record holds (not written yet: #2388).
///
/// # Errors
///
/// Always `Unsupported`.
pub fn from_record(_k: &Kernel, _rec: &ProceduralBlendRecord) -> Result<ProceduralBlend> {
    Err(SolidError::Unsupported {
        operation: OP_TO,
        reason: "not written yet",
    })
}

/// [`from_record`] and the net re-certified (not written yet: #2388).
///
/// # Errors
///
/// Always `Unsupported`.
pub fn rebuild(k: &Kernel, rec: &ProceduralBlendRecord) -> Result<ProceduralBlend> {
    let pb = from_record(k, rec)?;
    pb.recertify(k)?;
    Ok(pb)
}
