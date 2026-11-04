//! The kernel's round refusals ([`solidcore::blend::round_room::RoundRefusal`]) as [`Error`]s.

use solidcore::blend::round_room::RoundRefusal;
use solidcore::{Brep, SolidError};

use crate::error::status;
use crate::Error;

/// The kernel's round refusal in `e` as an [`Error`], its vertex as `b`'s record index and its
/// station measured from the edge's end when `reversed`; any other kernel error by its type.
pub(crate) fn refusal(b: &Brep, e: SolidError, reversed: bool) -> Error {
    let SolidError::RoundRefused { refusal, .. } = e else {
        return status(e);
    };
    match refusal {
        RoundRefusal::RoundsCross { at, fit } => Error::RoundsCross {
            at: at.iter().map(|p| [p.x, p.y, p.z]).collect(),
            fit,
        },
        RoundRefusal::TaperEndUnsupported { vertex, lean } => {
            match b.vertex_list().iter().position(|&v| v == vertex) {
                Some(vertex) => Error::TaperEndUnsupported { vertex, lean },
                None => Error::Exception,
            }
        }
        RoundRefusal::FacesCloseAtBothEnds { .. } => Error::FacesCloseAtBothEnds,
        RoundRefusal::RoomRunsOut { at, largest, .. } => Error::RoomRunsOut {
            at: if reversed { 1.0 - at } else { at },
            largest,
        },
        _ => Error::Exception,
    }
}
