//! Why an operation failed, as data. An application words these for its users; the
//! [`std::fmt::Display`] here is a plain statement of the variant for logs and tests.

use solidcore::SolidError;

/// SMLib's `SM_ERR_INVALID_INPUT`, the status of input the kernel refuses.
pub const INVALID_INPUT: i32 = 1007;

/// Why an operation failed.
#[derive(Clone, Debug, PartialEq)]
pub enum Error {
    /// The kernel's own status (1001 `SM_ERR`, 1007 `SM_ERR_INVALID_INPUT`, …).
    Status(i32),
    /// An argument refused before the kernel saw it: not finite, out of range, or an index
    /// past its list.
    Argument,
    /// The kernel reported success but the result has no faces (for a boolean: it is empty).
    Empty,
    /// The kernel reported success, but a solid input gave a result that is not a closed solid.
    NotSolid,
    /// No room for the operation: an inward offset (a shell's wall, faces moved in) leaves none
    /// inside the body, or a blend is larger than the faces it rolls on allow.
    NoRoom,
    /// An offset whose walls would meet along a line (a shell's): the offsets of two faces
    /// touch without crossing.
    WallsTouch,
    /// A blend larger than its faces allow, refused with the largest size that builds it
    /// (the kernel's fit, `SolidError::FilletRefused`): the size, what it sizes, a point at the
    /// middle of the edge it was found for, and why there is no room, when the refusal names them.
    TooLarge {
        /// The largest size that builds it.
        size: f64,
        /// What the size sizes.
        kind: FitKind,
        /// A point at the middle of the edge it was found for.
        at: Option<[f64; 3]>,
        /// Why there is no room.
        why: Option<RoomReason>,
    },
    /// A shell, offset or thicken whose wall is thicker than a face curves (the kernel's
    /// `ShellRefusal::Fold`): the largest wall that works, when the kernel names one.
    WallFolds {
        /// The largest wall that works.
        largest: Option<f64>,
    },
    /// The kernel cannot build this exactly yet (a variable radius whose law has a corner).
    Unsupported,
    /// A blend's edges meet at a vertex of a kind no blend corner handles yet.
    CornerUnsupported,
    /// The kernel does not set back corner `corner` (an index into the call's corners): it sets
    /// back a corner only where every edge meeting it is blended, all convex or all concave.
    CornerRefused {
        /// The corner's index among the call's corners.
        corner: usize,
    },
    /// The kernel reported success but the operation did nothing: a blend that made no face, a
    /// face move or draft that moved no face.
    NoEffect,
    /// The answer is not exact in the input's geometry.
    NotExact,
    /// The kernel found its own state inconsistent part way through a call.
    Exception,
    /// A record the conversion cannot carry to the kernel, and why.
    Record(String),
    /// The kernel reported an error: the call, then its message.
    Kernel(String),
    /// The kernel panicked; the session caught it: the call, then the panic's message.
    Internal(String),
    /// A step the kernel has not ported yet, refused rather than built another way: what is
    /// missing.
    Unported(&'static str),
}

/// Why a blend has no room, as the kernel names its refusal (`NoRoomReason`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RoomReason {
    /// The rolling ball does not fit a face it rolls on.
    BallDoesNotFit,
    /// A contact line leaves the face it rolls on.
    RailLeavesFace,
    /// The blends take the whole of a face they roll on.
    FaceConsumed,
    /// One blend's contact line runs inside the band another takes off a face they share.
    BandsOverlap,
    /// The balls fold over one another along the edge.
    BallsFold,
    /// A face curves tighter than a chamfer's offset along it.
    OffsetFolds,
    /// The blend would take all of the part between its contact lines and cut it in pieces.
    CutsThePart,
    /// The blend's section would break through a wall into the hollow behind it.
    BreaksThrough,
    /// The ball is as wide as a face it rolls on, so its centre stays put.
    OffsetCollapses,
    /// The blend's corners take the whole of its edge.
    CornersTakeWholeEdge,
    /// A blend past its faces runs off a face into the material.
    RunsIntoMaterial,
}

/// What the largest size that builds a refused blend ([`Error::TooLarge`]) sizes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FitKind {
    /// A round's radius.
    Radius,
    /// An equal chamfer's distance.
    Distance,
    /// A factor on a chamfer's two distances.
    Scale,
    /// A chord round's width.
    Chord,
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::Status(INVALID_INPUT) => write!(f, "solidcore rejected the input (SM_ERR_INVALID_INPUT)"),
            Error::Status(s) => write!(f, "solidcore failed (status {s})"),
            Error::Argument => write!(f, "invalid argument"),
            Error::Empty => write!(f, "the result is empty"),
            Error::NotSolid => write!(f, "the result is not a closed solid"),
            Error::NoRoom => write!(f, "no room"),
            Error::WallsTouch => write!(f, "the offset's walls meet"),
            Error::TooLarge { size, kind, .. } => write!(f, "too large: the largest {kind:?} that fits is {size}"),
            Error::WallFolds { largest } => write!(f, "the wall folds (largest {largest:?})"),
            Error::Unsupported => write!(f, "not supported yet"),
            Error::CornerUnsupported => write!(f, "a blend corner of this kind is not supported yet"),
            Error::CornerRefused { corner } => write!(f, "corner {corner} is not set back"),
            Error::NoEffect => write!(f, "the operation did nothing"),
            Error::NotExact => write!(f, "not exact in the input's geometry"),
            Error::Exception => write!(f, "solidcore raised an internal exception"),
            Error::Record(w) => write!(f, "record: {w}"),
            Error::Kernel(e) => write!(f, "solidcore: {e}"),
            Error::Internal(e) => write!(f, "solidcore failed internally in {e}"),
            Error::Unported(what) => write!(f, "not supported yet: {what}"),
        }
    }
}

impl std::error::Error for Error {}

/// An operation's result.
pub type Result<T> = std::result::Result<T, Error>;

/// A kernel error as the operation's: the kernel's status where it stands for one, a
/// port gap as [`Error::Unsupported`], any other as [`Error::Kernel`].
#[must_use]
pub fn status(e: SolidError) -> Error {
    match e {
        SolidError::Upstream { code, .. } => Error::Status(code),
        SolidError::InvalidInput { .. } => Error::Status(INVALID_INPUT),
        e if e.is_unported() => Error::Unsupported,
        e => Error::Kernel(e.to_string()),
    }
}
