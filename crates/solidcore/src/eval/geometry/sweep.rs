//! Procedural sweep definitions and their display nets as plain data.
//!
//! These kernel DTOs mirror the portable record without depending on that layer.
//! Reconstruction and certification belong to the procedural sweep implementation.

/// A definition curve, with its original polynomial or rational coefficients.
#[derive(Clone, Debug, PartialEq)]
pub struct NurbsCurveRecord {
    /// Image dimension.
    pub dim: u32,
    /// Whether the stored coefficients include z.
    pub points_3d: bool,
    /// Polynomial degree.
    pub degree: u32,
    /// Full knot vector.
    pub knots: Vec<f64>,
    /// Weighted coordinates when rational, unweighted coordinates otherwise.
    pub points: Vec<[f64; 3]>,
    /// Rational weights, absent for a polynomial curve.
    pub weights: Option<Vec<f64>>,
}

/// The fitted surface net, with its original coefficients.
#[derive(Clone, Debug, PartialEq)]
pub struct NurbsSurfaceRecord {
    /// Degree along u.
    pub u_degree: u32,
    /// Degree along v.
    pub v_degree: u32,
    /// Control-point count along u.
    pub u_count: u32,
    /// Control-point count along v.
    pub v_count: u32,
    /// Full u knot vector.
    pub u_knots: Vec<f64>,
    /// Full v knot vector.
    pub v_knots: Vec<f64>,
    /// Coefficients in u-major order, weighted when rational.
    pub points: Vec<[f64; 3]>,
    /// Rational weights, absent for a polynomial surface.
    pub weights: Option<Vec<f64>>,
}

/// A piecewise linear radial displacement, held constant beyond its endpoints.
#[derive(Clone, Debug, PartialEq)]
pub struct RadialLawRecord {
    /// Unit displacement direction.
    pub direction: [f64; 3],
    /// Strictly ascending angles in radians.
    pub angles: Vec<f64>,
    /// Displacement at each angle.
    pub offsets: Vec<f64>,
}

/// A screw motion about an oriented axis.
#[derive(Clone, Debug, PartialEq)]
pub struct ScrewMotionRecord {
    /// A point on the axis.
    pub origin: [f64; 3],
    /// Unit axis direction.
    pub axis: [f64; 3],
    /// Signed rise per full turn.
    pub pitch: f64,
    /// Optional radial displacement law.
    pub law: Option<RadialLawRecord>,
}

/// The exact source of an approximated normal-offset section.
#[derive(Clone, Debug, PartialEq)]
pub struct ScrewOffsetRecord {
    /// Original affine profile, before the normal offset.
    pub source: NurbsCurveRecord,
    /// Signed normal-offset distance.
    pub signed_distance: f64,
    /// Certified section approximation error.
    pub approximation_bound: f64,
}

/// The construction defining a sweep; each motion has its own data.
#[derive(Clone, Debug, PartialEq)]
pub enum SweepDefinitionRecord {
    /// A profile carried by a screw motion.
    Screw {
        /// Profile used to construct the fitted net.
        profile: NurbsCurveRecord,
        /// Screw motion.
        motion: ScrewMotionRecord,
        /// Exact offset source, when the profile approximates an offset section.
        offset: Option<ScrewOffsetRecord>,
    },
}

/// A sweep's definition, fitted net and certified relationship, in the face's parameters.
#[derive(Clone, Debug, PartialEq)]
pub struct ProceduralSweepRecord {
    /// Exact defining construction.
    pub definition: SweepDefinitionRecord,
    /// Fitted display and intersection net.
    pub net: NurbsSurfaceRecord,
    /// Net stations along the turn, in radians.
    pub stations: Vec<f64>,
    /// Half angle of the net's spans.
    pub half: f64,
    /// Certified net-to-definition distance; the kernel verifies it on reconstruction.
    pub bound: f64,
    /// Affine parameter-map origins.
    pub map_a: [f64; 2],
    /// Affine parameter-map scales.
    pub map_b: [f64; 2],
    /// Whether the face's u and v exchange roles in the map.
    pub swapped: bool,
}
