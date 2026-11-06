//! Writing ISO 10303-21: a DATA section built entity by entity ([`Writer`]), and the Part 42
//! curves and surfaces of [`super::geom`] as ISO 10303-42 entities ([`GeomWriter`]). Identical
//! entities are written once and shared, as the reader shares them; an entity whose identity
//! matters is written on its own ([`Writer::add_distinct`]).
//!
//! The writer writes what the record holds, and computes nothing: the record keeps a file's
//! values as the file gives them (directions not normalised, a placement's axes optional and
//! not made orthogonal, a `VECTOR` as its orientation and magnitude), so writing them needs no
//! normalising, no axis chosen and no length taken. The only work is formatting: a real in the
//! shortest form that reads back to the same `f64` ([`real`]), a string with Part 21's escapes
//! ([`string`]), and entity ids counted from `#1` in the order entities are added. The record's
//! lengths are millimetres and its angles radians, so the representation context a caller writes
//! around this geometry must say so; that context, the topology and the product structure are
//! not written here (#2579: the topology writer and the AP242 product structure come with the
//! kernel's and `solidcore-ops`' parts of the export).
//!
//! Deterministic: the same calls give the same text, byte for byte. Shared entities are found
//! by their text in a `BTreeMap`, which is only looked up, never iterated; the text is written
//! from a `Vec` in insertion order. Every refusal is an error message, never a panic: a value
//! Part 21 cannot write (a non-finite real, a malformed keyword), geometry nested deeper than the
//! reader reads ([`MAX_GEOM_DEPTH`]), and data that cannot be written in millimetres and radians
//! without a computation (below).
//!
//! Origin: SolidCraft's `crates/step/src/write.rs` at `a5b71ca2` (MIT OR Apache-2.0), without
//! its `topology` module (a schema 4 record's solids: brep to STEP geometry with transformations,
//! a computation, which goes to the kernel and `solidcore-ops`).
//!
//! | Here | SolidCraft (`crates/step/src/write.rs`) |
//! |---|---|
//! | [`real`] | `real` (`None` for a non-finite real) |
//! | [`string`] | `string` (control characters escaped, not dropped) |
//! | [`Writer`] | `Writer` (fallible; `BTreeMap` instead of `HashMap`) |
//! | [`GeomWriter`] | `GeomWriter` (fallible) |
//! | [`GeomWriter::point`], [`GeomWriter::direction`], [`GeomWriter::placement`], [`GeomWriter::axis1`] | the same, with the record's optional axes |
//! | [`GeomWriter::curve`], [`GeomWriter::surface`] | the same, with the record's [`Vector`] |
//! | [`GeomWriter::pcurve`] | (new) the `PCURVE` branch of `curve_in` |
//!
//! Deviations from SolidCraft, none of them a computation:
//! - Every write returns a result (was: an id, whatever the data). Refused: a non-finite real
//!   (was: written as `0.`), an entity, enumeration or typed-value name that is not a Part 21
//!   keyword, a binary that is not hexadecimal, a count that does not fit an `INTEGER`, a 2D
//!   point or direction with a z component other than 0 or a 2D placement with an axis (was:
//!   z dropped), a curve that has no 2D form in a parameter space (`SURFACE_CURVE` and its
//!   subtypes, `OFFSET_CURVE_3D`), and a [`Vector`] whose unit factors are not all 1 (below).
//! - A placement writes `$` for a missing axis or reference direction, an `AXIS1_PLACEMENT` for
//!   a missing direction; a [`Vector`] writes its orientation and magnitude as the record holds
//!   them (was: a unit direction and the vector's length). Neither is normalised, because the
//!   record keeps the file's ratios.
//! - A [`Vector`] read in other units than millimetres and radians keeps the factors of its
//!   space (`scale`), and its components in millimetres are `scale[i] · magnitude ·
//!   orientation[i] / |orientation|`. Writing it in millimetres and radians would take that
//!   product, and a norm where the factors differ, so the writer refuses it: converting it is
//!   the kernel's (#2578, slice 4). Geometry an export takes from the kernel comes back to the
//!   record in millimetres, its factors 1.
//! - [`string`] writes an ASCII control character (`U+0000` to `U+001F`, and `U+007F`) as a
//!   `\X2\00hh\X0\` run (was: dropped), so every string reads back as it was. Not as `\X\hh`:
//!   the public descriptions of Part 21's string encoding give `\X\hh` as an 8-bit code of
//!   ISO 8859-1 and print the graphic characters (`U+0020` to `U+007E`) as themselves, and none
//!   we may read settles whether `\X\` covers the control codes below `U+0020` (the standard's
//!   own text is not freely available). A `\X2\` run is a sequence of ISO 10646 code units,
//!   which covers every code, so that form is valid whichever way `\X\` is read.
//! - Nesting is bounded as the reader bounds it: a curve or surface more than
//!   [`MAX_GEOM_DEPTH`] links below the entity asked for (a trim's or an offset's basis, a
//!   surface curve's 3D curve or associated geometry, a swept surface's profile, an offset or
//!   trimmed surface's basis; a pcurve's surface and 2D curve at the pcurve's own depth) is
//!   refused with the reader's message (was: unbounded recursion, a stack overflow for a deep
//!   enough record). So everything the reader reads in one call is written, and an in-memory
//!   record deeper than that is an error, not a crash (on wasm32, which has no stack guard,
//!   silent corruption).
//! - [`Writer::file`] leaves `FILE_NAME`'s `time_stamp` empty (as SolidCraft did), so the file
//!   is the same bytes for the same entities: a library that must be deterministic and build
//!   for wasm32 reads no clock. A caller that needs a stamp (a strict validator) is the export
//!   that writes the header for a product, which comes with the AP242 product structure (#2579).
//! - Shared entities are found in a `BTreeMap` (was a `HashMap`, seeded at random; only looked
//!   up there too). Ids are counted with checked arithmetic.
//! - [`Writer::file`] names `solidcore` as the originating system (was `SolidCraft`).

use std::collections::BTreeMap;
use std::fmt::Write as _;

use super::geom::{
    Associated, Axis1, BSplineCurve, BSplineSurface, Curve, CurveForm, KnotSpec, Logical, P3,
    Pcurve, Placement, PreferredRep, Surface, SurfaceCurveKind, SurfaceForm, Trim, TrimPreference,
    Vector,
};
use super::p21::Param;
use super::read::MAX_GEOM_DEPTH;

/// A result whose error is a message.
pub type R<T> = Result<T, String>;

/// A real as Part 21 writes it: always with a decimal point, the exponent as `E`, in the
/// shortest form that reads back to the same `f64` (Rust's `{:?}`), so writing and reading a
/// real loses nothing: `1.0`, `0.1`, `1.E-7`, `-0.0`. `None` for a NaN or an infinity, which
/// Part 21 cannot write.
pub fn real(x: f64) -> Option<String> {
    if !x.is_finite() {
        return None;
    }
    let s = format!("{x:?}");
    Some(match s.split_once('e') {
        Some((m, e)) if m.contains('.') => format!("{m}E{e}"),
        Some((m, e)) => format!("{m}.E{e}"),
        None if s.contains('.') => s,
        None => format!("{s}."),
    })
}

/// A string literal: the apostrophe doubled, the backslash doubled, and every other character
/// that is not printable ASCII (a control character, or one outside ASCII) as a `\X2\…\X0\` run
/// of its UTF-16 code units (a surrogate pair for one beyond the basic plane), so the text is
/// plain ASCII and reads back as it was. A control character is `\X2\000A\X0\`, not `\X\0A`:
/// the module's deviations say why.
pub fn string(s: &str) -> String {
    let mut o = String::from("'");
    for c in s.chars() {
        match c {
            '\'' => o.push_str("''"),
            '\\' => o.push_str("\\\\"),
            c if c.is_ascii() && !c.is_ascii_control() => o.push(c),
            c => {
                o.push_str("\\X2\\");
                let mut buf = [0u16; 2];
                for u in c.encode_utf16(&mut buf) {
                    let _ = write!(o, "{u:04X}");
                }
                o.push_str("\\X0\\");
            }
        }
    }
    o.push('\'');
    o
}

/// Whether `s` is a Part 21 keyword: a letter, then letters, digits and `_`, upper case (an
/// entity name, an enumeration value or a typed value's name).
fn keyword(s: &str) -> R<&str> {
    let mut chars = s.chars();
    let first = chars.next().is_some_and(|c| c.is_ascii_uppercase());
    if first && chars.all(|c| c.is_ascii_uppercase() || c.is_ascii_digit() || c == '_') {
        Ok(s)
    } else {
        Err(format!("`{s}` is not a Part 21 keyword"))
    }
}

fn param(p: &Param, o: &mut String) -> R<()> {
    match p {
        Param::Int(i) => {
            let _ = write!(o, "{i}");
        }
        Param::Real(x) => o.push_str(&real(*x).ok_or_else(|| format!("{x} is not a real"))?),
        Param::Str(s) => o.push_str(&string(s)),
        Param::Enum(e) => {
            let _ = write!(o, ".{}.", keyword(e)?);
        }
        Param::Ref(r) => {
            let _ = write!(o, "#{r}");
        }
        Param::List(v) => {
            o.push('(');
            params(v, o)?;
            o.push(')');
        }
        Param::Typed(name, v) => {
            o.push_str(keyword(name)?);
            o.push('(');
            params(v, o)?;
            o.push(')');
        }
        Param::Binary(b) => {
            if !b.chars().all(|c| c.is_ascii_hexdigit()) {
                return Err(format!("binary `{b}` is not hexadecimal"));
            }
            let _ = write!(o, "\"{b}\"");
        }
        Param::Null => o.push('$'),
        Param::Derived => o.push('*'),
    }
    Ok(())
}

fn params(v: &[Param], o: &mut String) -> R<()> {
    for (i, p) in v.iter().enumerate() {
        if i > 0 {
            o.push(',');
        }
        param(p, o)?;
    }
    Ok(())
}

/// A DATA section under construction. Entities are numbered from `#1` in the order they are
/// added; adding an entity identical to an earlier shared one returns the earlier one's id.
#[derive(Debug, Default)]
pub struct Writer {
    /// The entity lines, `#id=…;`, in id order.
    lines: Vec<String>,
    /// The shared entities by their text (looked up only, never iterated).
    shared: BTreeMap<String, u64>,
}

impl Writer {
    /// An empty DATA section.
    pub fn new() -> Self {
        Self::default()
    }

    fn push(&mut self, body: String) -> R<u64> {
        if let Some(id) = self.shared.get(&body) {
            return Ok(*id);
        }
        let id = self.push_distinct(&body)?;
        self.shared.insert(body, id);
        Ok(id)
    }

    fn push_distinct(&mut self, body: &str) -> R<u64> {
        let id = u64::try_from(self.lines.len())
            .ok()
            .and_then(|n| n.checked_add(1))
            .ok_or("too many entities")?;
        self.lines.push(format!("#{id}={body};"));
        Ok(id)
    }

    fn simple(name: &str, p: &[Param]) -> R<String> {
        let mut body = format!("{}(", keyword(name)?);
        params(p, &mut body)?;
        body.push(')');
        Ok(body)
    }

    /// A simple entity instance, shared with an identical earlier one.
    pub fn add(&mut self, name: &str, p: &[Param]) -> R<u64> {
        let body = Self::simple(name, p)?;
        self.push(body)
    }

    /// A simple entity instance of its own, even when an identical one exists: an entity whose
    /// identity matters (two vertices at one point are two vertices; a product is one product).
    pub fn add_distinct(&mut self, name: &str, p: &[Param]) -> R<u64> {
        let body = Self::simple(name, p)?;
        self.push_distinct(&body)
    }

    /// A complex entity instance, shared with an identical earlier one; its records are written
    /// in alphabetical order of their names, as Part 21 asks.
    pub fn add_complex(&mut self, records: &[(&str, Vec<Param>)]) -> R<u64> {
        let mut sorted: Vec<&(&str, Vec<Param>)> = records.iter().collect();
        sorted.sort_by(|a, b| a.0.cmp(b.0));
        let mut body = String::from("(");
        for (name, p) in sorted {
            body.push_str(&Self::simple(name, p)?);
        }
        body.push(')');
        self.push(body)
    }

    /// The number of entities added so far.
    pub fn len(&self) -> usize {
        self.lines.len()
    }

    /// Whether no entity has been added.
    pub fn is_empty(&self) -> bool {
        self.lines.is_empty()
    }

    /// The DATA section's entity lines, one per line, in id order.
    pub fn data(&self) -> String {
        let mut o = String::new();
        for l in &self.lines {
            o.push_str(l);
            o.push('\n');
        }
        o
    }

    /// A whole exchange file around the entities written so far: the HEADER with the file's
    /// `name` and `schema`, and the DATA section. `FILE_NAME`'s `time_stamp` is empty, so the
    /// text depends on the entities alone (the module's deviations).
    pub fn file(&self, name: &str, schema: &str) -> String {
        format!(
            "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION((''),'2;1');\nFILE_NAME({},'',(''),(''),'solidcore','solidcore','');\nFILE_SCHEMA(({}));\nENDSEC;\nDATA;\n{}ENDSEC;\nEND-ISO-10303-21;\n",
            string(name),
            string(schema),
            self.data()
        )
    }
}

/// The dimension geometry is written in: model space, or a surface's parameter space (and a
/// 2D placement).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Dim {
    Two,
    Three,
}

fn s() -> Param {
    Param::Str(String::new())
}
fn r(id: u64) -> Param {
    Param::Ref(id)
}
fn f(x: f64) -> Param {
    Param::Real(x)
}
fn e(name: &str) -> Param {
    Param::Enum(name.into())
}
fn list(v: Vec<Param>) -> Param {
    Param::List(v)
}
fn reals(v: &[f64]) -> Param {
    list(v.iter().map(|x| f(*x)).collect())
}
fn count(n: usize) -> R<Param> {
    i64::try_from(n)
        .map(Param::Int)
        .map_err(|_| format!("{n} does not fit a Part 21 integer"))
}
fn counts(v: &[usize]) -> R<Param> {
    Ok(list(v.iter().map(|x| count(*x)).collect::<R<_>>()?))
}
fn logical(l: Logical) -> Param {
    e(match l {
        Logical::True => "T",
        Logical::False => "F",
        Logical::Unknown => "U",
    })
}
fn boolean(b: bool) -> Param {
    e(if b { "T" } else { "F" })
}

fn curve_form(f: CurveForm) -> Param {
    e(match f {
        CurveForm::PolylineForm => "POLYLINE_FORM",
        CurveForm::CircularArc => "CIRCULAR_ARC",
        CurveForm::EllipticArc => "ELLIPTIC_ARC",
        CurveForm::ParabolicArc => "PARABOLIC_ARC",
        CurveForm::HyperbolicArc => "HYPERBOLIC_ARC",
        CurveForm::Unspecified => "UNSPECIFIED",
    })
}

fn surface_form(f: SurfaceForm) -> Param {
    e(match f {
        SurfaceForm::PlaneSurf => "PLANE_SURF",
        SurfaceForm::CylindricalSurf => "CYLINDRICAL_SURF",
        SurfaceForm::ConicalSurf => "CONICAL_SURF",
        SurfaceForm::SphericalSurf => "SPHERICAL_SURF",
        SurfaceForm::ToroidalSurf => "TOROIDAL_SURF",
        SurfaceForm::SurfOfRevolution => "SURF_OF_REVOLUTION",
        SurfaceForm::RuledSurf => "RULED_SURF",
        SurfaceForm::GeneralisedCone => "GENERALISED_CONE",
        SurfaceForm::QuadricSurf => "QUADRIC_SURF",
        SurfaceForm::SurfOfLinearExtrusion => "SURF_OF_LINEAR_EXTRUSION",
        SurfaceForm::Unspecified => "UNSPECIFIED",
    })
}

fn knot_spec(k: KnotSpec) -> Param {
    e(match k {
        KnotSpec::UniformKnots => "UNIFORM_KNOTS",
        KnotSpec::QuasiUniformKnots => "QUASI_UNIFORM_KNOTS",
        KnotSpec::PiecewiseBezierKnots => "PIECEWISE_BEZIER_KNOTS",
        KnotSpec::Unspecified => "UNSPECIFIED",
    })
}

/// The coordinates of `p` in `dim`: all three in model space, the first two in a parameter
/// space, whose z must be 0 (the record holds a 2D point with z = 0).
fn coords(p: P3, dim: Dim) -> R<Param> {
    match (dim, p) {
        (Dim::Three, p) => Ok(reals(&p)),
        (Dim::Two, [x, y, 0.0]) => Ok(reals(&[x, y])),
        (Dim::Two, _) => Err(format!("2D coordinates {p:?} with z other than 0")),
    }
}

/// Writes the record's geometry into a [`Writer`]: lengths in millimetres and angles in
/// radians, as the record holds them (the representation context the caller writes must say
/// so). Every curve and surface is written with the entities it refers to, each shared with an
/// identical one written before.
pub struct GeomWriter<'w> {
    /// The DATA section written into (the caller adds its topology and contexts there).
    pub w: &'w mut Writer,
    /// The 2D parameter-space context pcurves refer to, once written.
    pcurve_context: Option<u64>,
}

impl<'w> GeomWriter<'w> {
    /// A geometry writer adding to `w`.
    pub fn new(w: &'w mut Writer) -> Self {
        GeomWriter {
            w,
            pcurve_context: None,
        }
    }

    /// A `CARTESIAN_POINT` in model space.
    pub fn point(&mut self, p: P3) -> R<u64> {
        self.point_in(p, Dim::Three)
    }

    fn point_in(&mut self, p: P3, dim: Dim) -> R<u64> {
        let c = coords(p, dim)?;
        self.w.add("CARTESIAN_POINT", &[s(), c])
    }

    /// A `DIRECTION` in model space, its ratios as given (not normalised).
    pub fn direction(&mut self, d: P3) -> R<u64> {
        self.direction_in(d, Dim::Three)
    }

    fn direction_in(&mut self, d: P3, dim: Dim) -> R<u64> {
        let c = coords(d, dim)?;
        self.w.add("DIRECTION", &[s(), c])
    }

    /// A direction that may be absent: its reference, or `$`.
    fn optional_direction(&mut self, d: Option<P3>, dim: Dim) -> R<Param> {
        Ok(match d {
            Some(d) => r(self.direction_in(d, dim)?),
            None => Param::Null,
        })
    }

    /// A `VECTOR`: its orientation and magnitude as the record holds them. Refused when its
    /// unit factors are not all 1: the vector was read in other units than millimetres and
    /// radians, and writing it in those would be a computation (the module's deviations).
    fn vector_in(&mut self, v: &Vector, dim: Dim) -> R<u64> {
        if v.scale != [1.0; 3] {
            return Err(format!(
                "a VECTOR with unit factors {:?}: converting it to millimetres and radians is the kernel's",
                v.scale
            ));
        }
        let d = self.direction_in(v.orientation, dim)?;
        self.w.add("VECTOR", &[s(), r(d), f(v.magnitude)])
    }

    fn placement_in(&mut self, p: &Placement, dim: Dim) -> R<u64> {
        let o = self.point_in(p.origin, dim)?;
        match dim {
            Dim::Two => {
                if p.axis.is_some() {
                    return Err("a placement in a parameter space has no axis".into());
                }
                let x = self.optional_direction(p.ref_direction, dim)?;
                self.w.add("AXIS2_PLACEMENT_2D", &[s(), r(o), x])
            }
            Dim::Three => {
                let z = self.optional_direction(p.axis, dim)?;
                let x = self.optional_direction(p.ref_direction, dim)?;
                self.w.add("AXIS2_PLACEMENT_3D", &[s(), r(o), z, x])
            }
        }
    }

    /// An `AXIS2_PLACEMENT_3D`, with `$` for an axis or reference direction the record does not
    /// hold.
    pub fn placement(&mut self, p: &Placement) -> R<u64> {
        self.placement_in(p, Dim::Three)
    }

    /// An `AXIS1_PLACEMENT`, with `$` for a direction the record does not hold.
    pub fn axis1(&mut self, a: &Axis1) -> R<u64> {
        let o = self.point(a.origin)?;
        let d = self.optional_direction(a.direction, Dim::Three)?;
        self.w.add("AXIS1_PLACEMENT", &[s(), r(o), d])
    }

    /// A curve in model space.
    pub fn curve(&mut self, c: &Curve) -> R<u64> {
        self.curve_in(c, Dim::Three, 0)
    }

    /// A curve `depth` links below the entity the caller asked for, counted as the reader counts
    /// them ([`MAX_GEOM_DEPTH`]).
    fn curve_in(&mut self, c: &Curve, dim: Dim, depth: usize) -> R<u64> {
        if depth > MAX_GEOM_DEPTH {
            return Err("curve definitions nested too deeply".into());
        }
        match c {
            Curve::Line { origin, vector } => {
                let p = self.point_in(*origin, dim)?;
                let v = self.vector_in(vector, dim)?;
                self.w.add("LINE", &[s(), r(p), r(v)])
            }
            Curve::Circle { position, radius } => {
                let p = self.placement_in(position, dim)?;
                self.w.add("CIRCLE", &[s(), r(p), f(*radius)])
            }
            Curve::Ellipse {
                position,
                semi_axis_1,
                semi_axis_2,
            } => {
                let p = self.placement_in(position, dim)?;
                self.w
                    .add("ELLIPSE", &[s(), r(p), f(*semi_axis_1), f(*semi_axis_2)])
            }
            Curve::Hyperbola {
                position,
                semi_axis,
                semi_imag_axis,
            } => {
                let p = self.placement_in(position, dim)?;
                self.w
                    .add("HYPERBOLA", &[s(), r(p), f(*semi_axis), f(*semi_imag_axis)])
            }
            Curve::Parabola {
                position,
                focal_dist,
            } => {
                let p = self.placement_in(position, dim)?;
                self.w.add("PARABOLA", &[s(), r(p), f(*focal_dist)])
            }
            Curve::BSpline(b) => self.bspline_curve(b, dim),
            Curve::Polyline(pts) => {
                let ids = pts
                    .iter()
                    .map(|p| Ok(r(self.point_in(*p, dim)?)))
                    .collect::<R<_>>()?;
                self.w.add("POLYLINE", &[s(), list(ids)])
            }
            Curve::Trimmed {
                basis,
                trim_1,
                trim_2,
                sense,
                preference,
            } => {
                let b = self.curve_in(basis, dim, depth + 1)?;
                let t1 = self.trims(trim_1, dim)?;
                let t2 = self.trims(trim_2, dim)?;
                let pref = e(match preference {
                    TrimPreference::Cartesian => "CARTESIAN",
                    TrimPreference::Parameter => "PARAMETER",
                    TrimPreference::Unspecified => "UNSPECIFIED",
                });
                self.w
                    .add("TRIMMED_CURVE", &[s(), r(b), t1, t2, boolean(*sense), pref])
            }
            Curve::SurfaceCurve {
                kind,
                curve_3d,
                associated,
                master,
            } if dim == Dim::Three => {
                let c3 = self.curve_in(curve_3d, Dim::Three, depth + 1)?;
                let items = associated
                    .iter()
                    .map(|a| {
                        Ok(r(match a {
                            Associated::Pcurve(pc) => self.pcurve_in(pc, depth + 1)?,
                            Associated::Surface(sf) => self.surface_in(sf, depth + 1)?,
                        }))
                    })
                    .collect::<R<_>>()?;
                let name = match kind {
                    SurfaceCurveKind::Surface => "SURFACE_CURVE",
                    SurfaceCurveKind::Seam => "SEAM_CURVE",
                    SurfaceCurveKind::Intersection => "INTERSECTION_CURVE",
                    SurfaceCurveKind::BoundedSurface => "BOUNDED_SURFACE_CURVE",
                };
                let m = e(match master {
                    PreferredRep::Curve3d => "CURVE_3D",
                    PreferredRep::PcurveS1 => "PCURVE_S1",
                    PreferredRep::PcurveS2 => "PCURVE_S2",
                });
                self.w.add(name, &[s(), r(c3), list(items), m])
            }
            Curve::OffsetCurve3d {
                basis,
                distance,
                self_intersect,
                ref_direction,
            } if dim == Dim::Three => {
                let b = self.curve_in(basis, Dim::Three, depth + 1)?;
                let d = self.direction(*ref_direction)?;
                self.w.add(
                    "OFFSET_CURVE_3D",
                    &[s(), r(b), f(*distance), logical(*self_intersect), r(d)],
                )
            }
            Curve::SurfaceCurve { .. } | Curve::OffsetCurve3d { .. } => {
                Err("a surface curve or a 3D offset curve has no form in a parameter space".into())
            }
        }
    }

    /// A `PCURVE`: its surface, its 2D curve in the surface's parameters (the record's
    /// millimetres and radians), the `DEFINITIONAL_REPRESENTATION` holding the curve, and the
    /// one 2D parametric context every pcurve of this writer shares.
    pub fn pcurve(&mut self, pc: &Pcurve) -> R<u64> {
        self.pcurve_in(pc, 0)
    }

    /// A pcurve `depth` links down: its surface and its 2D curve at the same depth, as the
    /// reader reads them.
    fn pcurve_in(&mut self, pc: &Pcurve, depth: usize) -> R<u64> {
        let surf = self.surface_in(&pc.surface, depth)?;
        let c2 = self.curve_in(&pc.curve, Dim::Two, depth)?;
        let ctx = self.pcurve_context()?;
        let rep = self.w.add(
            "DEFINITIONAL_REPRESENTATION",
            &[s(), list(vec![r(c2)]), r(ctx)],
        )?;
        self.w.add("PCURVE", &[s(), r(surf), r(rep)])
    }

    fn pcurve_context(&mut self) -> R<u64> {
        if let Some(c) = self.pcurve_context {
            return Ok(c);
        }
        let c = self.w.add_complex(&[
            ("GEOMETRIC_REPRESENTATION_CONTEXT", vec![Param::Int(2)]),
            ("PARAMETRIC_REPRESENTATION_CONTEXT", vec![]),
            (
                "REPRESENTATION_CONTEXT",
                vec![Param::Str("2D SPACE".into()), s()],
            ),
        ])?;
        self.pcurve_context = Some(c);
        Ok(c)
    }

    fn trims(&mut self, t: &[Trim], dim: Dim) -> R<Param> {
        Ok(list(
            t.iter()
                .map(|t| match t {
                    Trim::Point(p) => Ok(r(self.point_in(*p, dim)?)),
                    Trim::Parameter(x) => Ok(Param::Typed("PARAMETER_VALUE".into(), vec![f(*x)])),
                })
                .collect::<R<_>>()?,
        ))
    }

    fn bspline_curve(&mut self, b: &BSplineCurve, dim: Dim) -> R<u64> {
        let ctrl = list(
            b.points
                .iter()
                .map(|p| Ok(r(self.point_in(*p, dim)?)))
                .collect::<R<_>>()?,
        );
        let common = vec![
            count(b.degree)?,
            ctrl,
            curve_form(b.form),
            logical(b.closed),
            logical(b.self_intersect),
        ];
        let knots = vec![
            counts(&b.multiplicities)?,
            reals(&b.knots),
            knot_spec(b.knot_spec),
        ];
        match &b.weights {
            None => {
                let mut p = vec![s()];
                p.extend(common);
                p.extend(knots);
                self.w.add("B_SPLINE_CURVE_WITH_KNOTS", &p)
            }
            Some(w) => self.w.add_complex(&[
                ("BOUNDED_CURVE", vec![]),
                ("B_SPLINE_CURVE", common),
                ("B_SPLINE_CURVE_WITH_KNOTS", knots),
                ("CURVE", vec![]),
                ("GEOMETRIC_REPRESENTATION_ITEM", vec![]),
                ("RATIONAL_B_SPLINE_CURVE", vec![reals(w)]),
                ("REPRESENTATION_ITEM", vec![s()]),
            ]),
        }
    }

    /// A surface.
    pub fn surface(&mut self, sf: &Surface) -> R<u64> {
        self.surface_in(sf, 0)
    }

    /// A surface `depth` links down ([`MAX_GEOM_DEPTH`]).
    fn surface_in(&mut self, sf: &Surface, depth: usize) -> R<u64> {
        if depth > MAX_GEOM_DEPTH {
            return Err("surface definitions nested too deeply".into());
        }
        match sf {
            Surface::Plane { position } => {
                let p = self.placement(position)?;
                self.w.add("PLANE", &[s(), r(p)])
            }
            Surface::Cylinder { position, radius } => {
                let p = self.placement(position)?;
                self.w.add("CYLINDRICAL_SURFACE", &[s(), r(p), f(*radius)])
            }
            Surface::Cone {
                position,
                radius,
                semi_angle,
            } => {
                let p = self.placement(position)?;
                self.w
                    .add("CONICAL_SURFACE", &[s(), r(p), f(*radius), f(*semi_angle)])
            }
            Surface::Sphere { position, radius } => {
                let p = self.placement(position)?;
                self.w.add("SPHERICAL_SURFACE", &[s(), r(p), f(*radius)])
            }
            Surface::Torus {
                position,
                major_radius,
                minor_radius,
            } => {
                let p = self.placement(position)?;
                self.w.add(
                    "TOROIDAL_SURFACE",
                    &[s(), r(p), f(*major_radius), f(*minor_radius)],
                )
            }
            Surface::DegenerateTorus {
                position,
                major_radius,
                minor_radius,
                select_outer,
            } => {
                let p = self.placement(position)?;
                self.w.add(
                    "DEGENERATE_TOROIDAL_SURFACE",
                    &[
                        s(),
                        r(p),
                        f(*major_radius),
                        f(*minor_radius),
                        boolean(*select_outer),
                    ],
                )
            }
            Surface::BSpline(b) => self.bspline_surface(b),
            Surface::Revolution { profile, axis } => {
                let c = self.curve_in(profile, Dim::Three, depth + 1)?;
                let a = self.axis1(axis)?;
                self.w.add("SURFACE_OF_REVOLUTION", &[s(), r(c), r(a)])
            }
            Surface::Extrusion { profile, vector } => {
                let c = self.curve_in(profile, Dim::Three, depth + 1)?;
                let v = self.vector_in(vector, Dim::Three)?;
                self.w
                    .add("SURFACE_OF_LINEAR_EXTRUSION", &[s(), r(c), r(v)])
            }
            Surface::Offset {
                basis,
                distance,
                self_intersect,
            } => {
                let b = self.surface_in(basis, depth + 1)?;
                self.w.add(
                    "OFFSET_SURFACE",
                    &[s(), r(b), f(*distance), logical(*self_intersect)],
                )
            }
            Surface::RectangularTrimmed {
                basis,
                u1,
                u2,
                v1,
                v2,
                u_sense,
                v_sense,
            } => {
                let b = self.surface_in(basis, depth + 1)?;
                self.w.add(
                    "RECTANGULAR_TRIMMED_SURFACE",
                    &[
                        s(),
                        r(b),
                        f(*u1),
                        f(*u2),
                        f(*v1),
                        f(*v2),
                        boolean(*u_sense),
                        boolean(*v_sense),
                    ],
                )
            }
        }
    }

    fn bspline_surface(&mut self, b: &BSplineSurface) -> R<u64> {
        let ctrl = list(
            b.points
                .iter()
                .map(|row| {
                    Ok(list(
                        row.iter()
                            .map(|p| Ok(r(self.point(*p)?)))
                            .collect::<R<_>>()?,
                    ))
                })
                .collect::<R<_>>()?,
        );
        let common = vec![
            count(b.u_degree)?,
            count(b.v_degree)?,
            ctrl,
            surface_form(b.form),
            logical(b.u_closed),
            logical(b.v_closed),
            logical(b.self_intersect),
        ];
        let knots = vec![
            counts(&b.u_multiplicities)?,
            counts(&b.v_multiplicities)?,
            reals(&b.u_knots),
            reals(&b.v_knots),
            knot_spec(b.knot_spec),
        ];
        match &b.weights {
            None => {
                let mut p = vec![s()];
                p.extend(common);
                p.extend(knots);
                self.w.add("B_SPLINE_SURFACE_WITH_KNOTS", &p)
            }
            Some(w) => self.w.add_complex(&[
                ("BOUNDED_SURFACE", vec![]),
                ("B_SPLINE_SURFACE", common),
                ("B_SPLINE_SURFACE_WITH_KNOTS", knots),
                ("GEOMETRIC_REPRESENTATION_ITEM", vec![]),
                (
                    "RATIONAL_B_SPLINE_SURFACE",
                    vec![list(w.iter().map(|row| reals(row)).collect())],
                ),
                ("REPRESENTATION_ITEM", vec![s()]),
                ("SURFACE", vec![]),
            ]),
        }
    }
}

#[cfg(test)]
mod tests;
