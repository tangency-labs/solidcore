//! The units of a representation context (ISO 10303-41 measures and units, ISO 10303-43
//! representation contexts): `GLOBAL_UNIT_ASSIGNED_CONTEXT`'s length and plane angle units and
//! `GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT`'s distance uncertainty, read from a parsed
//! [`Exchange`]. A unit that cannot be read is reported as [`Unreadable`] (the entity where
//! reading stopped, and [`Why`]), never read as millimetres or radians: whether to refuse the
//! file or ask for its unit is the importer's decision in `solidcore-ops` (docs/layers.md §1).
//!
//! A unit is a factor and a name: millimetres per file length unit, radians per file plane angle
//! unit. An SI unit's factor is a literal of a table by prefix (a `METRE` is 1000 mm, a `MILLI`
//! `METRE` 1.0, a `NANO` `METRE` 1e-6; a `RADIAN` is 1 rad, a `MILLI` `RADIAN` 1e-3), never a
//! product; a `CONVERSION_BASED_UNIT` is the value of its measure times the factor of the
//! measure's own unit (an inch is `LENGTH_MEASURE(25.4)` of millimetres, a degree
//! `PLANE_ANGLE_MEASURE(0.01745…)` of radians). Only multiplications: the factor of a degree is
//! the number the file writes, not a π computed here (the record computes no geometric quantity,
//! docs/layers.md §4). The name (`'DEGREE'`) is kept next to it, so a caller can tell a degree
//! written as a shortened π/180 from another unit.
//!
//! Origin: SolidCraft's `crates/step/src/units.rs` and the `Units` type of
//! `crates/step/src/read.rs` at `a5b71ca2` (MIT OR Apache-2.0).
//!
//! Deviations from SolidCraft: [`Units`] lives here (it was in `read.rs`, which is not carried
//! over yet). A unit that cannot be read is [`Unreadable`] (SolidCraft reads it as mm or rad, and
//! an SI unit with an unknown base or prefix as factor 1); [`Units`] keeps each unit's name. A
//! `MEASURE_WITH_UNIT` is read by its record's name (SolidCraft takes the first record, the empty
//! one in the alphabetical complex form); a conversion's measure must be of the unit's own kind.
//! SI prefixes are literals (SolidCraft multiplies, so a `NANO` `METRE` is
//! 1.0000000000000002e-6 mm). Only a length measure is the distance uncertainty (SolidCraft takes
//! the last measure listed, an angle accuracy too). The nesting bound of a
//! `CONVERSION_BASED_UNIT` chain is the named `MAX_CONVERSION_DEPTH` (was the literal 4).

use super::p21::{Entity, Exchange, Param};

/// Most `CONVERSION_BASED_UNIT` links followed from one unit (a unit defined from a unit defined
/// from …); a deeper chain, or a cycle, is [`Why::TooDeep`].
const MAX_CONVERSION_DEPTH: usize = 4;

/// A unit read from a file.
#[derive(Clone, Debug, PartialEq)]
pub struct Unit {
    /// Millimetres per unit for a length unit, radians per unit for a plane angle unit: finite
    /// and positive.
    pub factor: f64,
    /// The unit's name: a `CONVERSION_BASED_UNIT`'s as the file writes it (`INCH`, `DEGREE`,
    /// `degree`; empty when it writes none), an SI unit's prefix and base joined (`MILLIMETRE`,
    /// `METRE`, `RADIAN`).
    pub name: String,
}

/// Why a unit cannot be read.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Why {
    /// The reference names no entity in the file.
    Missing,
    /// The context has no `GLOBAL_UNIT_ASSIGNED_CONTEXT` record with a list of units (a context
    /// written with Part 21 short names has `GUAC` instead).
    NoUnits,
    /// The context lists no unit of this kind.
    NotListed,
    /// The entity is not a `NAMED_UNIT` (a unit written with Part 21 short names is an `NMDUNT`).
    NotAUnit,
    /// A length or plane angle unit that is neither an `SI_UNIT` nor a `CONVERSION_BASED_UNIT`.
    NoDefinition,
    /// An SI unit whose name is not the base of its kind (`METRE` for a length, `RADIAN` for a
    /// plane angle).
    SiBase,
    /// An SI unit whose prefix is not one of ISO 10303-41 `si_prefix`.
    SiPrefix,
    /// A measure with unit that cannot be read: no `MEASURE_WITH_UNIT`, a value that is not a
    /// finite number, or a unit that is not a reference.
    Measure,
    /// A measure whose unit is not of the kind it defines (an inch defined in radians).
    KindMismatch,
    /// A factor or uncertainty that is zero, negative or not finite.
    NotPositive,
    /// A `CONVERSION_BASED_UNIT` chain longer than `MAX_CONVERSION_DEPTH` links, or a cycle.
    TooDeep,
}

/// A unit that cannot be read: the entity where reading stopped (the context, a unit or a
/// measure) and why. Reported, not replaced by a default.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Unreadable {
    /// The entity id.
    pub entity: u64,
    /// Why.
    pub why: Why,
}

/// The units of a representation context.
#[derive(Clone, Debug, PartialEq)]
pub struct Units {
    /// The length unit.
    pub len: Result<Unit, Unreadable>,
    /// The plane angle unit.
    pub ang: Result<Unit, Unreadable>,
}

/// What a representation context says about its numbers.
#[derive(Clone, Debug, PartialEq)]
pub struct Context {
    /// The length and plane angle units.
    pub units: Units,
    /// The distance uncertainty in millimetres of the last length `UNCERTAINTY_MEASURE_WITH_UNIT`
    /// the context lists, or why it cannot be read; `None` when it lists none.
    pub uncertainty: Option<Result<f64, Unreadable>>,
}

/// The SI prefixes of ISO 10303-41 `si_prefix`: each prefix's name, the millimetres in its metre
/// and its factor (the radians in its radian), as literals. A literal is the double nearest the
/// exact value; a product such as `1e-9 * 1000.0` is not (it is 1.0000000000000002e-6).
const SI_PREFIXES: [(&str, f64, f64); 16] = [
    ("EXA", 1e21, 1e18),
    ("PETA", 1e18, 1e15),
    ("TERA", 1e15, 1e12),
    ("GIGA", 1e12, 1e9),
    ("MEGA", 1e9, 1e6),
    ("KILO", 1e6, 1e3),
    ("HECTO", 1e5, 1e2),
    ("DECA", 1e4, 1e1),
    ("DECI", 100.0, 1e-1),
    ("CENTI", 10.0, 1e-2),
    ("MILLI", 1.0, 1e-3),
    ("MICRO", 1e-3, 1e-6),
    ("NANO", 1e-6, 1e-9),
    ("PICO", 1e-9, 1e-12),
    ("FEMTO", 1e-12, 1e-15),
    ("ATTO", 1e-15, 1e-18),
];

/// What a unit measures.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum UnitKind {
    /// A `LENGTH_UNIT`.
    Length,
    /// A `PLANE_ANGLE_UNIT`.
    Angle,
}

/// The value and unit of a measure with unit, by name: its `MEASURE_WITH_UNIT` record (in the
/// alphabetical complex form `(LENGTH_MEASURE_WITH_UNIT()MEASURE_WITH_UNIT(…)…)` the first record
/// is the empty subtype one), or the parameters of a simple instance of `MEASURE_WITH_UNIT` or a
/// subtype (`LENGTH_MEASURE_WITH_UNIT`, `UNCERTAINTY_MEASURE_WITH_UNIT`, …), whose first two are
/// the inherited value and unit.
fn measure(e: &Entity) -> Option<&[Param]> {
    e.record("MEASURE_WITH_UNIT").or_else(|| {
        (!e.is_complex() && e.name().ends_with("MEASURE_WITH_UNIT")).then(|| e.params())
    })
}

/// Whether a measure value is typed `LENGTH_MEASURE(…)`.
fn is_length_measure(v: Option<&Param>) -> bool {
    matches!(v, Some(Param::Typed(name, _)) if name == "LENGTH_MEASURE")
}

/// The kind of a named unit: `None` for a unit of another kind (solid angle, mass, …).
fn kind_of(e: &Entity) -> Result<Option<UnitKind>, Why> {
    if e.has("LENGTH_UNIT") {
        Ok(Some(UnitKind::Length))
    } else if e.has("PLANE_ANGLE_UNIT") {
        Ok(Some(UnitKind::Angle))
    } else if e.has("NAMED_UNIT") {
        Ok(None)
    } else {
        Err(Why::NotAUnit)
    }
}

/// The unit `#id`: its kind, factor and name, `None` for a unit that is neither a length nor a
/// plane angle unit, or why it cannot be read.
///
/// # Errors
///
/// [`Unreadable`] when `#id` is no entity or not a named unit, or its length or plane angle
/// unit cannot be read (see [`Why`]).
pub fn unit(ex: &Exchange, id: u64) -> Result<Option<(UnitKind, Unit)>, Unreadable> {
    read(ex, id, 0)
}

fn read(ex: &Exchange, id: u64, depth: usize) -> Result<Option<(UnitKind, Unit)>, Unreadable> {
    let e = ex.get(id).ok_or(Unreadable {
        entity: id,
        why: Why::Missing,
    })?;
    match kind_of(e) {
        Ok(Some(kind)) => read_kind(ex, e, id, kind, depth).map(|u| Some((kind, u))),
        Ok(None) => Ok(None),
        Err(why) => Err(Unreadable { entity: id, why }),
    }
}

/// The length or plane angle unit `e` (`#id`), `depth` conversion links from the unit asked for.
fn read_kind(
    ex: &Exchange,
    e: &Entity,
    id: u64,
    kind: UnitKind,
    depth: usize,
) -> Result<Unit, Unreadable> {
    let at = |entity| move |why| Unreadable { entity, why };
    let unit = if let Some(si) = e.record("SI_UNIT") {
        si_unit(kind, si).map_err(at(id))?
    } else if let Some(cb) = e.record("CONVERSION_BASED_UNIT") {
        if depth >= MAX_CONVERSION_DEPTH {
            return Err(at(id)(Why::TooDeep));
        }
        let mid = cb
            .get(1)
            .and_then(Param::as_ref_id)
            .ok_or_else(|| at(id)(Why::Measure))?;
        let m = ex.get(mid).ok_or_else(|| at(mid)(Why::Missing))?;
        let mp = measure(m).ok_or_else(|| at(mid)(Why::Measure))?;
        let (Some(v), Some(uid)) = (
            mp.first().and_then(Param::as_f64),
            mp.get(1).and_then(Param::as_ref_id),
        ) else {
            return Err(at(mid)(Why::Measure));
        };
        let of = match read(ex, uid, depth + 1)? {
            Some((k, u)) if k == kind => u,
            _ => return Err(at(mid)(Why::KindMismatch)),
        };
        Unit {
            factor: v * of.factor,
            name: cb
                .first()
                .and_then(Param::as_str)
                .unwrap_or_default()
                .to_owned(),
        }
    } else {
        return Err(at(id)(Why::NoDefinition));
    };
    if unit.factor.is_finite() && unit.factor > 0.0 {
        Ok(unit)
    } else {
        Err(at(id)(Why::NotPositive))
    }
}

/// An SI length or plane angle unit from its `SI_UNIT(prefix, name)` record.
fn si_unit(kind: UnitKind, si: &[Param]) -> Result<Unit, Why> {
    let ((UnitKind::Length, Some(base @ "METRE")) | (UnitKind::Angle, Some(base @ "RADIAN"))) =
        (kind, si.get(1).and_then(Param::as_enum))
    else {
        return Err(Why::SiBase);
    };
    let (prefix, mm, rad) = match si.first() {
        Some(Param::Null) => ("", 1000.0, 1.0),
        p => p
            .and_then(Param::as_enum)
            .and_then(|name| SI_PREFIXES.iter().find(|x| x.0 == name))
            .copied()
            .ok_or(Why::SiPrefix)?,
    };
    let factor = match kind {
        UnitKind::Length => mm,
        UnitKind::Angle => rad,
    };
    Ok(Unit {
        factor,
        name: format!("{prefix}{base}"),
    })
}

/// The units and uncertainty of the representation context `#id`.
///
/// Each length or plane angle unit `GLOBAL_UNIT_ASSIGNED_CONTEXT` lists sets that unit, read or
/// [`Unreadable`]; a later one wins. A kind the context lists no unit of is
/// [`Why::NotListed`], or, when it lists an entry whose kind cannot be told (no entity, or not a
/// named unit), that entry's [`Unreadable`]. Other units (solid angle, …) are skipped. A
/// context without `GLOBAL_UNIT_ASSIGNED_CONTEXT` (a `PARAMETRIC_REPRESENTATION_CONTEXT`, whose
/// numbers are parameters, or one written in Part 21 short names) is [`Why::NoUnits`].
///
/// Each length `UNCERTAINTY_MEASURE_WITH_UNIT` that `GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT` lists
/// (its unit a length unit, or its unit unreadable and its value a `LENGTH_MEASURE`) sets the
/// uncertainty to its value times its unit's factor, or to why that cannot be read; the last one
/// wins. Angle and solid angle accuracies, and entries that are not measures, are skipped.
pub fn context(ex: &Exchange, id: u64) -> Context {
    let at = |why| Unreadable { entity: id, why };
    let Some(e) = ex.get(id) else {
        return Context {
            units: Units {
                len: Err(at(Why::Missing)),
                ang: Err(at(Why::Missing)),
            },
            uncertainty: None,
        };
    };
    let units = match e
        .record("GLOBAL_UNIT_ASSIGNED_CONTEXT")
        .and_then(<[Param]>::first)
        .and_then(Param::as_list)
    {
        Some(list) => assigned(ex, list, at(Why::NotListed)),
        None => Units {
            len: Err(at(Why::NoUnits)),
            ang: Err(at(Why::NoUnits)),
        },
    };
    let uncertainty = e
        .record("GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT")
        .and_then(<[Param]>::first)
        .and_then(Param::as_list)
        .and_then(|list| {
            list.iter()
                .rev()
                .filter_map(Param::as_ref_id)
                .find_map(|uid| distance_uncertainty(ex, uid))
        });
    Context { units, uncertainty }
}

/// The units of a `GLOBAL_UNIT_ASSIGNED_CONTEXT` list (`none`: a kind not listed).
fn assigned(ex: &Exchange, list: &[Param], none: Unreadable) -> Units {
    let mut units = Units {
        len: Err(none),
        ang: Err(none),
    };
    let mut unknown = None;
    for uid in list.iter().filter_map(Param::as_ref_id) {
        let Some(e) = ex.get(uid) else {
            unknown.get_or_insert(Unreadable {
                entity: uid,
                why: Why::Missing,
            });
            continue;
        };
        match kind_of(e) {
            Ok(Some(kind)) => {
                let u = read_kind(ex, e, uid, kind, 0);
                match kind {
                    UnitKind::Length => units.len = u,
                    UnitKind::Angle => units.ang = u,
                }
            }
            Ok(None) => {}
            Err(why) => {
                unknown.get_or_insert(Unreadable { entity: uid, why });
            }
        }
    }
    if let Some(u) = unknown {
        for slot in [&mut units.len, &mut units.ang] {
            if *slot == Err(none) {
                *slot = Err(u);
            }
        }
    }
    units
}

/// The distance uncertainty `#id` lists, in millimetres, or why it cannot be read; `None` when
/// `#id` is not a length measure with unit.
fn distance_uncertainty(ex: &Exchange, id: u64) -> Option<Result<f64, Unreadable>> {
    let m = measure(ex.get(id)?)?;
    let uid = m.get(1).and_then(Param::as_ref_id);
    match uid.and_then(|u| ex.get(u)).map(kind_of) {
        Some(Ok(Some(UnitKind::Length))) => {}
        Some(Ok(_)) => return None,
        _ if is_length_measure(m.first()) => {}
        _ => return None,
    }
    let at = |why| Unreadable { entity: id, why };
    let value = || {
        let v = m
            .first()
            .and_then(Param::as_f64)
            .ok_or_else(|| at(Why::Measure))?;
        let uid = uid.ok_or_else(|| at(Why::Measure))?;
        let Some((UnitKind::Length, u)) = read(ex, uid, 0)? else {
            return Err(at(Why::KindMismatch));
        };
        let x = v * u.factor;
        if x.is_finite() && x > 0.0 {
            Ok(x)
        } else {
            Err(at(Why::NotPositive))
        }
    };
    Some(value())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::step::p21;
    use std::f64::consts::PI;
    use std::fmt::Write as _;

    /// Millimetres, radians and steradians as SolidCraft's writer and its topology tests
    /// (`crates/step/src/topology_tests.rs`, `file`) write them.
    const MM: &str = "#1=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));\
         #2=(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.));\
         #3=(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT());";

    /// The units of a CATIA V4 file (`CATIA SOLUTIONS V4 RELEASE 2.1`, the STEP Tools AP203
    /// corpus's `clip.stp` and `socks.stp`), written with Part 21 short names: an inch, a
    /// degree and a steradian in the context `#48`.
    const SHORT_NAMES: &str = "#26=DMNEXP(0.,0.,0.,0.,0.,0.,0.);\
         #27=DMNEXP(1.,0.,0.,0.,0.,0.,0.);\
         #28=(LNGUNT()NMDUNT(*)SUNT(.MILLI.,.METRE.));\
         #29=LMWU(LENGTH_MEASURE(25.4000000000),#28);\
         #30=(CNBSUN('INCH',#29)LNGUNT()NMDUNT(#27));\
         #31=(NMDUNT(*)PLANUN()SUNT($,.RADIAN.));\
         #32=PAMWU(PLANE_ANGLE_MEASURE(0.0174532925199),#31);\
         #33=(CNBSUN('DEGREE',#32)NMDUNT(#26)PLANUN());\
         #34=(NMDUNT(*)SUNT($,.STERADIAN.)SLANUN());\
         #48=(GMRPCN(3)GC((#49))GUAC((#30,#33,#34))RPRCNT('*MASTER','WORKSPACE'));";

    /// The mm per unit of `SI_UNIT(prefix, .METRE.)` and the rad per unit of
    /// `SI_UNIT(prefix, .RADIAN.)` for every prefix of ISO 10303-41 `si_prefix`, as literals,
    /// with the prefix's name.
    const SI_PREFIX_LITERALS: [(&str, &str, f64, f64); 17] = [
        ("$", "", 1000.0, 1.0),
        (".EXA.", "EXA", 1e21, 1e18),
        (".PETA.", "PETA", 1e18, 1e15),
        (".TERA.", "TERA", 1e15, 1e12),
        (".GIGA.", "GIGA", 1e12, 1e9),
        (".MEGA.", "MEGA", 1e9, 1e6),
        (".KILO.", "KILO", 1e6, 1e3),
        (".HECTO.", "HECTO", 1e5, 1e2),
        (".DECA.", "DECA", 1e4, 1e1),
        (".DECI.", "DECI", 100.0, 1e-1),
        (".CENTI.", "CENTI", 10.0, 1e-2),
        (".MILLI.", "MILLI", 1.0, 1e-3),
        (".MICRO.", "MICRO", 1e-3, 1e-6),
        (".NANO.", "NANO", 1e-6, 1e-9),
        (".PICO.", "PICO", 1e-9, 1e-12),
        (".FEMTO.", "FEMTO", 1e-12, 1e-15),
        (".ATTO.", "ATTO", 1e-15, 1e-18),
    ];

    fn parse(data: &str) -> Exchange {
        let t = format!(
            "ISO-10303-21;\nHEADER;FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));ENDSEC;\nDATA;\n{data}\nENDSEC;\nEND-ISO-10303-21;\n"
        );
        match p21::parse(&t) {
            Ok(ex) => ex,
            Err(e) => panic!("{e}\n{t}"),
        }
    }

    /// A context `#100` with the units `units` and the uncertainties `unc` (entity ids).
    fn ctx(data: &str, units: &str, unc: &str) -> Context {
        let c = format!(
            "#100=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT(({unc}))\
             GLOBAL_UNIT_ASSIGNED_CONTEXT(({units}))REPRESENTATION_CONTEXT('',''));"
        );
        context(&parse(&format!("{data}{c}")), 100)
    }

    fn unit_named(factor: f64, name: &str) -> Unit {
        Unit {
            factor,
            name: name.into(),
        }
    }

    fn unreadable<T>(entity: u64, why: Why) -> Result<T, Unreadable> {
        Err(Unreadable { entity, why })
    }

    #[test]
    fn solidcraft_file_context() {
        // SolidCraft's topology test file: mm, rad, sr and an uncertainty of 1e-7 mm.
        let c = ctx(
            &format!(
                "{MM}#4=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07),#1,'distance_accuracy_value','');"
            ),
            "#1,#2,#3",
            "#4",
        );
        assert_eq!(c.units.len, Ok(unit_named(1.0, "MILLIMETRE")));
        assert_eq!(c.units.ang, Ok(unit_named(1.0, "RADIAN")));
        assert_eq!(c.uncertainty, Some(Ok(1e-7)));
        // The steradian is a unit of neither kind.
        assert_eq!(unit(&parse(MM), 3), Ok(None));
    }

    #[test]
    fn si_prefixes_are_literals() {
        // The product the table used to compute is not the literal for a nano metre.
        assert_ne!(1e-9 * 1000.0, 1e-6);
        for (prefix, name, mm, rad) in SI_PREFIX_LITERALS {
            let ex = parse(&format!(
                "#1=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT({prefix},.METRE.));\
                 #2=(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT({prefix},.RADIAN.));"
            ));
            let metre =
                Ok(unit_named(mm, &format!("{name}METRE"))).map(|u| Some((UnitKind::Length, u)));
            let radian =
                Ok(unit_named(rad, &format!("{name}RADIAN"))).map(|u| Some((UnitKind::Angle, u)));
            assert_eq!(unit(&ex, 1), metre, "{prefix}");
            assert_eq!(unit(&ex, 2), radian, "{prefix}");
        }
    }

    #[test]
    fn inches_and_degrees() {
        // An inch as 25.4 mm, a degree as π/180 rad (written as the shortest decimal of
        // PI / 180, so it reads as exactly that double).
        let data = format!(
            "{MM}#10=LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4),#1);\
             #11=(CONVERSION_BASED_UNIT('INCH',#10)LENGTH_UNIT()NAMED_UNIT(#12));\
             #12=DIMENSIONAL_EXPONENTS(1.,0.,0.,0.,0.,0.,0.);\
             #13=PLANE_ANGLE_MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(0.017453292519943295),#2);\
             #14=(CONVERSION_BASED_UNIT('DEGREE',#13)NAMED_UNIT(#15)PLANE_ANGLE_UNIT());\
             #15=DIMENSIONAL_EXPONENTS(0.,0.,0.,0.,0.,0.,0.);"
        );
        let c = ctx(&data, "#11,#14,#3", "");
        assert_eq!(c.units.len, Ok(unit_named(25.4, "INCH")));
        assert_eq!(c.units.ang, Ok(unit_named(PI / 180.0, "DEGREE")));
        // A degree as most of the corpus writes it, a shortened π/180: the file's factor and the
        // name, so a caller can tell.
        let c = ctx(
            "#2=(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.));\
             #13=PLANE_ANGLE_MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(0.0174532925),#2);\
             #14=(CONVERSION_BASED_UNIT('degree',#13)NAMED_UNIT(*)PLANE_ANGLE_UNIT());",
            "#14",
            "",
        );
        assert_eq!(c.units.ang, Ok(unit_named(0.0174532925, "degree")));
        // An inch defined from metres: 0.0254 times 1000 mm.
        let c = ctx(
            "#1=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT($,.METRE.));\
             #10=LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(0.0254),#1);\
             #11=(CONVERSION_BASED_UNIT('INCH',#10)LENGTH_UNIT()NAMED_UNIT(*));",
            "#11",
            "",
        );
        assert_eq!(c.units.len, Ok(unit_named(0.0254 * 1000.0, "INCH")));
        // A foot defined from the inch: a chain of two conversions, 12 times 25.4 mm.
        let c = ctx(
            &format!(
                "{data}#20=LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(12.),#11);\
                 #21=(CONVERSION_BASED_UNIT('FOOT',#20)LENGTH_UNIT()NAMED_UNIT(*));"
            ),
            "#21",
            "",
        );
        assert_eq!(c.units.len, Ok(unit_named(12.0 * 25.4, "FOOT")));
    }

    #[test]
    fn measure_with_unit_in_every_form() {
        // An inch's measure as a simple LENGTH_MEASURE_WITH_UNIT, a simple MEASURE_WITH_UNIT, the
        // alphabetical complex form whose first record is the empty LENGTH_MEASURE_WITH_UNIT(),
        // and a complex measure representation item.
        for measure in [
            "LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4),#1)",
            "MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4),#1)",
            "(LENGTH_MEASURE_WITH_UNIT()MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4),#1))",
            "(LENGTH_MEASURE_WITH_UNIT()MEASURE_REPRESENTATION_ITEM()\
             MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4),#1)REPRESENTATION_ITEM('inch'))",
        ] {
            let c = ctx(
                &format!(
                    "{MM}#10={measure};\
                     #11=(CONVERSION_BASED_UNIT('INCH',#10)LENGTH_UNIT()NAMED_UNIT(*));"
                ),
                "#11,#2",
                "",
            );
            assert_eq!(c.units.len, Ok(unit_named(25.4, "INCH")), "{measure}");
        }
        // A degree in the complex form.
        let c = ctx(
            &format!(
                "{MM}#13=(MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(0.017453292519943295),#2)\
                 PLANE_ANGLE_MEASURE_WITH_UNIT());\
                 #14=(CONVERSION_BASED_UNIT('DEGREE',#13)NAMED_UNIT(*)PLANE_ANGLE_UNIT());"
            ),
            "#1,#14",
            "",
        );
        assert_eq!(c.units.ang, Ok(unit_named(PI / 180.0, "DEGREE")));
        // An uncertainty in the simple and in the complex form.
        for measure in [
            "UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07),#1,'distance_accuracy_value','')",
            "(LENGTH_MEASURE_WITH_UNIT()MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07),#1)\
             UNCERTAINTY_MEASURE_WITH_UNIT('distance_accuracy_value',''))",
        ] {
            let c = ctx(&format!("{MM}#20={measure};"), "#1,#2", "#20");
            assert_eq!(c.uncertainty, Some(Ok(1e-7)), "{measure}");
        }
        // A measure with neither record, or in short names, cannot be read.
        for measure in [
            "(LENGTH_MEASURE_WITH_UNIT()MEASURE_REPRESENTATION_ITEM())",
            "LMWU(LENGTH_MEASURE(25.4),#1)",
        ] {
            let c = ctx(
                &format!(
                    "{MM}#10={measure};\
                     #11=(CONVERSION_BASED_UNIT('INCH',#10)LENGTH_UNIT()NAMED_UNIT(*));"
                ),
                "#11",
                "",
            );
            assert_eq!(c.units.len, unreadable(10, Why::Measure), "{measure}");
        }
    }

    #[test]
    fn distance_uncertainty_is_a_length() {
        // A distance accuracy, then a solid angle and an angle accuracy listed last: the
        // distance uncertainty is the length one.
        let c = ctx(
            &format!(
                "{MM}#20=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07),#1,'distance_accuracy_value','');\
                 #21=UNCERTAINTY_MEASURE_WITH_UNIT(SOLID_ANGLE_MEASURE(1.E-05),#3,'solid_angle_accuracy_value','');\
                 #22=UNCERTAINTY_MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(1.E-03),#2,'angle_accuracy_value','');"
            ),
            "#1,#2,#3",
            "#20,#21,#22",
        );
        assert_eq!(c.uncertainty, Some(Ok(1e-7)));
        // An angle accuracy alone is no distance uncertainty.
        let c = ctx(
            &format!(
                "{MM}#22=UNCERTAINTY_MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(1.E-03),#2,'angle_accuracy_value','');"
            ),
            "#1,#2,#3",
            "#22",
        );
        assert_eq!(c.uncertainty, None);
    }

    #[test]
    fn short_name_units_are_unreadable() {
        // The CATIA V4 context: no GLOBAL_UNIT_ASSIGNED_CONTEXT record (it is GUAC), so its inch
        // is reported, not read as a millimetre.
        let ex = parse(SHORT_NAMES);
        let c = context(&ex, 48);
        assert_eq!(c.units.len, unreadable(48, Why::NoUnits));
        assert_eq!(c.units.ang, unreadable(48, Why::NoUnits));
        assert_eq!(c.uncertainty, None);
        // Its units are not named units (NMDUNT).
        assert_eq!(unit(&ex, 30), unreadable(30, Why::NotAUnit));
        // A long-name context listing the short-name units: the first one whose kind cannot be
        // told stands for both kinds.
        let c = ctx(SHORT_NAMES, "#30,#33,#34", "");
        assert_eq!(c.units.len, unreadable(30, Why::NotAUnit));
        assert_eq!(c.units.ang, unreadable(30, Why::NotAUnit));
    }

    #[test]
    fn broken_units_are_unreadable() {
        // No such context, or a context without units.
        let c = context(&parse(MM), 100);
        assert_eq!(c.units.len, unreadable(100, Why::Missing));
        assert_eq!(c.units.ang, unreadable(100, Why::Missing));
        assert_eq!(c.uncertainty, None);
        assert_eq!(
            context(&parse(MM), 1).units.len,
            unreadable(1, Why::NoUnits)
        );
        // A unit list that is not a list.
        let ex = parse(&format!(
            "{MM}#100=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNIT_ASSIGNED_CONTEXT(#1)REPRESENTATION_CONTEXT('',''));"
        ));
        assert_eq!(context(&ex, 100).units.len, unreadable(100, Why::NoUnits));
        // Each broken length unit `#1` and why.
        for (data, at, why) in [
            // An SI unit with an unknown base, the base of the other kind, or an unknown prefix.
            (
                "#1=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.FURLONG.));",
                1,
                Why::SiBase,
            ),
            (
                "#1=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT($,.RADIAN.));",
                1,
                Why::SiBase,
            ),
            (
                "#1=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.KIBI.,.METRE.));",
                1,
                Why::SiPrefix,
            ),
            // A length unit that is neither an SI unit nor a conversion.
            (
                "#1=(CONTEXT_DEPENDENT_UNIT('x')LENGTH_UNIT()NAMED_UNIT(*));",
                1,
                Why::NoDefinition,
            ),
            // A conversion whose measure is missing, unreadable, zero or infinite.
            (
                "#1=(CONVERSION_BASED_UNIT('INCH',#9)LENGTH_UNIT()NAMED_UNIT(*));",
                9,
                Why::Missing,
            ),
            (
                "#2=(NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.)LENGTH_UNIT());\
                 #3=LENGTH_MEASURE_WITH_UNIT($,#2);\
                 #1=(CONVERSION_BASED_UNIT('INCH',#3)LENGTH_UNIT()NAMED_UNIT(*));",
                3,
                Why::Measure,
            ),
            (
                "#2=(NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.)LENGTH_UNIT());\
                 #3=LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(0.),#2);\
                 #1=(CONVERSION_BASED_UNIT('INCH',#3)LENGTH_UNIT()NAMED_UNIT(*));",
                1,
                Why::NotPositive,
            ),
            (
                "#2=(NAMED_UNIT(*)SI_UNIT(.EXA.,.METRE.)LENGTH_UNIT());\
                 #3=LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E300),#2);\
                 #1=(CONVERSION_BASED_UNIT('BIG',#3)LENGTH_UNIT()NAMED_UNIT(*));",
                1,
                Why::NotPositive,
            ),
            // An inch defined in radians.
            (
                "#2=(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.));\
                 #3=LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4),#2);\
                 #1=(CONVERSION_BASED_UNIT('INCH',#3)LENGTH_UNIT()NAMED_UNIT(*));",
                3,
                Why::KindMismatch,
            ),
            // A conversion that refers to itself: the depth bound ends it.
            (
                "#3=LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(2.),#1);\
                 #1=(CONVERSION_BASED_UNIT('LOOP',#3)LENGTH_UNIT()NAMED_UNIT(*));",
                1,
                Why::TooDeep,
            ),
        ] {
            let c = ctx(data, "#1", "");
            assert_eq!(c.units.len, unreadable(at, why), "{data}");
            assert_eq!(c.units.ang, unreadable(100, Why::NotListed), "{data}");
        }
        // A negative degree.
        let c = ctx(
            "#2=(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.));\
             #3=PLANE_ANGLE_MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(-1.),#2);\
             #1=(CONVERSION_BASED_UNIT('DEGREE',#3)PLANE_ANGLE_UNIT()NAMED_UNIT(*));",
            "#1",
            "",
        );
        assert_eq!(c.units.ang, unreadable(1, Why::NotPositive));
        // Units that are neither lengths nor plane angles, or none: not listed.
        for units in ["#3", ""] {
            let c = ctx(MM, units, "");
            assert_eq!(c.units.len, unreadable(100, Why::NotListed), "{units}");
            assert_eq!(c.units.ang, unreadable(100, Why::NotListed), "{units}");
        }
        // Dangling references: the first stands for both kinds.
        let c = ctx("", "#7,#8", "");
        assert_eq!(c.units.len, unreadable(7, Why::Missing));
        assert_eq!(c.units.ang, unreadable(7, Why::Missing));
        // A broken length unit leaves the readable angle unit in place.
        let c = ctx(
            "#2=(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT(.MILLI.,.RADIAN.));\
             #1=(CONVERSION_BASED_UNIT('INCH',#9)LENGTH_UNIT()NAMED_UNIT(*));",
            "#1,#2",
            "",
        );
        assert_eq!(c.units.len, unreadable(9, Why::Missing));
        assert_eq!(c.units.ang, Ok(unit_named(1e-3, "MILLIRADIAN")));
        // The last length unit listed wins, read or not.
        let broken = "#4=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.FURLONG.));";
        let c = ctx(&format!("{MM}{broken}"), "#1,#4", "");
        assert_eq!(c.units.len, unreadable(4, Why::SiBase));
        let c = ctx(&format!("{MM}{broken}"), "#4,#1", "");
        assert_eq!(c.units.len, Ok(unit_named(1.0, "MILLIMETRE")));
    }

    #[test]
    fn conversion_depth_is_bounded() {
        // #1 is millimetres; #2..=#7 each 2 of the one before, so #(1+k) is 2^k mm and
        // MAX_CONVERSION_DEPTH links are followed from the unit asked for.
        let mut data = String::from("#1=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));");
        for k in 2..=7 {
            write!(
                data,
                "#{m}=LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(2.),#{p});\
                 #{k}=(CONVERSION_BASED_UNIT('U{k}',#{m})LENGTH_UNIT()NAMED_UNIT(*));",
                m = 100 + k,
                p = k - 1
            )
            .unwrap();
        }
        let ex = parse(&data);
        assert_eq!(
            unit(&ex, 5),
            Ok(unit_named(16.0, "U5")).map(|u| Some((UnitKind::Length, u)))
        );
        // From #6 the fifth link would be #2's.
        assert_eq!(unit(&ex, 6), unreadable(2, Why::TooDeep));
    }

    #[test]
    fn uncertainty_in_file_units() {
        let inch = format!(
            "{MM}#10=LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4),#1);\
             #11=(CONVERSION_BASED_UNIT('INCH',#10)LENGTH_UNIT()NAMED_UNIT(*));"
        );
        // 1e-4 inch is 1e-4 * 25.4 mm.
        let c = ctx(
            &format!("{inch}#20=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-04),#11,'d','');"),
            "#11,#2",
            "#20",
        );
        assert_eq!(c.units.len, Ok(unit_named(25.4, "INCH")));
        assert_eq!(c.uncertainty, Some(Ok(1e-4 * 25.4)));
        // A measure in metres in an inch context: its own unit, 1e-6 m = 1e-6 * 1000 mm.
        let c = ctx(
            &format!(
                "{inch}#12=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT($,.METRE.));\
                 #20=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-06),#12,'d','');"
            ),
            "#11",
            "#20",
        );
        assert_eq!(c.uncertainty, Some(Ok(1e-6 * 1000.0)));
        // A length measure whose unit cannot be read is reported, not read in the context's
        // length unit.
        let c = ctx(
            &format!("{inch}#20=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-04),$,'d','');"),
            "#11",
            "#20",
        );
        assert_eq!(c.uncertainty, Some(unreadable(20, Why::Measure)));
        let c = ctx(
            &format!("{inch}#20=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-04),#99,'d','');"),
            "#11",
            "#20",
        );
        assert_eq!(c.uncertainty, Some(unreadable(99, Why::Missing)));
        // The last one listed wins.
        let c = ctx(
            &format!(
                "{MM}#20=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-04),#1,'d','');\
                 #21=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(2.E-04),#1,'d','');"
            ),
            "#1",
            "#20,#21",
        );
        assert_eq!(c.uncertainty, Some(Ok(2e-4)));
    }

    #[test]
    fn uncertainty_zero_negative_or_infinite_is_unreadable() {
        for v in ["0.", "-1.E-04", "1.E300"] {
            // In exametres (1e21 mm): 0, negative, and 1e300 * 1e21 mm, which is infinite.
            let c = ctx(
                &format!(
                    "{MM}#12=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.EXA.,.METRE.));\
                     #20=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE({v}),#12,'d','');"
                ),
                "#1",
                "#20",
            );
            assert_eq!(c.uncertainty, Some(unreadable(20, Why::NotPositive)), "{v}");
        }
        // A bad one after a good one is reported (the last one wins).
        let good = "#20=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-04),#1,'d','');";
        let c = ctx(
            &format!("{MM}{good}#21=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(0.),#1,'d','');"),
            "#1",
            "#20,#21",
        );
        assert_eq!(c.uncertainty, Some(unreadable(21, Why::NotPositive)));
        // One without a readable value, in a length unit, too.
        let c = ctx(
            &format!("{MM}{good}#21=UNCERTAINTY_MEASURE_WITH_UNIT($,#1,'d','');"),
            "#1",
            "#20,#21",
        );
        assert_eq!(c.uncertainty, Some(unreadable(21, Why::Measure)));
        // A dangling reference, or an entity that is no measure, is skipped.
        let c = ctx(&format!("{MM}{good}"), "#1", "#20,#99,#2");
        assert_eq!(c.uncertainty, Some(Ok(1e-4)));
    }
}
