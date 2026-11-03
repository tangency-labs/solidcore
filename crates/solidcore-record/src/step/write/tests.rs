//! The writer, one entity kind at a time: record geometry written, parsed again with
//! [`p21::parse`] and read with [`Reader`] in millimetres and radians gives the same record,
//! bit for bit (compared with `==` and by the `{:?}` text, which tells every `f64` apart,
//! `-0.0` from `0.0` included). Plus the formatting of reals and strings, the sharing of
//! identical entities, the ids, and what the writer refuses.
//!
//! Origin: the writer cases of SolidCraft's `crates/step/src/geom_tests.rs` at `a5b71ca2` (a
//! write and re-read of each kind), on the record's types; the formatting, sharing and refusal
//! cases are new.

use std::sync::Arc;

use super::{GeomWriter, R, Writer, real, string};
use crate::step::geom::{
    Associated, Axis1, BSplineCurve, BSplineSurface, Curve, CurveForm, KnotSpec, Logical, P3,
    Pcurve, Placement, PreferredRep, Surface, SurfaceCurveKind, SurfaceForm, Trim, TrimPreference,
    Vector,
};
use crate::step::p21::{self, Exchange, Param};
use crate::step::read::Reader;
use crate::step::units::Units;

const SCHEMA: &str = "AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }";

/// Values whose shortest form is long, tiny, negative zero or has an exponent: every one must
/// come back bit for bit.
const AWKWARD: [f64; 6] = [0.1, 1.0 / 3.0, -0.0, 1e-300, 123456789.123456789, -2.5e-7];

fn parse(w: &Writer) -> Exchange {
    let t = w.file("roundtrip.stp", SCHEMA);
    match p21::parse(&t) {
        Ok(ex) => ex,
        Err(e) => panic!("{e}\n{t}"),
    }
}

/// `c` written, parsed and read back; the same record exactly.
fn curve_back(c: &Curve) -> Arc<Curve> {
    let mut w = Writer::new();
    let id = GeomWriter::new(&mut w).curve(c).unwrap();
    let ex = parse(&w);
    let back = Reader::new(&ex, Units::default()).curve(id).unwrap();
    assert_eq!(*back, *c, "\n{}", w.data());
    assert_eq!(format!("{back:?}"), format!("{c:?}"), "bits");
    back
}

/// `sf` written, parsed and read back; the same record exactly.
fn surface_back(sf: &Surface) -> Arc<Surface> {
    let mut w = Writer::new();
    let id = GeomWriter::new(&mut w).surface(sf).unwrap();
    let ex = parse(&w);
    let back = Reader::new(&ex, Units::default()).surface(id).unwrap();
    assert_eq!(*back, *sf, "\n{}", w.data());
    assert_eq!(format!("{back:?}"), format!("{sf:?}"), "bits");
    back
}

/// `pc` written, parsed and read back; the same record exactly.
fn pcurve_back(pc: &Pcurve) {
    let mut w = Writer::new();
    let id = GeomWriter::new(&mut w).pcurve(pc).unwrap();
    let ex = parse(&w);
    let back = Reader::new(&ex, Units::default()).pcurve(id).unwrap();
    assert_eq!(back, *pc, "\n{}", w.data());
    assert_eq!(format!("{back:?}"), format!("{pc:?}"), "bits");
}

fn curve_err(c: &Curve) -> R<u64> {
    GeomWriter::new(&mut Writer::new()).curve(c)
}

const WORLD: Placement = Placement {
    origin: [0.0, 0.0, 0.0],
    axis: Some([0.0, 0.0, 1.0]),
    ref_direction: Some([1.0, 0.0, 0.0]),
};

/// A placement with unnormalised ratios, an axis not orthogonal to its reference direction,
/// and awkward reals: written as held.
const SKEW: Placement = Placement {
    origin: [0.1, -0.0, 1e-300],
    axis: Some([0.0, 0.0, 2.0]),
    ref_direction: Some([3.0, 0.0, 1.0 / 3.0]),
};

/// Model-space vector in millimetres.
fn vector(orientation: P3, magnitude: f64) -> Vector {
    Vector {
        orientation,
        magnitude,
        scale: [1.0; 3],
    }
}

fn line(origin: P3, orientation: P3, magnitude: f64) -> Curve {
    Curve::Line {
        origin,
        vector: vector(orientation, magnitude),
    }
}

fn bspline_curve(weights: Option<Vec<f64>>, points: Vec<P3>) -> BSplineCurve {
    BSplineCurve {
        degree: 2,
        points,
        weights,
        knots: vec![0.0, 0.1, 1.0 / 3.0],
        multiplicities: vec![3, 1, 3],
        form: CurveForm::CircularArc,
        closed: Logical::False,
        self_intersect: Logical::Unknown,
        knot_spec: KnotSpec::Unspecified,
    }
}

const NET: [[P3; 2]; 2] = [
    [[0.0, 0.0, 0.0], [0.0, 1.0, 0.1]],
    [[1.0, 0.0, -0.0], [1.0, 1.0, 1.0 / 3.0]],
];

fn bspline_surface(weights: Option<Vec<Vec<f64>>>) -> BSplineSurface {
    BSplineSurface {
        u_degree: 1,
        v_degree: 1,
        points: NET.iter().map(|r| r.to_vec()).collect(),
        weights,
        u_knots: vec![0.0, 1.0],
        u_multiplicities: vec![2, 2],
        v_knots: vec![-0.5, 0.25],
        v_multiplicities: vec![2, 2],
        form: SurfaceForm::RuledSurf,
        u_closed: Logical::False,
        v_closed: Logical::True,
        self_intersect: Logical::Unknown,
        knot_spec: KnotSpec::PiecewiseBezierKnots,
    }
}

fn cylinder() -> Arc<Surface> {
    Arc::new(Surface::Cylinder {
        position: WORLD,
        radius: 5.0,
    })
}

// Formatting.

#[test]
fn reals_are_written_in_their_shortest_exact_form() {
    let cases = [
        (0.0, "0.0"),
        (-0.0, "-0.0"),
        (1.0, "1.0"),
        (100.0, "100.0"),
        (-3.0, "-3.0"),
        (0.1, "0.1"),
        (1e-7, "1.E-7"),
        (-2.5e-7, "-2.5E-7"),
        (1e16, "1.E16"),
        (1.5e300, "1.5E300"),
        (5e-324, "5.E-324"),
    ];
    for (x, t) in cases {
        assert_eq!(real(x).as_deref(), Some(t), "{x:?}");
    }
    for x in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(real(x), None, "{x:?}");
    }
    // Every one reads back to the same bits, through the Part 21 parser.
    let values: Vec<f64> = AWKWARD
        .into_iter()
        .chain(cases.map(|c| c.0))
        .chain([
            f64::MAX,
            f64::MIN,
            f64::MIN_POSITIVE,
            f64::EPSILON,
            0.1 + 0.2,
            std::f64::consts::PI,
        ])
        .collect();
    let mut w = Writer::new();
    let id = w
        .add(
            "X",
            &[Param::List(
                values.iter().map(|x| Param::Real(*x)).collect(),
            )],
        )
        .unwrap();
    let ex = parse(&w);
    let back: Vec<u64> = ex.get(id).unwrap().params()[0]
        .as_list()
        .unwrap()
        .iter()
        .map(|p| p.as_f64().unwrap().to_bits())
        .collect();
    assert_eq!(back, values.iter().map(|x| x.to_bits()).collect::<Vec<_>>());
}

#[test]
fn strings_escape_quotes_backslashes_and_non_ascii() {
    let cases = [
        ("", "''"),
        ("plain text 1", "'plain text 1'"),
        ("it's", "'it''s'"),
        ("a\\b", "'a\\\\b'"),
        ("\\X2\\", "'\\\\X2\\\\'"),
        ("é", "'\\X2\\00E9\\X0\\'"),
        ("Größe", "'Gr\\X2\\00F6\\X0\\\\X2\\00DF\\X0\\e'"),
        ("😀", "'\\X2\\D83DDE00\\X0\\'"),
        ("a\nb\t", "'a\\X\\0Ab\\X\\09'"),
    ];
    for (s, t) in cases {
        assert_eq!(string(s), t, "{s:?}");
    }
    // Each reads back as it was.
    let mut w = Writer::new();
    let ids: Vec<u64> = cases
        .iter()
        .map(|(s, _)| w.add_distinct("NAMED", &[Param::Str((*s).into())]).unwrap())
        .collect();
    let ex = parse(&w);
    for ((s, _), id) in cases.iter().zip(ids) {
        assert_eq!(ex.get(id).unwrap().params()[0].as_str(), Some(*s));
    }
    // The file's name and schema are strings too.
    let ex = p21::parse(&Writer::new().file("Größe 'A'.stp", SCHEMA)).unwrap();
    assert_eq!(ex.schemas, vec![SCHEMA.to_string()]);
}

#[test]
fn ids_are_dense_in_insertion_order_and_identical_entities_shared() {
    let mut w = Writer::new();
    let mut g = GeomWriter::new(&mut w);
    let p = [1.0, 2.0, 3.0];
    let q = [4.0, 5.0, 6.0];
    assert_eq!(g.curve(&Curve::Polyline(vec![p, q, p])).unwrap(), 3);
    assert_eq!(g.curve(&line(p, [1.0, 0.0, 0.0], 5.0)).unwrap(), 6);
    // The same line again is the same entity; the same point is the same point.
    assert_eq!(g.curve(&line(p, [1.0, 0.0, 0.0], 5.0)).unwrap(), 6);
    assert_eq!(g.point(q).unwrap(), 2);
    assert_eq!(
        w.data(),
        "#1=CARTESIAN_POINT('',(1.0,2.0,3.0));\n\
         #2=CARTESIAN_POINT('',(4.0,5.0,6.0));\n\
         #3=POLYLINE('',(#1,#2,#1));\n\
         #4=DIRECTION('',(1.0,0.0,0.0));\n\
         #5=VECTOR('',#4,5.0);\n\
         #6=LINE('',#1,#5);\n"
    );
    // A distinct entity is written again, and is not shared afterwards.
    let a = w.add_distinct("VERTEX_POINT", &[Param::Str(String::new()), Param::Ref(1)]);
    let b = w.add_distinct("VERTEX_POINT", &[Param::Str(String::new()), Param::Ref(1)]);
    assert_eq!((a, b), (Ok(7), Ok(8)));
    assert_eq!(w.len(), 8);
    // A complex instance writes its records in alphabetical order.
    let c = w
        .add_complex(&[("B", vec![Param::Int(2)]), ("A", vec![])])
        .unwrap();
    assert!(w.data().ends_with("#9=(A()B(2));\n"));
    assert_eq!(c, 9);
    // The same calls give the same text.
    let again = {
        let mut w2 = Writer::new();
        let mut g = GeomWriter::new(&mut w2);
        g.curve(&Curve::Polyline(vec![p, q, p])).unwrap();
        g.curve(&line(p, [1.0, 0.0, 0.0], 5.0)).unwrap();
        w2.data()
    };
    assert!(w.data().starts_with(&again));
}

// Curves.

#[test]
fn lines_and_conics() {
    curve_back(&line([0.1, -0.0, 1e-300], [3.0, 4.0, 0.0], 1.0 / 3.0));
    curve_back(&Curve::Circle {
        position: SKEW,
        radius: 0.1,
    });
    curve_back(&Curve::Ellipse {
        position: WORLD,
        semi_axis_1: 10.0,
        semi_axis_2: 1.0 / 3.0,
    });
    curve_back(&Curve::Hyperbola {
        position: SKEW,
        semi_axis: 2.0,
        semi_imag_axis: 3.0,
    });
    curve_back(&Curve::Parabola {
        position: WORLD,
        focal_dist: -0.5,
    });
    // Axes the file omits are written as `$` and read back as absent.
    let back = curve_back(&Curve::Circle {
        position: Placement {
            origin: [1.0, 2.0, 3.0],
            axis: None,
            ref_direction: None,
        },
        radius: 2.0,
    });
    assert!(matches!(*back, Curve::Circle { position, .. } if position.axis.is_none()));
}

#[test]
fn bspline_curves() {
    let pts = vec![
        [0.0, 0.0, 0.0],
        [0.1, 1.0 / 3.0, -0.0],
        [2.0, 0.0, 1e-300],
        [3.0, 1.0, 0.0],
    ];
    curve_back(&Curve::BSpline(Box::new(bspline_curve(None, pts.clone()))));
    curve_back(&Curve::BSpline(Box::new(bspline_curve(
        Some(vec![1.0, 0.1, 1.0 / 3.0, 2.0]),
        pts.clone(),
    ))));
    // The rational one is a complex instance, its records in alphabetical order.
    let mut w = Writer::new();
    GeomWriter::new(&mut w)
        .curve(&Curve::BSpline(Box::new(bspline_curve(
            Some(vec![1.0; 4]),
            pts,
        ))))
        .unwrap();
    assert!(
        w.data()
            .contains("=(BOUNDED_CURVE()B_SPLINE_CURVE(2,(#1,#2,#3,#4),.CIRCULAR_ARC.,.F.,.U.)B_SPLINE_CURVE_WITH_KNOTS((3,1,3),(0.0,0.1,0.3333333333333333),.UNSPECIFIED.)CURVE()GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_CURVE((1.0,1.0,1.0,1.0))REPRESENTATION_ITEM(''));"),
        "{}",
        w.data()
    );
    // Every knot type and form is written as its own enumeration value.
    for knot_spec in [
        KnotSpec::UniformKnots,
        KnotSpec::QuasiUniformKnots,
        KnotSpec::PiecewiseBezierKnots,
    ] {
        for form in [
            CurveForm::PolylineForm,
            CurveForm::EllipticArc,
            CurveForm::ParabolicArc,
            CurveForm::HyperbolicArc,
            CurveForm::Unspecified,
        ] {
            let mut b = bspline_curve(
                None,
                vec![[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
            );
            b.knot_spec = knot_spec;
            b.form = form;
            b.closed = Logical::True;
            curve_back(&Curve::BSpline(Box::new(b)));
        }
    }
}

#[test]
fn polylines_and_trimmed_curves() {
    curve_back(&Curve::Polyline(vec![
        [0.0; 3],
        [0.1, 0.2, 0.3],
        [1.0 / 3.0, -0.0, 1.0],
    ]));
    let circle = Arc::new(Curve::Circle {
        position: WORLD,
        radius: 2.0,
    });
    for preference in [
        TrimPreference::Cartesian,
        TrimPreference::Parameter,
        TrimPreference::Unspecified,
    ] {
        curve_back(&Curve::Trimmed {
            basis: circle.clone(),
            trim_1: vec![Trim::Point([2.0, 0.0, 0.0]), Trim::Parameter(0.0)],
            trim_2: vec![Trim::Parameter(1.0 / 3.0)],
            sense: false,
            preference,
        });
    }
    curve_back(&Curve::Trimmed {
        basis: Arc::new(line([0.0; 3], [1.0, 1.0, 0.0], 2.0)),
        trim_1: vec![Trim::Parameter(-0.0)],
        trim_2: vec![Trim::Point([0.1, 0.1, 0.0])],
        sense: true,
        preference: TrimPreference::Parameter,
    });
}

#[test]
fn surface_curves_with_pcurves_and_surfaces() {
    let cyl = cylinder();
    let up = vector([0.0, 1.0, 0.0], 1.0);
    let pc = |u: f64| {
        Associated::Pcurve(Pcurve {
            surface: cyl.clone(),
            curve: Arc::new(Curve::Line {
                origin: [u, 0.0, 0.0],
                vector: up,
            }),
        })
    };
    for (kind, master) in [
        (SurfaceCurveKind::Seam, PreferredRep::PcurveS1),
        (SurfaceCurveKind::Surface, PreferredRep::Curve3d),
        (SurfaceCurveKind::Intersection, PreferredRep::PcurveS2),
        (SurfaceCurveKind::BoundedSurface, PreferredRep::Curve3d),
    ] {
        curve_back(&Curve::SurfaceCurve {
            kind,
            curve_3d: Arc::new(line([5.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0)),
            associated: vec![pc(0.0), pc(std::f64::consts::TAU)],
            master,
        });
    }
    // A surface as associated geometry, and a surface curve without any.
    curve_back(&Curve::SurfaceCurve {
        kind: SurfaceCurveKind::Intersection,
        curve_3d: Arc::new(line([5.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0)),
        associated: vec![
            Associated::Surface(cyl.clone()),
            Associated::Surface(Arc::new(Surface::Plane { position: SKEW })),
        ],
        master: PreferredRep::Curve3d,
    });
    curve_back(&Curve::SurfaceCurve {
        kind: SurfaceCurveKind::Surface,
        curve_3d: Arc::new(line([0.0; 3], [1.0, 0.0, 0.0], 1.0)),
        associated: vec![],
        master: PreferredRep::Curve3d,
    });
    // Both pcurves share one surface entity and one parametric context.
    let mut w = Writer::new();
    GeomWriter::new(&mut w)
        .curve(&Curve::SurfaceCurve {
            kind: SurfaceCurveKind::Seam,
            curve_3d: Arc::new(line([5.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0)),
            associated: vec![pc(0.0), pc(1.0)],
            master: PreferredRep::PcurveS1,
        })
        .unwrap();
    let d = w.data();
    assert_eq!(d.matches("CYLINDRICAL_SURFACE").count(), 1, "{d}");
    assert_eq!(
        d.matches("PARAMETRIC_REPRESENTATION_CONTEXT").count(),
        1,
        "{d}"
    );
    assert_eq!(d.matches("=PCURVE(").count(), 2, "{d}");
}

#[test]
fn offset_curves() {
    for self_intersect in [Logical::True, Logical::False, Logical::Unknown] {
        curve_back(&Curve::OffsetCurve3d {
            basis: Arc::new(Curve::Circle {
                position: WORLD,
                radius: 3.0,
            }),
            distance: -0.1,
            self_intersect,
            ref_direction: [0.0, 0.0, 2.0],
        });
    }
}

// Surfaces.

#[test]
fn analytic_surfaces() {
    surface_back(&Surface::Plane { position: SKEW });
    surface_back(&Surface::Plane {
        position: Placement {
            origin: [1.0, 2.0, 3.0],
            axis: Some([0.0, 1.0, 1.0]),
            ref_direction: None,
        },
    });
    surface_back(&Surface::Cylinder {
        position: WORLD,
        radius: 1.0 / 3.0,
    });
    surface_back(&Surface::Cone {
        position: SKEW,
        radius: 0.0,
        semi_angle: std::f64::consts::FRAC_PI_4,
    });
    surface_back(&Surface::Cone {
        position: WORLD,
        radius: 2.5,
        semi_angle: -0.1,
    });
    surface_back(&Surface::Sphere {
        position: WORLD,
        radius: 0.1,
    });
    surface_back(&Surface::Torus {
        position: SKEW,
        major_radius: -2.0,
        minor_radius: 3.0,
    });
    for select_outer in [true, false] {
        surface_back(&Surface::DegenerateTorus {
            position: WORLD,
            major_radius: 1.0,
            minor_radius: 2.0,
            select_outer,
        });
    }
}

#[test]
fn bspline_surfaces() {
    surface_back(&Surface::BSpline(Box::new(bspline_surface(None))));
    surface_back(&Surface::BSpline(Box::new(bspline_surface(Some(vec![
        vec![1.0, 0.1],
        vec![1.0 / 3.0, 4.0],
    ])))));
    for form in [
        SurfaceForm::PlaneSurf,
        SurfaceForm::CylindricalSurf,
        SurfaceForm::ConicalSurf,
        SurfaceForm::SphericalSurf,
        SurfaceForm::ToroidalSurf,
        SurfaceForm::SurfOfRevolution,
        SurfaceForm::GeneralisedCone,
        SurfaceForm::QuadricSurf,
        SurfaceForm::SurfOfLinearExtrusion,
        SurfaceForm::Unspecified,
    ] {
        let mut b = bspline_surface(None);
        b.form = form;
        surface_back(&Surface::BSpline(Box::new(b)));
    }
}

#[test]
fn swept_offset_and_trimmed_surfaces() {
    let profile = Arc::new(Curve::Circle {
        position: Placement {
            origin: [10.0, 0.0, 0.0],
            axis: Some([0.0, 1.0, 0.0]),
            ref_direction: Some([1.0, 0.0, 0.0]),
        },
        radius: 2.0,
    });
    surface_back(&Surface::Revolution {
        profile: profile.clone(),
        axis: Axis1 {
            origin: [0.0; 3],
            direction: Some([0.0, 0.0, 3.0]),
        },
    });
    surface_back(&Surface::Revolution {
        profile: profile.clone(),
        axis: Axis1 {
            origin: [0.1, 0.0, 0.0],
            direction: None,
        },
    });
    surface_back(&Surface::Extrusion {
        profile: profile.clone(),
        vector: vector([0.0, 2.0, 2.0], 1.0 / 3.0),
    });
    surface_back(&Surface::Offset {
        basis: cylinder(),
        distance: -0.1,
        self_intersect: Logical::False,
    });
    for (u_sense, v_sense) in [(true, false), (false, true)] {
        surface_back(&Surface::RectangularTrimmed {
            basis: cylinder(),
            u1: 0.0,
            u2: 1.0 / 3.0,
            v1: -0.0,
            v2: 10.0,
            u_sense,
            v_sense,
        });
    }
}

// Pcurves.

#[test]
fn pcurves_in_parameter_space() {
    let on = |surface: Arc<Surface>, curve: Curve| Pcurve {
        surface,
        curve: Arc::new(curve),
    };
    let flat = Placement {
        origin: [0.1, 1.0 / 3.0, 0.0],
        axis: None,
        ref_direction: Some([2.0, 1.0, 0.0]),
    };
    pcurve_back(&on(
        cylinder(),
        Curve::Line {
            origin: [1.0 / 3.0, -0.0, 0.0],
            vector: vector([0.0, 1.0, 0.0], 2.0),
        },
    ));
    let plane = Arc::new(Surface::Plane { position: WORLD });
    pcurve_back(&on(
        plane.clone(),
        Curve::Circle {
            position: flat,
            radius: 0.1,
        },
    ));
    pcurve_back(&on(
        plane.clone(),
        Curve::Ellipse {
            position: Placement {
                ref_direction: None,
                ..flat
            },
            semi_axis_1: 2.0,
            semi_axis_2: 1.0,
        },
    ));
    pcurve_back(&on(
        Arc::new(Surface::BSpline(Box::new(bspline_surface(None)))),
        Curve::BSpline(Box::new(bspline_curve(
            Some(vec![1.0, 2.0, 1.0, 0.5]),
            vec![
                [0.0, 0.0, 0.0],
                [0.1, 0.2, 0.0],
                [1.0 / 3.0, 1.0, 0.0],
                [1.0, 1.0, 0.0],
            ],
        ))),
    ));
    pcurve_back(&on(
        plane.clone(),
        Curve::Polyline(vec![[0.0, 0.0, 0.0], [1.0, 0.1, 0.0]]),
    ));
    pcurve_back(&on(
        plane,
        Curve::Trimmed {
            basis: Arc::new(Curve::Line {
                origin: [0.0; 3],
                vector: vector([1.0, 0.0, 0.0], 1.0),
            }),
            trim_1: vec![Trim::Parameter(0.0), Trim::Point([0.0; 3])],
            trim_2: vec![Trim::Parameter(1.0 / 3.0)],
            sense: true,
            preference: TrimPreference::Cartesian,
        },
    ));
    // 2D points have two coordinates.
    let mut w = Writer::new();
    GeomWriter::new(&mut w)
        .pcurve(&on(
            cylinder(),
            Curve::Polyline(vec![[0.5, 2.0, 0.0], [1.0, 2.0, 0.0]]),
        ))
        .unwrap();
    assert!(
        w.data().contains("CARTESIAN_POINT('',(0.5,2.0))"),
        "{}",
        w.data()
    );
}

// Refusals.

#[test]
fn refusals_say_why() {
    let inch = Vector {
        orientation: [1.0, 0.0, 0.0],
        magnitude: 2.0,
        scale: [25.4; 3],
    };
    assert_eq!(
        curve_err(&Curve::Line {
            origin: [0.0; 3],
            vector: inch,
        }),
        Err("a VECTOR with unit factors [25.4, 25.4, 25.4]: converting it to millimetres and radians is the kernel's".into())
    );
    assert_eq!(
        GeomWriter::new(&mut Writer::new()).surface(&Surface::Extrusion {
            profile: Arc::new(line([0.0; 3], [1.0, 0.0, 0.0], 1.0)),
            vector: Vector {
                scale: [1.0, 1.0, 0.0174532925199433],
                ..inch
            },
        }),
        Err("a VECTOR with unit factors [1.0, 1.0, 0.0174532925199433]: converting it to millimetres and radians is the kernel's".into())
    );
    assert_eq!(
        curve_err(&Curve::Circle {
            position: WORLD,
            radius: f64::NAN,
        }),
        Err("NaN is not a real".into())
    );
    assert_eq!(
        curve_err(&line([f64::INFINITY, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0)),
        Err("inf is not a real".into())
    );
    let mut w = Writer::new();
    let mut g = GeomWriter::new(&mut w);
    assert_eq!(
        g.pcurve(&Pcurve {
            surface: cylinder(),
            curve: Arc::new(Curve::Polyline(vec![[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]])),
        }),
        Err("2D coordinates [0.0, 0.0, 1.0] with z other than 0".into())
    );
    assert_eq!(
        g.pcurve(&Pcurve {
            surface: cylinder(),
            curve: Arc::new(Curve::Circle {
                position: WORLD,
                radius: 1.0,
            }),
        }),
        Err("a placement in a parameter space has no axis".into())
    );
    assert_eq!(
        g.pcurve(&Pcurve {
            surface: cylinder(),
            curve: Arc::new(Curve::Trimmed {
                basis: Arc::new(Curve::OffsetCurve3d {
                    basis: Arc::new(line([0.0; 3], [1.0, 0.0, 0.0], 1.0)),
                    distance: 1.0,
                    self_intersect: Logical::False,
                    ref_direction: [0.0, 0.0, 1.0],
                }),
                trim_1: vec![Trim::Parameter(0.0)],
                trim_2: vec![Trim::Parameter(1.0)],
                sense: true,
                preference: TrimPreference::Parameter,
            }),
        }),
        Err("a surface curve or a 3D offset curve has no form in a parameter space".into())
    );
    let mut b = bspline_curve(None, vec![[0.0; 3]; 4]);
    b.degree = usize::MAX;
    assert_eq!(
        curve_err(&Curve::BSpline(Box::new(b))),
        Err(format!("{} does not fit a Part 21 integer", usize::MAX))
    );
    assert_eq!(
        w.add("lower_case", &[]),
        Err("`lower_case` is not a Part 21 keyword".into())
    );
    assert_eq!(
        w.add("X", &[Param::Enum(".T".into())]),
        Err("`.T` is not a Part 21 keyword".into())
    );
    assert_eq!(
        w.add("X", &[Param::Typed("1A".into(), vec![])]),
        Err("`1A` is not a Part 21 keyword".into())
    );
    assert_eq!(
        w.add("X", &[Param::Binary("0G".into())]),
        Err("binary `0G` is not hexadecimal".into())
    );
    assert_eq!(
        w.add_complex(&[("", vec![])]),
        Err("`` is not a Part 21 keyword".into())
    );
    // A refused entity leaves nothing behind but the entities it refers to.
    assert!(!w.data().contains("lower_case") && !w.data().contains("=X("));
}

#[test]
fn every_parameter_kind_is_written() {
    let mut w = Writer::new();
    let id = w
        .add(
            "X",
            &[
                Param::Int(-7),
                Param::Real(1.5),
                Param::Str("s".into()),
                Param::Enum("UNSPECIFIED".into()),
                Param::Ref(1),
                Param::List(vec![]),
                Param::Typed("LENGTH_MEASURE".into(), vec![Param::Real(25.4)]),
                Param::Binary("0F".into()),
                Param::Null,
                Param::Derived,
            ],
        )
        .unwrap();
    assert_eq!(
        w.data(),
        "#1=X(-7,1.5,'s',.UNSPECIFIED.,#1,(),LENGTH_MEASURE(25.4),\"0F\",$,*);\n"
    );
    let ex = parse(&w);
    assert_eq!(
        ex.get(id).unwrap().params(),
        &[
            Param::Int(-7),
            Param::Real(1.5),
            Param::Str("s".into()),
            Param::Enum("UNSPECIFIED".into()),
            Param::Ref(1),
            Param::List(vec![]),
            Param::Typed("LENGTH_MEASURE".into(), vec![Param::Real(25.4)]),
            Param::Binary("0F".into()),
            Param::Null,
            Param::Derived,
        ]
    );
}
