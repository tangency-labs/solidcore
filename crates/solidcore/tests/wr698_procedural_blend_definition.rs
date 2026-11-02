//! WR-698 (#2273, the fork's 014f66f, c6d9b9a and 0136260): a procedural blend evaluates
//! as its definition, a station solve and the section exact in t, with the kernel's
//! station memo; its file record is written and read back (type 41 only with the flags).
//!
//! The geometry: a ball of radius 2 between the plane y = 0 (side -1) and the plane
//! z = 0 (side +1), or a biquadratic bump over z = 0 (centre control point 0.6 high),
//! along the leader (0, 2, 2) -> (10, 2, 2). The values were checked against the fork
//! (tangency-labs/smlib@0136260, `SMLib` and `SM_API` built locally with the harness flags)
//! through `SmProceduralBlend::Create`, `EvaluateDefinition`, `GetCounts`, `WriteToDB` and
//! `ReadFromDB` on the same inputs: every partial (3564 values over 11 definitions, 12
//! points each, to order 2 in s and t), the counts and the binary records agreed. The
//! values quoted below are the fork's. The circle's points are also checked against the
//! closed form.

#![allow(clippy::unwrap_used)]

use solidcore::blend::fillet_section::section::{GnSectionFsg, GN_THUMBWEIGHT_EQUAL_AREA};
use solidcore::geom::vector3d::{Point3, Vec3};
use solidcore::kernel::Kernel;
use solidcore::math::Point2;
use solidcore::nurbs::b_spline_curve::{BSplineCurve, KnotType};
use solidcore::nurbs::b_spline_surface::{
    BSplineSurface, BSplineSurfaceForm, CanonicalSurfaceData,
};
use solidcore::surface::procedural_blend::{
    BlendCounts, BlendDefinitionInput, ProceduralBlend, PROCEDURAL_BLEND_TYPE,
};
use solidcore::surface::surface::Surface;
use solidcore::topology::database_io::{DatabaseReader, DatabaseWriter, FileType};
use solidcore::topology::geometry_write::write_surface;

const L: f64 = 10.0;
const R: f64 = 2.0;

fn kernel() -> Kernel {
    let mut k = Kernel::default();
    k.fixes.fillet_gn_section_generator = true;
    k.fixes.fillet_round_as_definition = true;
    k
}

fn bilinear(p: [Point3; 4]) -> BSplineSurface {
    BSplineSurface::create_bilinear_surface(p[0], p[1], p[2], p[3]).unwrap()
}

fn bump(h: f64) -> BSplineSurface {
    let mut cps = Vec::new();
    for i in 0..3 {
        for j in 0..3 {
            let z = if i == 1 && j == 1 { h } else { 0.0 };
            cps.push(Point3::new(
                f64::from(i) * L / 2.0,
                f64::from(j) * L / 2.0,
                z,
            ));
        }
    }
    BSplineSurface::create_canonical(&CanonicalSurfaceData {
        u_degree: 2,
        v_degree: 2,
        control_points: &cps,
        form: BSplineSurfaceForm::Unspecified,
        u_multiplicities: &[3, 3],
        v_multiplicities: &[3, 3],
        u_knots: &[0.0, 1.0],
        v_knots: &[0.0, 1.0],
        knot_type: KnotType::Unspecified,
        weights: None,
        uv_domain: None,
    })
    .unwrap()
}

fn line(a: Point3, b: Point3) -> BSplineCurve {
    BSplineCurve::create_line_segment(3, a, b).unwrap()
}

/// The definition: `bumpy` puts the bump in place of z = 0, `law` the radius
/// 1.5 + s (s over [0, 1]).
fn blend(k: &Kernel, bumpy: bool, order: u32, thumbweight: f64, law: bool) -> ProceduralBlend {
    let o = Point3::new(0.0, 0.0, 0.0);
    let a = Surface::from(bilinear([
        o,
        Point3::new(L, 0.0, 0.0),
        Point3::new(0.0, 0.0, L),
        Point3::new(L, 0.0, L),
    ]));
    let b = if bumpy {
        Surface::from(bump(0.6))
    } else {
        Surface::from(bilinear([
            o,
            Point3::new(L, 0.0, 0.0),
            Point3::new(0.0, L, 0.0),
            Point3::new(L, L, 0.0),
        ]))
    };
    let net = bilinear([
        o,
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
    ]);
    let leader = line(Point3::new(0.0, R, R), Point3::new(L, R, R));
    let seed = line(Point3::new(0.0, R / L, 0.0), Point3::new(1.0, R / L, 0.0));
    let radius = line(Point3::new(1.5, 0.0, 0.0), Point3::new(2.5, 0.0, 0.0));
    let input = BlendDefinitionInput {
        supports: [&a, &b],
        sides: [-1.0, 1.0],
        offsets: [R, R],
        law: law.then_some(&radius),
        order,
        thumbweight,
        leader: &leader,
        seeds: [&seed, &seed],
        bound: 1e-3,
    };
    ProceduralBlend::create(k, &net, &input, |n, w| {
        GnSectionFsg::is_valid(n as usize, w)
    })
    .unwrap()
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9 * 1f64.max(a.abs()).max(b.abs())
}

fn assert_vec(got: Vec3, want: [f64; 3], what: &str) {
    assert!(
        close(got.x, want[0]) && close(got.y, want[1]) && close(got.z, want[2]),
        "{what}: {got:?} != {want:?}"
    );
}

/// At (s, t) = (0.37, 0.6): the point, d^2 S / ds dt and d^4 S / ds^2 dt^2, the fork's.
#[test]
fn definitions_evaluate_as_the_fork_s() {
    let k = kernel();
    let ea = GN_THUMBWEIGHT_EQUAL_AREA;
    #[allow(clippy::type_complexity)]
    let cases: [(&str, bool, u32, f64, bool, [[f64; 3]; 3]); 7] = [
        (
            "plane n 2 w 0",
            false,
            2,
            0.0,
            false,
            [
                [3.6999999999999988, 0.92571428571428493, 0.18285714285714266],
                [0.0; 3],
                [0.0; 3],
            ],
        ),
        (
            "plane n 3 w -1",
            false,
            3,
            -1.0,
            false,
            [
                [3.6999999999999997, 1.2525773195876284, 0.18851251840942546],
                [0.0; 3],
                [0.0; 3],
            ],
        ),
        (
            "plane n 2 equal area",
            false,
            2,
            ea,
            false,
            [
                [3.6999999999999993, 1.0735924639518113, 0.23047636434570989],
                [0.0; 3],
                [0.0; 3],
            ],
        ),
        (
            "bump n 2 equal area",
            true,
            2,
            ea,
            false,
            [
                [3.7081710618578949, 1.1216620061581504, 0.29794471822651397],
                [
                    -0.33027757576901201,
                    0.16993205953095944,
                    0.05794413006238934,
                ],
                [0.13167003776583264, 6.7676147994552807, -3.8411749725938713],
            ],
        ),
        (
            "bump n 3 w 0",
            true,
            3,
            0.0,
            false,
            [
                [3.7067830252568141, 1.0441037518960723, 0.19005330121057595],
                [
                    -0.30262813328989235,
                    0.17384926344700599,
                    0.11123189916022572,
                ],
                [
                    0.098355439055517191,
                    8.6358244659062446,
                    -2.5725191740645639,
                ],
            ],
        ),
        (
            "bump n 1",
            true,
            1,
            0.0,
            false,
            [
                [3.7166223003140217, 0.87103605536466744, 0.43990796265579407],
                [
                    -0.14685906050567388,
                    0.10685015971002275,
                    0.012807829671101237,
                ],
                [
                    0.09921158399671709,
                    -0.14682964597773029,
                    -2.3572242028696042,
                ],
            ],
        ),
        (
            "bump n 1 law",
            true,
            1,
            0.0,
            true,
            [
                [3.7927539614653258, 0.81750870003218679, 0.41335652579648763],
                [
                    -0.16045050294581342,
                    1.4665607799716496,
                    -0.93913899925185074,
                ],
                [1.0741110766659812, -0.11777436551962284, -2.195182493991612],
            ],
        ),
    ];
    for (name, bumpy, n, w, law, want) in cases {
        let pb = blend(&k, bumpy, n, w, law);
        let d = pb.evaluate_definition(&k, 0.37, 0.6, 2, 2).unwrap();
        assert_vec(d[0][0], want[0], name);
        // (the fork's d[4] and d[8] beyond the point are compared where they are not
        // rounding noise: the planes' are zero to 1e-14)
        if bumpy {
            assert_vec(d[1][1], want[1], name);
            assert_vec(d[2][2], want[2], name);
        } else {
            assert!(
                d[1][1].length() < 1e-12 && d[2][2].length() < 1e-12,
                "{name}"
            );
        }
    }
}

/// The circle between the planes: about (x, 2, 2), radius 2, from the foot on y = 0
/// (t = 0) to the foot on z = 0 (t = 1), x = 10 s (closed form).
#[test]
fn the_circle_between_planes_is_the_ball_s() {
    let k = kernel();
    let pb = blend(&k, false, 1, 0.0, false);
    for s in [0.0, 0.37, 1.0] {
        for t in [0.0, 0.25, 0.6, 1.0] {
            let d = pb.evaluate_definition(&k, s, t, 1, 1).unwrap();
            let p = d[0][0];
            assert!((p.x - L * s).abs() < 1e-12, "{p:?}");
            let r = ((p.y - R).powi(2) + (p.z - R).powi(2)).sqrt();
            assert!((r - R).abs() < 1e-12, "{s} {t}: {r}");
            assert!((d[1][0] - Vec3::new(L, 0.0, 0.0)).length() < 1e-12);
        }
        let d0 = pb.evaluate_definition(&k, s, 0.0, 0, 0).unwrap();
        let d1 = pb.evaluate_definition(&k, s, 1.0, 0, 0).unwrap();
        assert!((d0[0][0] - Vec3::new(L * s, 0.0, R)).length() < 1e-12);
        assert!((d1[0][0] - Vec3::new(L * s, R, 0.0)).length() < 1e-12);
    }
}

/// The station memo: points at one s solve once (the fork's counts: 12 points at 3 s
/// values, 3 solves; on the bump 9 Newton steps, at most 3 in one solve).
#[test]
fn a_row_of_points_solves_its_station_once() {
    let k = kernel();
    let pb = blend(&k, true, 2, GN_THUMBWEIGHT_EQUAL_AREA, false);
    k.blend_counts.set(BlendCounts::default());
    let mut first = Vec::new();
    for s in [0.0, 0.37, 1.0] {
        for t in [0.0, 0.25, 0.6, 1.0] {
            first.push(pb.evaluate_definition(&k, s, t, 2, 2).unwrap());
        }
    }
    let c = k.blend_counts.get();
    assert_eq!(
        (c.evaluations, c.solves, c.newton_steps, c.most_newton_steps),
        (12, 3, 9, 3)
    );
    // a hit gives what a fresh solve gives: a copy (a new serial) solves afresh
    let copy = pb.copy(&k);
    let mut i = 0;
    for s in [0.0, 0.37, 1.0] {
        for t in [0.0, 0.25, 0.6, 1.0] {
            let d = copy.evaluate_definition(&k, s, t, 2, 2).unwrap();
            for (a, b) in d.iter().flatten().zip(first[i].iter().flatten()) {
                assert_eq!(a, b);
            }
            i += 1;
        }
    }
    assert_eq!(k.blend_counts.get().solves, 6);
}

fn fnv1a(bytes: &[u8]) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325_u64;
    for &b in bytes {
        h = (h ^ u64::from(b)).wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

/// The binary record is the fork's (its length and FNV-1a hash), and reads back.
#[test]
fn the_record_is_the_fork_s_and_reads_back() {
    let k = kernel();
    for (bumpy, n, w, law, len, hash) in [
        (false, 1, 0.0, false, 988, 0xa632_9df0_1772_45bd_u64),
        (
            true,
            2,
            GN_THUMBWEIGHT_EQUAL_AREA,
            false,
            1140,
            0x59b3_c770_4a8d_1f3a,
        ),
        (true, 1, 0.0, true, 1250, 0x695b_ae87_4d74_044a),
    ] {
        let surface = Surface::from(blend(&k, bumpy, n, w, law));
        let want = surface.evaluate_point(&k, Point2::new(0.37, 0.6)).unwrap();
        for ft in [FileType::Binary, FileType::Ascii] {
            let mut db = DatabaseWriter::open_for_write(ft);
            write_surface(&surface, &mut db, 30).unwrap();
            let bytes = db.into_bytes();
            if ft == FileType::Binary {
                assert_eq!((bytes.len(), fnv1a(&bytes)), (len, hash));
            }
            let mut r = DatabaseReader::new();
            r.open_for_read(bytes.clone(), ft);
            let back = Surface::read_from_db(&k, PROCEDURAL_BLEND_TYPE, &mut r, 30).unwrap();
            let p = back.evaluate_point(&k, Point2::new(0.37, 0.6)).unwrap();
            // binary is exact; ASCII rounds the doubles it writes
            let tol = if ft == FileType::Binary { 0.0 } else { 1e-9 };
            assert!((p - want).length() <= tol, "{ft:?}: {p:?} != {want:?}");

            // without the flag (or with only one of the two) the type is refused, as the
            // pinned reader refuses it
            for (gn, def) in [(false, false), (true, false), (false, true)] {
                let mut off = Kernel::default();
                off.fixes.fillet_gn_section_generator = gn;
                off.fixes.fillet_round_as_definition = def;
                let mut r = DatabaseReader::new();
                r.open_for_read(bytes.clone(), ft);
                assert!(Surface::read_from_db(&off, PROCEDURAL_BLEND_TYPE, &mut r, 30).is_err());
            }
        }
    }
}

/// A hostile record is refused: a singular map, a number not finite, a negative order.
#[test]
fn hostile_records_are_refused() {
    let k = kernel();
    let surface = Surface::from(blend(&k, false, 2, 0.0, false));
    let mut db = DatabaseWriter::open_for_write(FileType::Binary);
    write_surface(&surface, &mut db, 30).unwrap();
    let good = db.into_bytes();
    // the record opens with the order (8 bytes), two booleans (4 each), then the ten
    // doubles: sides, offsets, thumbweight, bound, a0, b0, a1, b1
    let at = |i: usize| 16 + 8 * i;
    for (offset, value) in [
        (at(7), 0.0_f64),
        (at(5), f64::NAN),
        (at(2), f64::INFINITY),
        (at(0), 0.5),
    ] {
        let mut bad = good.clone();
        bad[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
        let mut r = DatabaseReader::new();
        r.open_for_read(bad, FileType::Binary);
        assert!(Surface::read_from_db(&k, PROCEDURAL_BLEND_TYPE, &mut r, 30).is_err());
    }
    let mut bad = good.clone();
    bad[0..8].copy_from_slice(&(-1_i64).to_le_bytes());
    let mut r = DatabaseReader::new();
    r.open_for_read(bad, FileType::Binary);
    assert!(Surface::read_from_db(&k, PROCEDURAL_BLEND_TYPE, &mut r, 30).is_err());
    // a truncated record
    let mut r = DatabaseReader::new();
    r.open_for_read(good[..good.len() / 2].to_vec(), FileType::Binary);
    assert!(Surface::read_from_db(&k, PROCEDURAL_BLEND_TYPE, &mut r, 30).is_err());
}

/// Partials past what the series hold are refused, not overrun.
#[test]
fn partials_past_the_series_are_refused() {
    let k = kernel();
    let pb = blend(&k, false, 4, 0.0, false);
    assert!(pb.evaluate_definition(&k, 0.5, 0.5, 8, 0).is_err());
    assert!(pb.evaluate_definition(&k, 0.5, 0.5, 0, 65).is_err());
    let ea = blend(&k, false, 2, GN_THUMBWEIGHT_EQUAL_AREA, false);
    assert!(ea.evaluate_definition(&k, 0.5, 0.5, 5, 0).is_err());
    assert!(ea.evaluate_definition(&k, 0.5, 0.5, 4, 3).is_ok());
}
