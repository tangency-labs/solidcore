//! Issue #1999: shallow crease fillets return success with invalid geometry.
//! The opt-in repair must construct valid rounds and preserve already valid results.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers

use std::f64::consts::PI;

use solidcore::blend::fillet_executive::api::fillet_edges;
use solidcore::curve::curve::Curve;
use solidcore::math::Vec3;
use solidcore::primitive::api_curves::create_line_segment;
use solidcore::sweep::api_brep::create_taper_extrude;
use solidcore::topology::brep::p2::is_manifold_solid;
use solidcore::topology::edge;
use solidcore::topology::EdgeId;
use solidcore::{Brep, Kernel};

fn line(a: (f64, f64), b: (f64, f64)) -> Curve {
    Curve::from(create_line_segment(Vec3::new(a.0, a.1, 0.0), Vec3::new(b.0, b.1, 0.0)).unwrap())
}

/// The rise of the crease above y = 1 for a turn of `turn` degrees.
fn rise(turn: f64) -> f64 {
    3.0 * (0.5 * turn * PI / 180.0).tan()
}

/// The 6 x 1 roof profile, extruded 2 along z.
fn block(k: &mut Kernel, turn: f64) -> Brep {
    block_scaled(k, turn, 1.0)
}

fn block_scaled(k: &mut Kernel, turn: f64, scale: f64) -> Brep {
    let line =
        |a: (f64, f64), b: (f64, f64)| line((a.0 * scale, a.1 * scale), (b.0 * scale, b.1 * scale));
    let h = rise(turn);
    let curves = [
        line((0.0, 0.0), (6.0, 0.0)),
        line((6.0, 0.0), (6.0, 1.0)),
        line((6.0, 1.0), (3.0, 1.0 + h)),
        line((3.0, 1.0 + h), (0.0, 1.0)),
        line((0.0, 1.0), (0.0, 0.0)),
    ];
    let refs: Vec<Option<&Curve>> = curves.iter().map(Some).collect();
    create_taper_extrude(k, &refs, 2.0 * scale, 0.0, 3, false).unwrap()
}

/// The crease: the edge along z at (3, 1 + rise).
fn crease(k: &Kernel, b: &Brep, turn: f64) -> EdgeId {
    crease_scaled(k, b, turn, 1.0)
}

fn crease_scaled(k: &Kernel, b: &Brep, turn: f64, scale: f64) -> EdgeId {
    let h = rise(turn);
    let on = |p: Vec3| {
        (p.x - 3.0 * scale).abs() < 1e-9 * scale && (p.y - (1.0 + h) * scale).abs() < 1e-9 * scale
    };
    b.edge_list()
        .iter()
        .copied()
        .find(|&e| {
            let curve = edge::curve(b, e).unwrap().unwrap();
            let ivl = edge::interval(b, e).unwrap();
            let p = curve.evaluate_point(k, ivl.min()).unwrap();
            let q = curve.evaluate_point(k, ivl.max()).unwrap();
            on(p) && on(q) && (p.z - q.z).abs() > scale
        })
        .expect("the crease edge")
}

/// The debug dump of `b` without the parts a copy never carries: notification
/// identities (issued to the objects themselves) and surface cache slots (caches only
/// speed up queries).
fn dump(b: &impl std::fmt::Debug) -> String {
    let s = format!("{b:?}");
    let mut out = String::new();
    let mut rest = s.as_str();
    loop {
        let next = ["IdentitySlot(", "CacheSlot("]
            .iter()
            .filter_map(|key| rest.find(key).map(|at| (at, key.len())))
            .min();
        let Some((at, len)) = next else {
            out.push_str(rest);
            return out;
        };
        out.push_str(&rest[..at]);
        out.push('_');
        let mut depth = 1;
        let mut end = at + len;
        for (i, c) in rest[at + len..].char_indices() {
            depth += match c {
                '(' => 1,
                ')' => -1,
                _ => 0,
            };
            if depth == 0 {
                end = at + len + i + 1;
                break;
            }
        }
        rest = &rest[end..];
    }
}

fn valid(k: &mut Kernel, b: &mut Brep) {
    use solidcore::topology::brep::validation::{check_brep, CheckLevel};
    let check = check_brep(k, b, CheckLevel::Level2).unwrap();
    assert_eq!(check.pointers_status, 1000, "{check:?}");
    assert!(check.is_valid(), "{:?}", check.lines());
}

#[test]
fn default_reproduces_the_native_shallow_crease_gap() {
    use solidcore::topology::brep::validation::{check_brep, CheckLevel};
    let mut k = Kernel::default();
    assert!(!k.fixes.fillet_repair_split_vertices);
    let mut b = block(&mut k, 2.5);
    valid(&mut k, &mut b);
    let e = crease(&k, &b, 2.5);
    fillet_edges(&mut k, &mut b, &[e], 0.12, 1, 1, 1.0).unwrap();
    assert!(is_manifold_solid(&b).unwrap());
    assert_eq!(b.face_list().len(), 8);
    let check = check_brep(&mut k, &mut b, CheckLevel::Level2).unwrap();
    assert_eq!(check.assert_valid, Some(false));
    let gaps: Vec<_> = check
        .failures()
        .filter(|r| r.reporting_type == 16021 && r.test_index == 1)
        .collect();
    assert_eq!(gaps.len(), 2);
    for gap in gaps {
        assert!((gap.value - 2.8563548e-5).abs() < 2e-12, "{gap:?}");
        assert_eq!(gap.tol3d, 2e-5);
    }
}

#[test]
fn misplaced_split_vertices_are_repaired_successfully() {
    for (half_degree, turn, radius) in [
        (false, 2.5, 0.12),
        (false, 3.0, 0.12),
        (true, 0.6, 2.0),
        (true, 1.5, 0.5),
    ] {
        for manifold in [false, true] {
            let mut k = Kernel::default();
            k.fixes.fillet_half_degree_crease = half_degree;
            k.fixes.fillet_repair_split_vertices = true;
            k.fixes.fillet_refuse_non_manifold_result = manifold;
            let mut b = block(&mut k, turn);
            valid(&mut k, &mut b);
            let e = crease(&k, &b, turn);
            fillet_edges(&mut k, &mut b, &[e], radius, 1, 1, 1.0)
                .unwrap_or_else(|e| panic!("turn={turn} radius={radius}: {e:?}"));
            valid(&mut k, &mut b);
            assert!(is_manifold_solid(&b).unwrap());
            assert_eq!(b.face_list().len(), 8);
        }
    }
}

#[test]
fn valid_rounds_are_identical_with_the_check_enabled() {
    for (half_degree, turn, radius) in [
        (false, 2.5, 2.0),
        (false, 3.0, 2.0),
        (true, 1.0, 2.0),
        (true, 1.0, 20.0),
    ] {
        let mut outputs = Vec::new();
        for enabled in [false, true] {
            let mut k = Kernel::default();
            k.fixes.fillet_half_degree_crease = half_degree;
            k.fixes.fillet_repair_split_vertices = enabled;
            let mut b = block(&mut k, turn);
            let e = crease(&k, &b, turn);
            fillet_edges(&mut k, &mut b, &[e], radius, 1, 1, 1.0).unwrap();
            outputs.push(dump(&b));
            valid(&mut k, &mut b);
        }
        assert_eq!(outputs[0], outputs[1], "turn={turn} radius={radius}");
    }
}

#[test]
fn all_three_edge_fillet_wrappers_validate_the_final_result() {
    use solidcore::blend::fillet_executive::api::{fillet_edges_per_edge, variable_radius_fillet};
    use solidcore::blend::fillet_solver::p1::FilletSurfaceGeneratorType;
    for per_edge in [false, true] {
        for radius in [0.12, 2.0] {
            let mut k = Kernel::default();
            k.fixes.fillet_repair_split_vertices = true;
            let mut b = block(&mut k, 2.5);
            let e = crease(&k, &b, 2.5);
            let result = if per_edge {
                fillet_edges_per_edge(
                    &mut k,
                    &mut b,
                    &[e],
                    &[radius],
                    &[FilletSurfaceGeneratorType::Circular],
                    1,
                    1.0,
                )
            } else {
                variable_radius_fillet(&mut k, &mut b, &[e], radius, radius, 1, 1)
            };
            result.unwrap_or_else(|e| panic!("per_edge={per_edge} radius={radius}: {e:?}"));
            valid(&mut k, &mut b);
        }
    }
}

#[test]
fn multi_edge_corner_is_refused_even_above_the_single_crease_gap_boundary() {
    use solidcore::topology::brep::validation::{check_brep, CheckLevel};
    for enabled in [false, true] {
        let mut k = Kernel::default();
        k.fixes.fillet_repair_split_vertices = enabled;
        let mut b = block(&mut k, 2.5);
        let ridge = Vec3::new(3.0, 1.0 + rise(2.5), 0.0);
        let edges: Vec<_> = b
            .edge_list()
            .iter()
            .copied()
            .filter(|&e| {
                let c = edge::curve(&b, e).unwrap().unwrap();
                let iv = edge::interval(&b, e).unwrap();
                [iv.min(), iv.max()]
                    .into_iter()
                    .any(|t| c.evaluate_point(&k, t).unwrap().distance(ridge) < 1e-9)
            })
            .collect();
        assert_eq!(edges.len(), 3);
        let before = dump(&b);
        let result = fillet_edges(&mut k, &mut b, &edges, 1.0, 1, 1, 1.0);
        if enabled {
            assert!(
                matches!(
                    result,
                    Err(solidcore::SolidError::Upstream {
                        operation: "blend::fillet_executive::api::invalid_result",
                        code: 1001,
                        entity: None
                    })
                ),
                "{result:?}"
            );
            assert_eq!(dump(&b), before);
            valid(&mut k, &mut b);
        } else {
            result.unwrap();
            assert!(!check_brep(&mut k, &mut b, CheckLevel::Level2)
                .unwrap()
                .is_valid());
        }
    }
}

/// Optional artifact for independent validation with the unmodified C++ reference.
#[test]
#[ignore = "exports repaired models to SOLIDCORE_1999_EXPORT_DIR"]
fn export_repaired_models_for_native_validation() {
    use solidcore::geom::attribute::{AttributeStore, ObjectKey};
    use solidcore::topology::brep_data_write::FromBrepEnv;
    use solidcore::topology::topology::NotifyRecords;
    let directory = std::path::PathBuf::from(std::env::var("SOLIDCORE_1999_EXPORT_DIR").unwrap());
    for (turn, radius) in [(2.5, 0.12), (0.6, 2.0), (1.5, 0.5)] {
        let mut k = Kernel::default();
        k.fixes.fillet_repair_split_vertices = true;
        k.fixes.fillet_half_degree_crease = true;
        let mut b = block(&mut k, turn);
        let e = crease(&k, &b, turn);
        fillet_edges(&mut k, &mut b, &[e], radius, 1, 1, 1.0).unwrap();
        valid(&mut k, &mut b);
        let attrs = AttributeStore::new();
        let records = NotifyRecords::new();
        let own = |e| b.owner_identity(e);
        let env = FromBrepEnv {
            attributes: &attrs,
            key_of: &own,
            records: &records,
            curve_key: &|_| ObjectKey(0),
            uv_curve_key: &|_| ObjectKey(0),
            surface_key: &|_| ObjectKey(0),
        };
        let bytes =
            solidcore::topology::api_files::write_brep_to_file(&k, Some(&b), &env, true).unwrap();
        std::fs::write(directory.join(format!("crease-{turn}-{radius}.smb")), bytes).unwrap();
    }
}

#[test]
fn repair_scales_with_model_tolerance() {
    for scale in [0.1, 10.0, 1000.0] {
        let mut k = Kernel::default();
        k.tolerances.zone_tol_3d *= scale;
        k.fixes.fillet_repair_split_vertices = true;
        let mut b = block_scaled(&mut k, 2.5, scale);
        valid(&mut k, &mut b);
        let e = crease_scaled(&k, &b, 2.5, scale);
        fillet_edges(&mut k, &mut b, &[e], 0.12 * scale, 1, 1, 1.0)
            .unwrap_or_else(|e| panic!("scale={scale}: {e:?}"));
        valid(&mut k, &mut b);
        assert!(is_manifold_solid(&b).unwrap());
    }
}

#[test]
fn repair_preserves_curves_topology_and_tolerances() {
    let run = |repair| {
        let mut k = Kernel::default();
        k.fixes.fillet_repair_split_vertices = repair;
        let mut b = block(&mut k, 2.5);
        let e = crease(&k, &b, 2.5);
        fillet_edges(&mut k, &mut b, &[e], 0.12, 1, 1, 1.0).unwrap();
        b
    };
    let original = run(false);
    let repaired = run(true);
    assert_eq!(original.tolerance(), repaired.tolerance());
    assert_eq!(original.face_list(), repaired.face_list());
    assert_eq!(original.edge_list(), repaired.edge_list());
    for (id, v) in original.vertices().iter() {
        let fixed = repaired.vertex(id).unwrap();
        assert_eq!(v.tolerance, fixed.tolerance);
        assert_eq!(v.vertexuses, fixed.vertexuses);
    }
    for (id, e) in original.edges().iter() {
        let fixed = repaired.edge(id).unwrap();
        assert_eq!(dump(e), dump(fixed), "edge {id:?}");
    }
    for (id, f) in original.faces().iter() {
        let fixed = repaired.face(id).unwrap();
        assert_eq!(dump(f), dump(fixed), "face {id:?}");
    }
    let moved = original
        .vertices()
        .iter()
        .filter(|(id, v)| v.point != repaired.vertex(*id).unwrap().point)
        .count();
    assert_eq!(moved, 2);
}
