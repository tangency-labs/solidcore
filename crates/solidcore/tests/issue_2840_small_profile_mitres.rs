//! Small circular profiles must produce valid mitred elbows at every tested scale.
#![allow(clippy::unwrap_used)]

use solidcore::curve::curve::Curve;
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::geom::axis2_placement::Axis2Placement;
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::math::Vec3;
use solidcore::nurbs::b_spline_curve::p1::{BSplineCurve, NurbCircleParam};
use solidcore::primitive::primitive_creation::{model::PrimitiveModel, p1::PrimitiveCreation};
use solidcore::topology::brep::validation::{check_brep, CheckLevel};
use solidcore::{Brep, Kernel};

fn elbow(k: &mut Kernel, radius: f64, leg: f64, tol: f64) -> Brep {
    let profile = Curve::from(
        BSplineCurve::create_circle_segment(
            3,
            &Axis2Placement::new(Vec3::ZERO, Vec3::X, Vec3::Y),
            radius,
            0.,
            360.,
            NurbCircleParam::Quadratic,
        )
        .unwrap(),
    );
    let bend = Vec3::new(0., 0., leg);
    let line = |a, b| Curve::from(BSplineCurve::create_line_segment(3, a, b).unwrap());
    let paths = [line(Vec3::ZERO, bend), line(bend, Vec3::new(leg, 0., leg))];
    let (_, mut body) = PrimitiveModel::build_new(k, |_, model| {
        let id = PrimitiveCreation::create_sweep_along_planar_path(
            model,
            PrimitiveModel::MAIN,
            &[profile],
            &paths,
            tol,
            true,
            true,
            true,
        )?;
        model.take_brep(id)
    })
    .unwrap();
    body.refresh_cached_topology();
    body
}

fn check_elbow(radius: f64, leg: f64, tol: f64) {
    let mut k = Kernel::default();
    let mut body = elbow(&mut k, radius, leg, tol);
    let check = check_brep(&mut k, &mut body, CheckLevel::Level2).unwrap();
    assert!(
        check.is_valid(),
        "radius={radius}, leg={leg}, tol={tol}: {:?}",
        check.lines()
    );
    let store = AttributeStore::default();
    let volume = brep_volume(
        &k,
        &mut body,
        1e-12,
        &PropertyAttributes {
            store: &store,
            key_of: &|_| ObjectKey(0),
        },
    )
    .unwrap();
    let exact = std::f64::consts::PI * radius * radius * 2. * leg;
    assert!((volume - exact).abs() < exact * 1e-4, "{volume} != {exact}");
}

#[test]
fn small_radius_short_legs() {
    check_elbow(1e-3, 1e-2, 1e-7);
}

#[test]
fn small_radius_long_legs() {
    check_elbow(1e-3, 1., 1e-7);
}

#[test]
fn intermediate_radius() {
    check_elbow(1e-2, 1e-1, 1e-5);
}

#[test]
fn larger_radius_control() {
    check_elbow(1e-1, 1., 1e-7);
}

/// Opt-in diagnostic output for the before/after render; normal tests write no files.
#[test]
#[ignore]
fn export_elbow_for_diagnostics() {
    use solidcore::topology::api_files::write_brep_to_file;
    use solidcore::topology::brep_data_write::FromBrepEnv;
    use solidcore::topology::topology::NotifyRecords;
    let out = std::env::var_os("SOLIDCORE_ELBOW_ARTIFACTS").unwrap();
    std::fs::create_dir_all(&out).unwrap();
    for (name, radius, leg) in [("small", 1e-3, 1e-2), ("control", 1e-1, 1.)] {
        let mut k = Kernel::default();
        let mut body = elbow(&mut k, radius, leg, 1e-7);
        for vertex in solidcore::topology::brep::p1::vertices(&body, None).unwrap() {
            let point = body.vertex(vertex).unwrap().point;
            let gap = solidcore::topology::vertex::corner::is_within_xsect_tol_3d_of_faces(
                &k, &mut body, vertex, true,
            )
            .unwrap();
            eprintln!("{name} vertex={vertex:?} point={point:?} gap={gap:?}");
        }
        let store = AttributeStore::default();
        let records = NotifyRecords::new();
        let env = FromBrepEnv {
            attributes: &store,
            key_of: &|e| body.owner_identity(e),
            records: &records,
            curve_key: &|_| ObjectKey(0),
            uv_curve_key: &|_| ObjectKey(0),
            surface_key: &|_| ObjectKey(0),
        };
        let bytes = write_brep_to_file(&k, Some(&body), &env, false).unwrap();
        std::fs::write(
            std::path::Path::new(&out).join(format!("{name}.smb")),
            bytes,
        )
        .unwrap();
    }
}
