//! tangency-labs/solidcore#2358 (the fork's smlib#124): the junction rig's boss, a cylinder
//! R6 standing 2 into a 60 x 60 x 10 plate, tilted about y, its root rounded, with
//! `UpstreamFixes::fillet_blend_corner` (WR-697).
//!
//! Tilted 20 degrees the cylinder's bottom cap rises out of the plate on the -x side, by
//! 6 sin 20 - 2 cos 20 = 0.173 at its top: a wedge of air under the cap. The root edge's
//! round ends at the two vertices where the root ellipse meets the cap's chord on the
//! plate, Ends of mixed convexity, on the cap: its rail on the wall must end on the rim,
//! which it does only while r <= r* = (6 sin 20 - 2 cos 20) / (1 - sin 20) = 0.262524174.
//! Past it the round's rail passes over the rim: a concave overflow, refused
//! (`SM_FILERR_NO_ROOM`) naming the largest radius whose rail still ends on the rim, less
//! the part's tolerance (1e-5): 0.262514174 (the fork's 54a4ee7, `BallTouchingAt` and
//! `RefuseEndOverflowingSideEdge`).
//!
//! Inputs and expectations are the fork's capi test `scenario_junction_boss_tilt_20_root_round`
//! (ce52279); the values are the fork at 18e96f1 (smlib#124 at 619c2f0) run on the same
//! inputs (`SmApiCreateBox`, `SmApiCreateCylinder` rotated and translated,
//! `SmApiBooleanUnion`, `SmFilletExecutive` with a circular section at constant radius,
//! continuity 2, keep features): at r 1.5 status 1001, `SM_FILERR_NO_ROOM`, fit
//! 0.262514174 and that reason; at 0.25 and 0.2625 built (11 faces, -0.576293097954 and
//! -0.635789172127 removed: the round adds material).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers

use solidcore::blend::fillet_executive::api::fillet_edges;
use solidcore::boolean::api_brep::boolean_union;
use solidcore::error::SolidError;
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::kernel::{Kernel, KernelOptions, UpstreamFixes};
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::math::{Point3, Vec3};
use solidcore::primitive::api_primitives::{create_box, create_cylinder};
use solidcore::topology::proto_topology::p1::api_transform::{
    api_rotate, api_translate, TransformObject,
};
use solidcore::topology::{edge, Brep, EdgeId};

/// Every opt-in fix on, with the blend corner as given, except the two the fork does not
/// have that move this body's measure: the exact trims of the cylinder's sections by the
/// plate (`cylinder_plane_section_trims_exact`) and the integrator's split at surface knots
/// (`precise_properties_split_at_surface_knots`). With them on the union measures 3.6e-5
/// less (38601.117929445747 against the fork's 38601.117965259655), and the removed
/// volumes move by 1.7e-4.
fn kernel(blend_corner: bool) -> Kernel {
    let mut k = Kernel::new(&KernelOptions::default());
    let mut f = UpstreamFixes::default();
    for &name in UpstreamFixes::NAMES {
        *f.flag_mut(name).unwrap() = true;
    }
    f.fillet_blend_corner = blend_corner;
    f.cylinder_plane_section_trims_exact = false;
    f.precise_properties_split_at_surface_knots = false;
    k.fixes = f;
    k
}

fn volume(k: &Kernel, b: &mut Brep) -> f64 {
    let store = AttributeStore::default();
    let attrs = PropertyAttributes {
        store: &store,
        key_of: &|_| ObjectKey(0),
    };
    brep_volume(k, b, 1e-10, &attrs).unwrap()
}

/// The plate with the boss tilted `deg` about y, its foot 2 into the plate.
fn boss(k: &mut Kernel, deg: f64) -> Brep {
    let a = deg.to_radians();
    let plate = create_box(k, Point3::new(-30.0, -30.0, -10.0), 60.0, 60.0, 10.0).unwrap();
    let mut cyl = create_cylinder(k, Point3::new(0.0, 0.0, 0.0), 6.0, 25.0).unwrap();
    let origin = Point3::new(0.0, 0.0, 0.0);
    api_rotate(
        k,
        Some(TransformObject::Brep(&mut cyl)),
        origin,
        Vec3::Y,
        deg,
    )
    .unwrap();
    api_translate(
        k,
        Some(TransformObject::Brep(&mut cyl)),
        Vec3::new(-2.0 * a.sin(), 0.0, -2.0 * a.cos()),
    )
    .unwrap();
    boolean_union(k, plate, cyl).unwrap()
}

/// The root ellipse's edges: on the plate's top, their middles on the cylinder.
fn root(k: &Kernel, b: &Brep, deg: f64) -> Vec<EdgeId> {
    let a = deg.to_radians();
    let axis = Vec3::new(a.sin(), 0.0, a.cos());
    let foot = Point3::new(-2.0 * a.sin(), 0.0, -2.0 * a.cos());
    let mut out = Vec::new();
    for &e in b.edge_list() {
        let Some(c) = edge::curve(b, e).unwrap() else {
            continue;
        };
        let iv = edge::interval(b, e).unwrap();
        let p0 = c.evaluate_point(k, iv.min()).unwrap();
        let p1 = c.evaluate_point(k, iv.max()).unwrap();
        let pm = c.evaluate_point(k, 0.5 * (iv.min() + iv.max())).unwrap();
        if p0.z.abs() > 1e-6 || p1.z.abs() > 1e-6 || pm.z.abs() > 1e-6 {
            continue;
        }
        let d = pm - foot;
        let q = d - axis * d.dot(axis);
        if (q.length() - 6.0).abs() < 1e-6 {
            out.push(e);
        }
    }
    out
}

/// The root of the boss tilted `deg` rounded `r`: the volume it removes and the face count.
fn rounded(blend_corner: bool, deg: f64, r: f64) -> Result<(f64, usize), SolidError> {
    let mut k = kernel(blend_corner);
    let mut b = boss(&mut k, deg);
    let before = volume(&k, &mut b);
    let es = root(&k, &b, deg);
    assert_eq!(es.len(), 2);
    fillet_edges(&mut k, &mut b, &es, r, 1, 2, 1.0)?;
    Ok((before - volume(&k, &mut b), b.face_list().len()))
}

/// The bound r*, closed form.
fn bound() -> f64 {
    let a = 20.0_f64.to_radians();
    (6.0 * a.sin() - 2.0 * a.cos()) / (1.0 - a.sin())
}

#[test]
fn the_bound_is_the_closed_form() {
    assert!((bound() - 0.262_524_174).abs() <= 1e-9, "{}", bound());
}

/// Tilted 20 degrees and rounded 1.5: refused, naming r* less the part's tolerance.
#[test]
fn tilt_20_is_refused_naming_the_largest_radius() {
    let err = rounded(true, 20.0, 1.5).unwrap_err();
    let SolidError::FilletRefused { code, fit, .. } = err else {
        panic!("{err:?}");
    };
    assert_eq!(code, 1001);
    assert!((fit.size - 0.262_514_174).abs() <= 1e-9, "{}", fit.size);
    assert!(fit.size < bound() && fit.size >= bound() - 0.02);
}

/// The radii below the bound build, as the fork's: the removed volumes agree to 2.5e-9
/// (-0.576293096943 and -0.635789174601 here; the union to 4e-8 of 38601.12).
#[test]
fn tilt_20_below_the_bound_builds() {
    for (r, removed) in [(0.25, -0.576_293_097_954), (0.2625, -0.635_789_172_127)] {
        let (got, faces) = rounded(true, 20.0, r).unwrap();
        assert!((got - removed).abs() <= 5e-9, "{r}: {got}");
        assert_eq!(faces, 11);
    }
}
