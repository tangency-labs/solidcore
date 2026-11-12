//! The searches a solid primitive runs for its infinite region, as a count.
//!
//! A box, a cylinder or a cone is made face by face and stitched. The stitch's last glue
//! closes the shell, which splits the B-rep's one region in two, and the glue then searches
//! for the infinite one of the two: it picks a ray to an edge and classifies it against the
//! shell (`FindInfiniteRegion`). The stitch ends with the same search over the whole B-rep,
//! and that second answer is the one the primitive keeps. [`Kernel::work_counters`] counts
//! each search (`infinite_region_searches`), the same on every machine.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers

use std::f64::consts::PI;

use solidcore::kernel::UpstreamFixes;
use solidcore::mass::integrator::properties::brep_volume_default;
use solidcore::math::Point3;
use solidcore::primitive::api_primitives::{create_box, create_cone, create_cylinder};
use solidcore::topology::brep::p1;
use solidcore::topology::Brep;
use solidcore::{Kernel, Result};

fn kernel(solidcraft: bool) -> Kernel {
    let mut fixes = UpstreamFixes::default();
    if solidcraft {
        fixes = UpstreamFixes::all();
        fixes.fillet_refuse_rounds_breaking_through = false;
        fixes.fillet_round_ends_on_outline = false;
    }
    Kernel {
        fixes,
        ..Kernel::default()
    }
}

/// Makes one primitive, checks it is the solid asked for (two regions, the infinite one
/// void and first on the list, the other not void, `faces` faces and the exact `volume`)
/// and returns the searches for the infinite region its making ran.
fn searches(
    solidcraft: bool,
    make: impl Fn(&mut Kernel) -> Result<Brep>,
    faces: usize,
    volume: f64,
) -> u64 {
    let mut k = kernel(solidcraft);
    k.work_counters.reset();
    let mut b = make(&mut k).unwrap();
    let n = k.work_counters.infinite_region_searches();
    let regions = p1::regions(&b, None).unwrap();
    assert_eq!(regions.len(), 2);
    let infinite = b.infinite_region().unwrap();
    assert_eq!(regions[0], infinite);
    assert_eq!(b.region(infinite).unwrap().is_void, 1);
    assert_eq!(b.region(regions[1]).unwrap().is_void, 0);
    assert_eq!(p1::faces(&b, None).unwrap().len(), faces);
    let measured = brep_volume_default(&k, &mut b, 1e-6).unwrap();
    assert!(
        (measured - volume).abs() <= 1e-6 * volume,
        "volume {measured}, exact {volume}"
    );
    n
}

/// A box or a cylinder made alone in its B-rep, on a unit frame, with sides, height and
/// radius beyond the stitch's reach and of a size its tolerance serves, is one closed
/// solid: the stitch can glue nothing but the edges its faces share, and nothing after the
/// glue that closes it, so the search at the stitch's end is the only one whose answer is
/// read. One search for each; with the search at the closing glue too there were two.
#[test]
fn a_primitive_searches_for_its_infinite_region_once() {
    for solidcraft in [true, false] {
        let p = Point3::new(10.0, 20.0, -1.0);
        let counts = (
            searches(solidcraft, |k| create_box(k, p, 8.0, 14.0, 3.0), 6, 336.0),
            searches(
                solidcraft,
                |k| create_cylinder(k, p, 4.0, 14.0),
                3,
                PI * 16.0 * 14.0,
            ),
        );
        assert_eq!(counts, (1, 1), "SolidCraft's fixes: {solidcraft}");
    }
}

/// Where the premise does not hold the primitive is stitched as every stitch is, with a
/// search at each shell a glue closes, two in all:
/// - a cylinder 5e-4 high, inside the stitch's reach of 8e-4, which lets the stitch glue
///   more than the shared edges (1e-7 high it comes out as three regions);
/// - a cylinder 600 from the origin, a size the kernel's table gives a tolerance of 1e-4
///   where the B-rep's is 1e-5;
/// - a cone, with an apex or without.
#[test]
fn a_thin_far_or_conic_primitive_searches_at_each_closure() {
    let p = Point3::new(10.0, 20.0, -1.0);
    let thin = searches(
        true,
        |k| create_cylinder(k, p, 4.0, 5e-4),
        3,
        PI * 16.0 * 5e-4,
    );
    let far = searches(
        true,
        |k| create_cylinder(k, Point3::new(600.0, 0.0, 0.0), 4.0, 14.0),
        3,
        PI * 16.0 * 14.0,
    );
    let pointed = searches(
        true,
        |k| create_cone(k, p, 5.0, 0.0, 7.0),
        2,
        PI * 7.0 * 25.0 / 3.0,
    );
    let frustum = searches(
        true,
        |k| create_cone(k, p, 5.0, 2.0, 7.0),
        3,
        PI * 7.0 * (25.0 + 10.0 + 4.0) / 3.0,
    );
    assert_eq!((thin, far, pointed, frustum), (2, 2, 2, 2));
}
