//! Interval enclosures in the projected queries' double-cross metric.
//!
//! A projected query measures a difference `V` as `|N x V x N|`, which is
//! `(N . N) |P(V)|` for the orthogonal projection `P` onto the plane normal to
//! `N`. `N` need not be unit length, so tree bounds must use the same scaled
//! metric as the answers they are compared with.

use crate::blend::fillet_section::interval::Ivl;
use crate::blend::fillet_section::series::V3;
use crate::geom::extent3d::Extent3d;
use crate::math::Vec3;

fn cross(a: V3<Ivl>, b: V3<Ivl>) -> V3<Ivl> {
    V3 {
        x: a.y * b.z - a.z * b.y,
        y: a.z * b.x - a.x * b.z,
        z: a.x * b.y - a.y * b.x,
    }
}

/// `N x V x N`, the projected difference in the queries' metric.
pub(super) fn projected(v: V3<Ivl>, normal: V3<Ivl>) -> V3<Ivl> {
    cross(cross(normal, v), normal)
}

fn square(x: Ivl) -> Ivl {
    Ivl::checked(x.mig() * x.mig(), x.mag() * x.mag())
}

pub(super) fn norm_squared(v: V3<Ivl>) -> Ivl {
    (square(v.x) + square(v.y)) + square(v.z)
}

/// A conservative range of the projected distance between any point of `a` and
/// any point of `b`, in the double-cross metric. A non-finite or inverted box
/// gives the whole line, which can never exclude a branch.
pub(super) fn projected_box_distance(a: &Extent3d, b: &Extent3d, normal: Vec3) -> Ivl {
    if [a, b].iter().any(|extent| {
        let (lo, hi) = (extent.min(), extent.max());
        ![lo.x, lo.y, lo.z, hi.x, hi.y, hi.z]
            .iter()
            .all(|x| x.is_finite())
            || lo.x > hi.x
            || lo.y > hi.y
            || lo.z > hi.z
    }) {
        return Ivl::whole();
    }
    let delta = V3 {
        x: Ivl::new(b.min().x, b.max().x) - Ivl::new(a.min().x, a.max().x),
        y: Ivl::new(b.min().y, b.max().y) - Ivl::new(a.min().y, a.max().y),
        z: Ivl::new(b.min().z, b.max().z) - Ivl::new(a.min().z, a.max().z),
    };
    norm_squared(projected(delta, V3::from_vec3(normal))).sqrt()
}
