// Audit-only bridge for the three continuous fixture nets; no production edits.
// With exact normal Z, the double-cross squared norm is exactly dx^2+dy^2.
// Point enclosures use native defining NURBS evaluation, not sampled geometry.
#[path = "projected_residual_audit_copy.rs"]
mod residual;
use solidcore::blend::fillet_section::{
    certificate::nurbs_partials_over_piece,
    interval::Ivl,
    station::{as_b_spline, Nurbs},
};
fn knot_cuts(
    knots: &[f64],
    p: usize,
    n: usize,
    lo: f64,
    hi: f64,
    _continuity: usize,
) -> solidcore::Result<Vec<(f64, f64, usize)>> {
    assert!(lo >= knots[p] && hi <= knots[n] && lo <= hi);
    if lo == hi {
        let span = if lo == knots[n] {
            n - 1
        } else {
            (p..n)
                .find(|&i| knots[i] <= lo && lo < knots[i + 1])
                .unwrap()
        };
        return Ok(vec![(lo, hi, span)]);
    }
    Ok((p..n)
        .filter_map(|i| {
            let a = lo.max(knots[i]);
            let b = hi.min(knots[i + 1]);
            (a < b).then_some((a, b, i))
        })
        .collect())
}
fn own_metric(net: &Nurbs, uv: Vec2, q: Vec3) -> solidcore::Result<Ivl> {
    let ku = knot_cuts(&net.u, net.p, net.nu, uv.x, uv.x, 0)?[0].2;
    let kv = knot_cuts(&net.v, net.q, net.nv, uv.y, uv.y, 0)?[0].2;
    let a = nurbs_partials_over_piece(net, ku, kv, uv.x, uv.x, uv.y, uv.y, 0)?.at(0, 0);
    let x = a.x - Ivl::from(q.x);
    let y = a.y - Ivl::from(q.y);
    // Ivl::checked rounds both computed square endpoints outward; see
    // blend/fillet_section/interval.rs:77 (round_down/up at lines 37/44).
    // Subsequent additions and sqrt also use checked outward endpoints.
    let sq = |z: Ivl| Ivl::checked(z.mig() * z.mig(), z.mag() * z.mag());
    Ok((sq(x) + sq(y) + sq(Ivl::from(0.))).sqrt())
}

use solidcore::curve::curve::Curve;
use solidcore::geom::{axis2_placement::Axis2Placement, extent2d::Extent2d, vector2d::Vec2};
use solidcore::nurbs::{
    b_spline_curve::p1::BSplineCurve,
    crv_basic::Curve as NlCurve,
    math::{CPoint, NOW},
};
use solidcore::surface::{
    cone::Cone,
    sphere::{Sphere, SphereOptions},
    surface::Surface,
};
use solidcore::{math::Vec3, Kernel};
fn main() -> solidcore::Result<()> {
    let k = Kernel::default();
    let curve = Curve::from(BSplineCurve::from_nurb(
        NlCurve::new(
            vec![CPoint::new(-1., 0., 0., NOW), CPoint::new(1., 0., 0., NOW)],
            1,
            vec![0., 0., 1., 1.],
        )?,
        3,
        false,
    ));
    let tangent = (30. * std::f64::consts::PI / 180.).tan();
    let cone = Surface::from(Cone::create_canonical(
        &k,
        &Axis2Placement::new(Vec3::new(-0.5, 0., 0.), Vec3::Y, Vec3::Z),
        0.5 * tangent,
        1.5 * tangent,
        1.,
        false,
        false,
    )?);
    let mut surfaces = vec![("cone_x", cone, Vec2::new(0., 1.), 3_f64.sqrt())];
    for (name, lon, lat, u) in [
        ("sphere_band", 360., 60., 2.3748193148539247e-9),
        ("sphere_quarter", 180., 45., 2.0053231314788843e-9),
    ] {
        let mut sphere = Sphere::new(
            &k,
            Vec3::new(0., 0., -1.),
            Vec3::X,
            Vec3::Y,
            Extent2d::from_bounds(0., -lat, lon, lat),
            0.5,
            SphereOptions::default(),
        )?;
        sphere.make_nurb(&k)?;
        surfaces.push((
            name,
            Surface::from(sphere),
            Vec2::new(u, 0.5000000000000001),
            1.5,
        ));
    }
    let mut torus = solidcore::surface::torus::Torus::new(
        &k,
        Vec3::new(0., 0., -2.),
        Vec3::Y,
        Vec3::Z,
        Extent2d::from_bounds(0., 0., 180., 360.),
        1.,
        0.25,
        solidcore::surface::torus::TorusOptions::default(),
    )?;
    torus.make_nurb(&k)?;
    surfaces.push((
        "torus_half_x",
        Surface::from(torus),
        Vec2::new(0., 0.12499999813735485),
        2_f64.sqrt() + 0.25,
    ));
    for (name, surface, uv, ideal) in surfaces {
        let domain = surface.natural_uv_domain_with_nurb(&k)?;
        assert!(domain.contains_point(uv, 0.));
        let c = curve.evaluate_point(&k, 0.)?;
        let p = surface.evaluate_point(&k, uv)?;
        let objective = Vec3::Z.cross(p - c).cross(Vec3::Z).length();
        let net = Nurbs::of(&as_b_spline(&k, &surface)?.unwrap())?;
        let own = own_metric(&net, uv, c)?;
        assert!(own.lo <= objective && objective <= own.hi);
        let left = residual::maximum_is_bounded(
            &k,
            &net,
            &domain,
            Vec3::new(-1., 0., 0.),
            Vec3::Z,
            own.hi,
        )?;
        let right = residual::maximum_is_bounded(
            &k,
            &net,
            &domain,
            Vec3::new(1., 0., 0.),
            Vec3::Z,
            own.hi,
        )?;
        eprintln!(
            "BOUND {name} own=[{:.17e},{:.17e}] left={left} right={right}",
            own.lo, own.hi
        );
        assert!(left && right);
        let negative =
            residual::maximum_is_bounded(&k, &net, &domain, c, Vec3::Z, objective.next_down())?;
        assert!(!negative);
        println!("{{\"proof\":\"{name}\",\"own\":[{:.17e},{:.17e}],\"left_bounded\":{left},\"right_bounded\":{right},\"reduced_bound_rejected\":{}}}",own.lo,own.hi,!negative);
        let boundary = surface.evaluate_point(
            &k,
            if name == "cone_x" {
                Vec2::new(0., 1.)
            } else {
                Vec2::new(0., 0.5)
            },
        )?;
        for entry in [6, 7] {
            println!("{{\"case\":\"analytic_face_{}_{entry}_15\",\"parameters\":[0,{:.17e},{:.17e}],\"point\":[{:.17e},{:.17e},{:.17e}],\"curve_point\":[{:.17e},{:.17e},{:.17e}],\"native_objective\":{:.17e},\"ideal_maximum\":{:.17e},\"boundary_point\":[{:.17e},{:.17e},{:.17e}],\"admissible\":true}}",name,uv.x,uv.y,p.x,p.y,p.z,c.x,c.y,c.z,objective,ideal,boundary.x,boundary.y,boundary.z);
        }
        eprintln!("{name} natural domain: {domain:?}");
    }
    Ok(())
}
