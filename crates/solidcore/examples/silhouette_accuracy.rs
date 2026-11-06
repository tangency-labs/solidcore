//! Diagnostic for #2967: silhouette samples versus the analytic sphere equator.
//! Sampling diagnoses error; it does not certify a maximum or a global bound.
use solidcore::curve::crv_on_surf::CrvOnSurf;
use solidcore::curve::curve::{p1::GlobalPointSolveOptions, Curve};
use solidcore::geom::extent2d::Extent2d;
use solidcore::geom::global_solver::{SolutionRequested, SolverOperation};
use solidcore::geom::vector2d::Vec2;
use solidcore::geom::vector3d::Vec3;
use solidcore::kernel::Kernel;
use solidcore::surface::sphere::{Sphere, SphereOptions};
use solidcore::surface::surface::Surface;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let k = Kernel::default();
    let mut sphere = Sphere::new(
        &k,
        Vec3::new(0.0, 0.0, -1.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Extent2d::from_bounds(0.0, -60.0, 360.0, 60.0),
        0.5,
        SphereOptions::default(),
    )?;
    sphere.make_nurb(&k)?;
    let surface = Surface::from(sphere);
    let domain = surface.natural_uv_domain_with_nurb(&k)?;
    let (curves, uv_curves) = surface.create_silhouette_curves(
        &k,
        &domain,
        Vec3::new(0.0, 0.0, 1.0),
        false,
        Some(1e-4),
        Some(20.0_f64.to_radians()),
    )?;
    let query = Vec3::new(0.2, 0.3, 2.0);
    let exact = 0.5 + (query.x * query.x + query.y * query.y).sqrt();
    eprintln!("curves={} exact_max={exact:.17e}", curves.len());
    println!("curve,t,x,y,z,sx,sy,sz,distance,surface_distance");
    let mut sampled_max = f64::NEG_INFINITY;
    let mut sampled_surface_max = f64::NEG_INFINITY;
    let mut sampled_gap = 0.0_f64;
    for (index, (curve, uv_curve)) in curves.iter().zip(&uv_curves).enumerate() {
        let lifted = Curve::from(CrvOnSurf::new(
            &k,
            uv_curve.clone(),
            surface.clone(),
            Some(domain),
            0,
        )?);
        for (label, c) in [("approximation", curve), ("on_surface", &lifted)] {
            let normals = [Vec3::new(0.0, 0.0, 1.0)];
            let sols = c.global_point_solve::<solidcore::surface::surface::p4::NoObject>(
                &k,
                &c.natural_interval(),
                SolverOperation::ProjectedMaximize,
                query,
                1e-5,
                SolutionRequested::Single,
                &GlobalPointSolveOptions {
                    opt_vectors: Some(&normals),
                    ..Default::default()
                },
            )?;
            for sol in sols.as_slice() {
                eprintln!(
                    "query={label} value={:.17e} t={:.17e}",
                    sol.start.solution_value,
                    sol.start.param(0)?
                );
            }
            eprintln!("query={label} count={}", sols.len());
        }
        let interval = curve.natural_interval();
        for step in 0..=4000 {
            let t = interval.evaluate(f64::from(step) / 4000.0);
            let p = curve.evaluate_point(&k, t)?;
            let uv = uv_curve.evaluate_point(&k, t)?;
            let s = surface.evaluate_point(&k, Vec2::new(uv.x, uv.y))?;
            let dist = ((p.x - query.x).powi(2) + (p.y - query.y).powi(2)).sqrt();
            let surface_dist = ((s.x - query.x).powi(2) + (s.y - query.y).powi(2)).sqrt();
            sampled_max = sampled_max.max(dist);
            sampled_surface_max = sampled_surface_max.max(surface_dist);
            sampled_gap = sampled_gap.max(p.distance(s));
            println!("{index},{t:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{dist:.17e},{surface_dist:.17e}", p.x,p.y,p.z,s.x,s.y,s.z);
        }
    }
    eprintln!("sampled_curve_max={sampled_max:.17e} sampled_surface_max={sampled_surface_max:.17e} sampled_curve_surface_gap={sampled_gap:.17e}");
    Ok(())
}
