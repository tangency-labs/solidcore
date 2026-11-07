use solidcore::geom::extent2d::Extent2d;
use solidcore::surface::{
    surface::Surface,
    torus::{Torus, TorusOptions},
};
use solidcore::{
    math::{Vec2, Vec3},
    Kernel,
};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let k = Kernel::default();
    {
        let mut torus = Torus::new(
            &k,
            Vec3::new(0.0, 0.0, -2.0),
            Vec3::Y,
            Vec3::Z,
            Extent2d::from_bounds(0.0, 0.0, 360.0, 360.0),
            1.0,
            0.25,
            TorusOptions::default(),
        )?;
        torus.make_nurb(&k)?;
        let face = Surface::from(torus);
        let uv = Vec2::new(0.5, 0.12494065250506982);
        let domain = face.natural_uv_domain_with_nurb(&k)?;
        let p = face.evaluate_point(&k, uv)?;
        let line = Vec3::new(-1.0 + 2.0 * 0.0, 0.0, 0.0);
        let measured = Vec3::Z.cross(line - p).cross(Vec3::Z).length();
        let bound = 2.0_f64.sqrt() + 0.25;
        println!("case=analytic_face_torus_x_6_15 admissible={} point={p:?} reported={:.17} measured={measured:.17} report_error={:.17e} analytic_bound={bound:.17} above_bound={:.17e} actual_gap={:.17e}",domain.contains_point(uv,0.0),1.6642154866673036,1.6642154866673036-measured,1.6642154866673036-bound,bound-measured);
    }
    {
        let mut torus = Torus::new(
            &k,
            Vec3::new(0.0, 0.0, -2.0),
            Vec3::Y,
            Vec3::Z,
            Extent2d::from_bounds(0.0, 0.0, 360.0, 360.0),
            1.0,
            0.25,
            TorusOptions::default(),
        )?;
        torus.make_nurb(&k)?;
        let face = Surface::from(torus);
        let uv = Vec2::new(0.5, 0.12494065250506982);
        let domain = face.natural_uv_domain_with_nurb(&k)?;
        let p = face.evaluate_point(&k, uv)?;
        let line = Vec3::new(-1.0 + 2.0 * 0.0, 0.0, 0.0);
        let measured = Vec3::Z.cross(line - p).cross(Vec3::Z).length();
        let bound = 2.0_f64.sqrt() + 0.25;
        println!("case=analytic_face_torus_x_7_15 admissible={} point={p:?} reported={:.17} measured={measured:.17} report_error={:.17e} analytic_bound={bound:.17} above_bound={:.17e} actual_gap={:.17e}",domain.contains_point(uv,0.0),1.6642154866673036,1.6642154866673036-measured,1.6642154866673036-bound,bound-measured);
    }
    {
        let mut torus = Torus::new(
            &k,
            Vec3::new(0.0, 0.0, -2.0),
            Vec3::Y,
            Vec3::Z,
            Extent2d::from_bounds(0.0, 0.0, 180.0, 360.0),
            1.0,
            0.25,
            TorusOptions::default(),
        )?;
        torus.make_nurb(&k)?;
        let face = Surface::from(torus);
        let uv = Vec2::new(0.0, 0.12499999813735485);
        let domain = face.natural_uv_domain_with_nurb(&k)?;
        let p = face.evaluate_point(&k, uv)?;
        let line = Vec3::new(-1.0 + 2.0 * 0.0, 0.0, 0.0);
        let measured = Vec3::Z.cross(line - p).cross(Vec3::Z).length();
        let bound = 2.0_f64.sqrt() + 0.25;
        println!("case=analytic_face_torus_half_x_6_15 admissible={} point={p:?} reported={:.17} measured={measured:.17} report_error={:.17e} analytic_bound={bound:.17} above_bound={:.17e} actual_gap={:.17e}",domain.contains_point(uv,0.0),1.664213562373095,1.664213562373095-measured,1.664213562373095-bound,bound-measured);
    }
    {
        let mut torus = Torus::new(
            &k,
            Vec3::new(0.0, 0.0, -2.0),
            Vec3::Y,
            Vec3::Z,
            Extent2d::from_bounds(0.0, 0.0, 180.0, 360.0),
            1.0,
            0.25,
            TorusOptions::default(),
        )?;
        torus.make_nurb(&k)?;
        let face = Surface::from(torus);
        let uv = Vec2::new(0.0, 0.12499999813735485);
        let domain = face.natural_uv_domain_with_nurb(&k)?;
        let p = face.evaluate_point(&k, uv)?;
        let line = Vec3::new(-1.0 + 2.0 * 0.0, 0.0, 0.0);
        let measured = Vec3::Z.cross(line - p).cross(Vec3::Z).length();
        let bound = 2.0_f64.sqrt() + 0.25;
        println!("case=analytic_face_torus_half_x_7_15 admissible={} point={p:?} reported={:.17} measured={measured:.17} report_error={:.17e} analytic_bound={bound:.17} above_bound={:.17e} actual_gap={:.17e}",domain.contains_point(uv,0.0),1.664213562373095,1.664213562373095-measured,1.664213562373095-bound,bound-measured);
    }
    Ok(())
}
