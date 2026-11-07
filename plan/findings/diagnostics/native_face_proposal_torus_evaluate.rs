// Standalone native evaluator for the frozen torus records. No solver is called.
// Build against the audited kernel rlib; pass major-angle, u, v.
use solidcore::geom::{extent2d::Extent2d, vector2d::Vec2};
use solidcore::kernel::Kernel;
use solidcore::math::Vec3;
use solidcore::surface::{
    surface::Surface,
    torus::{Torus, TorusOptions},
};
fn main() {
    let args: Vec<f64> = std::env::args()
        .skip(1)
        .map(|s| s.parse().unwrap())
        .collect();
    assert_eq!(args.len(), 3);
    let k = Kernel::default();
    let mut torus = Torus::new(
        &k,
        Vec3::new(0.0, 0.0, -2.0),
        Vec3::Y,
        Vec3::Z,
        Extent2d::from_bounds(0.0, 0.0, args[0], 360.0),
        1.0,
        0.25,
        TorusOptions::default(),
    )
    .unwrap();
    torus.make_nurb(&k).unwrap();
    let surface = Surface::from(torus);
    let domain = surface.natural_uv_domain_with_nurb(&k).unwrap();
    let uv = Vec2::new(args[1], args[2]);
    assert!(domain.contains_point(uv, 0.0));
    let p = surface.evaluate_point(&k, uv).unwrap();
    let value = Vec3::Z
        .cross(p - Vec3::new(0.2, 0.3, 2.0))
        .cross(Vec3::Z)
        .length();
    println!("{{\"point\":[{:.17e},{:.17e},{:.17e}],\"value\":{:.17e},\"domain\":[{:.17e},{:.17e},{:.17e},{:.17e}]}}",
        p.x, p.y, p.z, value, domain.min().x, domain.min().y, domain.max().x, domain.max().y);
}
