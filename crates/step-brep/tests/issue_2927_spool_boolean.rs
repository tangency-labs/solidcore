//! Analytical and invariant verification of issue #2927:
//! `boolean_difference` and `boolean_intersection` on the imported 92-face
//! `flanged_spool.step` across half-space box cuts along X, Y and Z.
//!
//! Expected volumes are computed directly from the analytic surface and curve
//! definitions of `flanged_spool.step` (closed-form cross-sectional areas integrated in
//! `z` with 16-node Gauss-Legendre quadrature split at every profile break and chord
//! tangency kink, smoothed at endpoints by Hermite cubic substitution so `sqrt` chord
//! singularities become analytic).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use core::f64::consts::{FRAC_PI_2, FRAC_PI_3, FRAC_PI_6, PI};
use solidcore::boolean::api_brep::{boolean_difference, boolean_intersection};
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::mass::integrator::properties::{
    brep_properties, brep_volume, MassPropertyBehavior, PropertyAttributes,
};
use solidcore::primitive::api_primitives::{create_box, create_cylinder};
use solidcore::topology::brep::p1::faces;
use solidcore::topology::brep::validation::{check_brep, CheckLevel};
use solidcore::{Brep, Kernel, Point3};

/// 16-node Gauss-Legendre nodes and weights on `[-1, 1]`.
const GL16: [(f64, f64); 16] = [
    (-0.989_400_934_991_649_9, 0.027_152_459_411_754_095),
    (-0.944_575_023_073_232_6, 0.062_253_523_938_647_894),
    (-0.865_631_202_387_831_8, 0.095_158_511_682_492_78),
    (-0.755_404_408_355_003, 0.124_628_971_255_533_87),
    (-0.617_876_244_402_643_8, 0.149_595_988_816_576_73),
    (-0.458_016_777_657_227_37, 0.169_156_519_395_002_54),
    (-0.281_603_550_779_258_9, 0.182_603_415_044_923_58),
    (-0.095_012_509_837_637_44, 0.189_450_610_455_068_5),
    (0.095_012_509_837_637_44, 0.189_450_610_455_068_5),
    (0.281_603_550_779_258_9, 0.182_603_415_044_923_58),
    (0.458_016_777_657_227_37, 0.169_156_519_395_002_54),
    (0.617_876_244_402_643_8, 0.149_595_988_816_576_73),
    (0.755_404_408_355_003, 0.124_628_971_255_533_87),
    (0.865_631_202_387_831_8, 0.095_158_511_682_492_78),
    (0.944_575_023_073_232_6, 0.062_253_523_938_647_894),
    (0.989_400_934_991_649_9, 0.027_152_459_411_754_095),
];

/// Integrates `f(z)` on `[a, b]` using 16-node Gauss-Legendre after the smooth Hermite
/// substitution `z(u) = a + (b - a) * (3 u^2 - 2 u^3)` (`u in [0, 1]`), which regularizes
/// endpoint `sqrt(z - a)` and `sqrt(b - z)` chord-entry singularities to analytic terms.
fn integrate_smooth_span(a: f64, b: f64, f: impl Fn(f64) -> f64) -> f64 {
    let span = b - a;
    if span <= 1e-14 {
        return 0.0;
    }
    let mut sum = 0.0;
    for &(xi, wi) in &GL16 {
        let u = 0.5 * (xi + 1.0);
        let phi = u * u * (3.0 - 2.0 * u);
        let dphi = 6.0 * u * (1.0 - u);
        let z = a + span * phi;
        sum += 0.5 * wi * f(z) * (span * dphi);
    }
    sum
}

/// Area of a circular disk of radius `r` centered at the origin inside `{coord <= t}`.
fn clipped_disk_area(r: f64, t: f64) -> f64 {
    if r <= 0.0 || t <= -r {
        0.0
    } else if t >= r {
        PI * r * r
    } else {
        let ratio = (-t / r).clamp(-1.0, 1.0);
        r * r * ratio.acos() + t * (r * r - t * t).max(0.0).sqrt()
    }
}

/// Line integral `int_{t_start}^{t_end} min(x(t), d) dy(t)` along the circular arc
/// `(cx + r cos t, cy + r sin t)`. Summing around a closed CCW Jordan curve yields
/// `Area(P intersect {x <= d})` by Green's theorem.
fn arc_clipped_green_x(cx: f64, _cy: f64, r: f64, t_start: f64, t_end: f64, d: f64) -> f64 {
    let sign = if t_end >= t_start { 1.0 } else { -1.0 };
    let t_min = t_start.min(t_end);
    let t_max = t_start.max(t_end);
    if t_max - t_min <= 1e-15 {
        return 0.0;
    }

    let mut breaks = vec![t_min, t_max];
    let u = (d - cx) / r;
    if u > -1.0 && u < 1.0 {
        let alpha = u.acos();
        for k in -2..=3 {
            let shift = 2.0 * PI * f64::from(k);
            for cand in [alpha + shift, -alpha + shift] {
                if cand > t_min + 1e-13 && cand < t_max - 1e-13 {
                    breaks.push(cand);
                }
            }
        }
        breaks.sort_by(f64::total_cmp);
    }

    let mut total = 0.0;
    for pair in breaks.windows(2) {
        let a = pair[0];
        let b = pair[1];
        if b - a <= 1e-15 {
            continue;
        }
        let mid = 0.5 * (a + b);
        if cx + r * mid.cos() <= d {
            total += cx * r * (b.sin() - a.sin())
                + 0.5 * r * r * ((b - a) + 0.5 * ((2.0 * b).sin() - (2.0 * a).sin()));
        } else {
            total += d * r * (b.sin() - a.sin());
        }
    }
    sign * total
}

/// Green's-theorem contribution of a circular arc to `Area(P intersect {coord[axis] <= d})`.
/// For `axis == 1` (`y <= d`), applies the orientation-preserving `-90 deg` rotation
/// `(X, Y) = (y, -x)` (`dX ^ dY = dx ^ dy`, `{y <= d} == {X <= d}`).
fn arc_clipped_green(
    axis: usize,
    cx: f64,
    cy: f64,
    r: f64,
    t_start: f64,
    t_end: f64,
    d: f64,
) -> f64 {
    if axis == 0 {
        arc_clipped_green_x(cx, cy, r, t_start, t_end, d)
    } else {
        arc_clipped_green_x(cy, -cx, r, t_start - FRAC_PI_2, t_end - FRAC_PI_2, d)
    }
}

/// Outer revolved profile radius `r_out(z)` on `z in [0, 128]` from `flanged_spool.step`:
/// - `[0, 1]`: bottom chamfer cone `64 + z` (STEP Face 72; bottom face `z = 0` is Face 71)
/// - `[1, 4]`: bottom flange cylinder `65` (STEP Face 87)
/// - `[4, 5]`: bottom flange upper chamfer cone `69 - z` (STEP Face 90; top face `z = 5` is Face 88)
/// - `[5, 15]`: lower fillet torus `45 - sqrt(100 - (15 - z)^2)` (STEP Face 91)
/// - `[15, 98]`: barrel cylinder `35` (STEP Face 89)
/// - `[98, 108]`: middle lower fillet torus `45 - sqrt(100 - (z - 98)^2)` (STEP Face 86)
/// - `[108, 109]`: middle flange lower chamfer cone `64 + (z - 108)` (STEP Face 85; Face 69 at `z = 108`)
/// - `[109, 112]`: middle flange cylinder `65` (STEP Face 70)
/// - `[112, 113]`: middle flange upper chamfer cone `65 - (z - 112)` (STEP Face 56; Face 54 at `z = 113`)
/// - `[113, 123]`: upper fillet torus `45 - sqrt(100 - (123 - z)^2)` (STEP Face 52)
/// - `[123, 128]`: neck fillet torus `40 - sqrt(25 - (z - 123)^2)` (STEP Face 19)
fn revolved_outer_radius(z: f64) -> f64 {
    if z <= 1.0 {
        64.0 + z
    } else if z <= 4.0 {
        65.0
    } else if z <= 5.0 {
        69.0 - z
    } else if z <= 15.0 {
        let dz = 15.0 - z;
        45.0 - (100.0 - dz * dz).max(0.0).sqrt()
    } else if z <= 98.0 {
        35.0
    } else if z <= 108.0 {
        let dz = z - 98.0;
        45.0 - (100.0 - dz * dz).max(0.0).sqrt()
    } else if z <= 109.0 {
        64.0 + (z - 108.0)
    } else if z <= 112.0 {
        65.0
    } else if z <= 113.0 {
        65.0 - (z - 112.0)
    } else if z <= 123.0 {
        let dz = 123.0 - z;
        45.0 - (100.0 - dz * dz).max(0.0).sqrt()
    } else {
        let dz = z - 123.0;
        40.0 - (25.0 - dz * dz).max(0.0).sqrt()
    }
}

/// Area of the 12-arc scalloped top flange outer boundary `P(z)` (`z in [128, 140]`)
/// clipped to `{coord[axis] <= d}`:
/// - `[128, 129]`: lower chamfer cones, `delta = 129 - z` (STEP Faces 4..=15; Face 2 at `z = 128`)
/// - `[129, 139]`: 6 lobe cylinders `R1 = 10` + 6 concave scallop cylinders `R2 = 40` (STEP Faces 40..=51)
/// - `[139, 140]`: upper chamfer cones, `delta = z - 139` (STEP Faces 22..=33; Face 3 at `z = 140`)
fn scalloped_outer_clipped_area(z: f64, axis: usize, d: f64) -> f64 {
    let delta = if z < 129.0 {
        129.0 - z
    } else if z > 139.0 {
        z - 139.0
    } else {
        0.0
    };
    let r1 = 10.0 - delta;
    let r2 = 40.0 + delta;
    let d_scallop = 0.5 * (45.0 * 3.0_f64.sqrt() + 5.0 * 319.0_f64.sqrt());
    let s0x = d_scallop * FRAC_PI_6.cos();
    let s0y = d_scallop * FRAC_PI_6.sin();
    let alpha = s0y.atan2(s0x - 45.0);

    let mut area = 0.0;
    for m in 0..6 {
        let theta = f64::from(m) * FRAC_PI_3;
        let psi = theta + FRAC_PI_6;
        let cx = 45.0 * theta.cos();
        let cy = 45.0 * theta.sin();
        let sx = d_scallop * psi.cos();
        let sy = d_scallop * psi.sin();

        // Convex lobe arc around C_m (CCW from theta - alpha to theta + alpha).
        area += arc_clipped_green(axis, cx, cy, r1, theta - alpha, theta + alpha, d);
        // Concave scallop arc around S_m (CW from theta + pi + alpha down to theta + 4*pi/3 - alpha).
        area += arc_clipped_green(
            axis,
            sx,
            sy,
            r2,
            theta + PI + alpha,
            theta + 4.0 * FRAC_PI_3 - alpha,
            d,
        );
    }
    area
}

/// Exact cross-sectional area of `flanged_spool.step` at height `z in [0, 140]` inside
/// `{coord[axis] <= d}` (`d = +infinity` gives the full cross-section `A_full(z)`).
fn cross_section_area(z: f64, axis: usize, d: f64) -> f64 {
    // Inner bore: r = 31 on [0, 130] (STEP Face 55), step plane at z = 130 (Face 53),
    // r = 35 on [130, 140] (Face 34).
    let r_in = if z < 130.0 { 31.0 } else { 35.0 };

    if z < 128.0 {
        let r_out = revolved_outer_radius(z);
        let mut area = clipped_disk_area(r_out, d) - clipped_disk_area(r_in, d);
        // 12 bolt holes of radius 4 at R = 55 on bottom [0, 5] (Faces 73..=84) and
        // middle [108, 113] (Faces 57..=68) flanges.
        if z <= 5.0 || (108.0..=113.0).contains(&z) {
            for k in 0..12 {
                let ang = f64::from(k) * FRAC_PI_6;
                let c = if axis == 0 {
                    55.0 * ang.cos()
                } else {
                    55.0 * ang.sin()
                };
                area -= clipped_disk_area(4.0, d - c);
            }
        }
        area
    } else {
        let mut area = scalloped_outer_clipped_area(z, axis, d) - clipped_disk_area(r_in, d);
        // 6 countersunk bolt holes at R = 45: cylinder r = 5 on [128, 139]
        // (STEP Faces 0, 16..=18, 20..=21) and cone r = 5 + (z - 139) on [139, 140]
        // (STEP Faces 1, 35..=39).
        let r_hole = if z <= 139.0 { 5.0 } else { 5.0 + (z - 139.0) };
        for m in 0..6 {
            let ang = f64::from(m) * FRAC_PI_3;
            let c = if axis == 0 {
                45.0 * ang.cos()
            } else {
                45.0 * ang.sin()
            };
            area -= clipped_disk_area(r_hole, d - c);
        }
        area
    }
}

/// Exact analytical volume of `A \ {coord[axis] >= at}` (`A intersect {coord[axis] <= at}`)
/// computed from `flanged_spool.step`'s analytic definitions.
fn expected_cut_volume(axis: usize, at: f64) -> f64 {
    if axis == 2 {
        // Horizontal cut z <= at: integrate full cross-section A_full(z) from 0 to clamp(at, 0, 140).
        let z_hi = at.clamp(0.0, 140.0);
        let base_breaks = [
            0.0, 1.0, 4.0, 5.0, 15.0, 98.0, 108.0, 109.0, 112.0, 113.0, 123.0, 128.0, 129.0, 130.0,
            139.0, 140.0,
        ];
        let mut z_breaks = vec![0.0, z_hi];
        for &b in &base_breaks {
            if b > 0.0 && b < z_hi {
                z_breaks.push(b);
            }
        }
        z_breaks.sort_by(f64::total_cmp);
        let mut vol = 0.0;
        for pair in z_breaks.windows(2) {
            let a = pair[0];
            let b = pair[1];
            if b - a > 1e-14 {
                vol += integrate_smooth_span(a, b, |z| cross_section_area(z, 0, 1e6));
            }
        }
        return vol;
    }

    // Vertical cut x <= at (axis == 0) or y <= at (axis == 1): split at every profile
    // break and every chord-tangency kink where r(z) == |at - c|.
    let mut z_breaks = vec![
        0.0, 1.0, 4.0, 5.0, 15.0, 98.0, 108.0, 109.0, 112.0, 113.0, 123.0, 128.0, 129.0, 130.0,
        139.0, 140.0,
    ];
    let ad = at.abs();
    if ad > 64.0 && ad < 65.0 {
        z_breaks.push(ad - 64.0);
        z_breaks.push(69.0 - ad);
        z_breaks.push(108.0 + (ad - 64.0));
        z_breaks.push(112.0 + (65.0 - ad));
    }
    if ad > 35.0 && ad < 45.0 {
        let dz = (100.0 - (45.0 - ad) * (45.0 - ad)).max(0.0).sqrt();
        z_breaks.push(15.0 - dz);
        z_breaks.push(98.0 + dz);
        z_breaks.push(123.0 - dz);
    }
    if ad > 35.0 && ad < 40.0 {
        let dz = (25.0 - (40.0 - ad) * (40.0 - ad)).max(0.0).sqrt();
        z_breaks.push(123.0 + dz);
    }
    // Top scalloped flange chamfer kinks on [128, 129] and [139, 140]:
    let d_scallop = 0.5 * (45.0 * 3.0_f64.sqrt() + 5.0 * 319.0_f64.sqrt());
    for m in 0..6 {
        let theta = f64::from(m) * FRAC_PI_3;
        let psi = theta + FRAC_PI_6;
        let c_lobe = if axis == 0 {
            45.0 * theta.cos()
        } else {
            45.0 * theta.sin()
        };
        let c_scal = if axis == 0 {
            d_scallop * psi.cos()
        } else {
            d_scallop * psi.sin()
        };
        // Countersink cone r_h(z) = 5 + (z - 139) == |at - c_lobe| on [139, 140]
        let dh = (at - c_lobe).abs();
        if dh > 5.0 && dh < 6.0 {
            z_breaks.push(139.0 + (dh - 5.0));
        }
        // Lobe chamfer cone R1(z) = 10 - delta(z) == |at - c_lobe| (R1 in (9, 10))
        if dh > 9.0 && dh < 10.0 {
            let delta = 10.0 - dh;
            z_breaks.push(129.0 - delta);
            z_breaks.push(139.0 + delta);
        }
        // Scallop chamfer cone R2(z) = 40 + delta(z) == |at - c_scal| (R2 in (40, 41))
        let ds = (at - c_scal).abs();
        if ds > 40.0 && ds < 41.0 {
            let delta = ds - 40.0;
            z_breaks.push(129.0 - delta);
            z_breaks.push(139.0 + delta);
        }
    }
    // Subdivide long cylinder barrel [15, 98] and [129, 139] for uniform high precision.
    for z_mid in [35.0, 56.0, 77.0, 134.0] {
        z_breaks.push(z_mid);
    }

    z_breaks.sort_by(f64::total_cmp);
    let mut vol = 0.0;
    for pair in z_breaks.windows(2) {
        let a = pair[0];
        let b = pair[1];
        if b - a > 1e-12 {
            let mid = 0.5 * (a + b);
            vol += integrate_smooth_span(a, mid, |z| cross_section_area(z, axis, at));
            vol += integrate_smooth_span(mid, b, |z| cross_section_area(z, axis, at));
        }
    }
    vol
}

/// Pure closed-form antiderivative of the full unclipped spool volume `V(A)` as an
/// independent check on the quadrature routine:
/// - Cone frustum `[0, h]` with linear radius `r0 -> r1`: `pi * h / 3 * (r0^2 + r0 r1 + r1^2)`
/// - Concave torus span of height `H` (`H = rho` quarter-torus, `r(u) = R0 - sqrt(rho^2 - u^2)`):
///   `pi * int_0^rho (R0^2 - 2 R0 sqrt(rho^2 - u^2) + rho^2 - u^2) du`
///   `= pi * (R0^2 rho - pi/2 R0 rho^2 + 2/3 rho^3)`.
fn closed_form_full_volume() -> f64 {
    let frustum = |r0: f64, r1: f64, h: f64| PI * h / 3.0 * (r0 * r0 + r0 * r1 + r1 * r1);
    let quarter_concave_torus = |r0: f64, rho: f64| {
        PI * (r0 * r0 * rho - 0.5 * PI * r0 * rho * rho + 2.0 / 3.0 * rho.powi(3))
    };

    // Outer solid on [0, 128]:
    let mut v = 0.0;
    v += frustum(64.0, 65.0, 1.0); // [0, 1]
    v += PI * 65.0 * 65.0 * 3.0; // [1, 4]
    v += frustum(65.0, 64.0, 1.0); // [4, 5]
    v += quarter_concave_torus(45.0, 10.0); // [5, 15]
    v += PI * 35.0 * 35.0 * 83.0; // [15, 98]
    v += quarter_concave_torus(45.0, 10.0); // [98, 108]
    v += frustum(64.0, 65.0, 1.0); // [108, 109]
    v += PI * 65.0 * 65.0 * 3.0; // [109, 112]
    v += frustum(65.0, 64.0, 1.0); // [112, 113]
    v += quarter_concave_torus(45.0, 10.0); // [113, 123]
    v += quarter_concave_torus(40.0, 5.0); // [123, 128]
                                           // Subtract bore r = 31 on [0, 128] and 2 * 12 bolt holes (r = 4, h = 5) on [0, 5] & [108, 113]:
    v -= PI * 31.0 * 31.0 * 128.0;
    v -= 24.0 * PI * 4.0 * 4.0 * 5.0;

    // Top scalloped flange on [128, 140]:
    // Notice Area(P(delta)) is a quadratic polynomial in delta (`[0, 1]`), so Simpson's
    // rule `(A(0) + 4 A(0.5) + A(1)) / 6` is algebraically exact for `int_0^1 Area(P(delta)) d(delta)`!
    let a0 = scalloped_outer_clipped_area(134.0, 0, 1e6); // delta = 0
    let a_half = scalloped_outer_clipped_area(128.5, 0, 1e6); // delta = 0.5
    let a1 = scalloped_outer_clipped_area(128.0, 0, 1e6); // delta = 1.0
    let chamfer_vol = (a0 + 4.0 * a_half + a1) / 6.0;
    v += 10.0 * a0 + 2.0 * chamfer_vol;
    // Subtract upper bore (`r = 31` on `[128, 130]`, `r = 35` on `[130, 140]`) and 6 countersunk holes:
    v -= PI * 31.0 * 31.0 * 2.0 + PI * 35.0 * 35.0 * 10.0;
    v -= 6.0 * (PI * 5.0 * 5.0 * 11.0 + frustum(5.0, 6.0, 1.0));
    v
}

fn volume(k: &Kernel, b: &Brep) -> f64 {
    if faces(b, None).unwrap().is_empty() {
        return 0.0;
    }
    let store = AttributeStore::default();
    let attrs = PropertyAttributes {
        store: &store,
        key_of: &|_| ObjectKey(0),
    };
    let mut copy = b.clone_standing_for_self().unwrap();
    brep_volume(k, &mut copy, 1e-8, &attrs).unwrap()
}

fn volume_about_zero(k: &Kernel, b: &Brep) -> f64 {
    if faces(b, None).unwrap().is_empty() {
        return 0.0;
    }
    let store = AttributeStore::default();
    let attrs = PropertyAttributes {
        store: &store,
        key_of: &|_| ObjectKey(0),
    };
    let mut copy = b.clone_standing_for_self().unwrap();
    brep_properties(
        k,
        &mut copy,
        1e-8,
        Point3::ZERO,
        MassPropertyBehavior::default(),
        &attrs,
    )
    .unwrap()
    .volume
}

fn assert_valid_brep(k: &mut Kernel, b: &Brep, label: &str) {
    let mut copy = b.clone_standing_for_self().unwrap();
    let chk = check_brep(k, &mut copy, CheckLevel::Level2).unwrap();
    assert!(
        chk.pointer_reports.is_empty(),
        "{label}: pointer errors: {:?}",
        chk.pointer_reports
    );
    assert_eq!(
        chk.assert_valid,
        Some(true),
        "{label}: assert_valid failed: {:?}",
        chk.reports.reports()
    );
    assert!(
        chk.reports.reports().is_empty(),
        "{label}: validation reports: {:?}",
        chk.reports.reports()
    );
}

fn halfspace_box(k: &mut Kernel, axis: usize, at: f64) -> Brep {
    let mut o = [-1000.0; 3];
    let mut s = [2000.0; 3];
    o[axis] = at;
    s[axis] = 1000.0 - at;
    create_box(k, Point3::new(o[0], o[1], o[2]), s[0], s[1], s[2]).unwrap()
}

/// Maximum B-rep integration tolerance (mm^3) against exact analytic volume:
/// `Plane::drop_curve` fits cubic B-spline UV trims to circular/conic/toric cap edges
/// within `approx_tol = 5e-6 mm`; over a ~600 mm cap perimeter and 140 mm part height
/// (`brep_volume` integrates divergence from each body's bounding-box midpoint), that
/// `5e-6 mm` trim-curve fit accounts for at most `~0.35 mm^3` (`1.4e-6` relative to
/// `V(A) = 248052.03 mm^3`). We bound two-sided error at `0.45 mm^3` (`< 2e-6 * V(A)`).
const VOLUME_TOL_MM3: f64 = 0.45;

fn import_verified_spool(k: &mut Kernel) -> (Brep, f64, f64, f64) {
    let v_exact_full = closed_form_full_volume();
    let v_quad_full = expected_cut_volume(2, 140.0);
    assert!(
        (v_exact_full - v_quad_full).abs() < 1e-8,
        "quadrature {v_quad_full} vs closed form {v_exact_full}"
    );
    assert!(
        (expected_cut_volume(0, 0.0) - 0.5 * v_exact_full).abs() < 1e-8,
        "x <= 0 half-volume mismatch"
    );
    assert!(
        (expected_cut_volume(1, 0.0) - 0.5 * v_exact_full).abs() < 1e-8,
        "y <= 0 half-volume mismatch"
    );

    let bytes = include_bytes!("fixtures/flanged_spool.step");
    let report = step_brep::import_with_report(bytes, &step_import::Limits::default()).unwrap();
    assert_eq!(report.bodies.len(), 1);
    assert!(
        report.face_failures.is_empty(),
        "unexpected STEP face failures: {:?}",
        report.face_failures
    );
    let spool = report.bodies.into_iter().next().unwrap();
    assert_eq!(faces(&spool, None).unwrap().len(), 92);

    assert_valid_brep(k, &spool, "imported spool");
    let v_orig = volume(k, &spool);
    let v_orig_zero = volume_about_zero(k, &spool);
    assert!(
        (v_orig - v_exact_full).abs() <= VOLUME_TOL_MM3,
        "imported spool volume {v_orig} vs exact {v_exact_full}"
    );
    (spool, v_exact_full, v_orig, v_orig_zero)
}

fn verify_axis_sweep(
    k: &mut Kernel,
    spool: &Brep,
    v_exact_full: f64,
    v_orig: f64,
    v_orig_zero: f64,
    axis: usize,
    positions: &[i64],
    symmetric_ds: &[i64],
    max_area: f64,
) {
    let ax = ["x", "y", "z"][axis];
    let mut measured: std::collections::BTreeMap<i64, (f64, f64)> =
        std::collections::BTreeMap::new();

    for &d_int in positions {
        let at = d_int as f64;
        let label_diff = format!("diff {ax} > {at}");
        let label_inter = format!("inter {ax} > {at}");

        let tool_diff = halfspace_box(k, axis, at);
        let mut diff = boolean_difference(k, spool.clone_standing_for_self().unwrap(), tool_diff)
            .unwrap_or_else(|e| panic!("{label_diff} returned Err: {e:?}"));
        diff.refresh_cached_topology();

        let tool_inter = halfspace_box(k, axis, at);
        let mut inter =
            boolean_intersection(k, spool.clone_standing_for_self().unwrap(), tool_inter)
                .unwrap_or_else(|e| panic!("{label_inter} returned Err: {e:?}"));
        inter.refresh_cached_topology();

        let diff_empty = faces(&diff, None).unwrap().is_empty();
        let inter_empty = faces(&inter, None).unwrap().is_empty();
        if !diff_empty {
            assert_valid_brep(k, &diff, &label_diff);
        }
        if !inter_empty {
            assert_valid_brep(k, &inter, &label_inter);
        }

        let v_diff = volume(k, &diff);
        let v_inter = volume(k, &inter);
        let v_exp_diff = expected_cut_volume(axis, at);
        let v_exp_inter = v_exact_full - v_exp_diff;

        // Range invariant: 0 <= V(A \ H) <= V(A), empty iff the half-space clears the body.
        if v_exp_diff <= 1e-12 {
            assert!(
                diff_empty,
                "{label_diff}: expected empty B-rep, got {v_diff}"
            );
        } else if (v_exp_diff - v_exact_full).abs() <= 1e-10 {
            assert!(
                inter_empty,
                "{label_inter}: expected empty B-rep, got {v_inter}"
            );
        } else {
            assert!(
                v_diff > 0.0 && v_diff < v_orig,
                "{label_diff}: {v_diff} out of (0, {v_orig})"
            );
            assert!(
                v_inter > 0.0 && v_inter < v_orig,
                "{label_inter}: {v_inter} out of (0, {v_orig})"
            );
        }

        // Two-sided agreement with closed-form / quadrature reference.
        assert!(
            (v_diff - v_exp_diff).abs() <= VOLUME_TOL_MM3,
            "{label_diff}: got {v_diff:.6}, expected {v_exp_diff:.6} (err {:+.6})",
            v_diff - v_exp_diff
        );
        assert!(
            (v_inter - v_exp_inter).abs() <= VOLUME_TOL_MM3,
            "{label_inter}: got {v_inter:.6}, expected {v_exp_inter:.6} (err {:+.6})",
            v_inter - v_exp_inter
        );

        // Additivity invariant: V(A \ H) + V(A intersect H) == V(A) (both default-origin and
        // common-origin `Point3::ZERO`, where shared planar cap trim approximations cancel).
        assert!(
            (v_diff + v_inter - v_orig).abs() <= VOLUME_TOL_MM3,
            "additivity {ax} > {at}: diff {v_diff:.6} + inter {v_inter:.6} - orig {v_orig:.6} = {:+.6}",
            (v_diff + v_inter) - v_orig
        );
        let sum_zero = volume_about_zero(k, &diff) + volume_about_zero(k, &inter);
        assert!(
            (sum_zero - v_orig_zero).abs() <= 1e-2,
            "common-origin additivity {ax} > {at}: sum {sum_zero:.6} vs orig {v_orig_zero:.6}"
        );

        measured.insert(d_int, (v_diff, v_inter));
    }

    // Mirror symmetry invariants (`x -> -x` and `y -> -y`).
    for &d in symmetric_ds {
        let (v_diff_pos, v_inter_pos) = measured[&d];
        let (v_diff_neg, v_inter_neg) = measured[&-d];
        assert!(
            (v_diff_pos - v_inter_neg).abs() <= 0.05,
            "{ax} mirror symmetry at d={d}: diff(+d)={v_diff_pos:.6} vs inter(-d)={v_inter_neg:.6}"
        );
        assert!(
            (v_diff_neg - v_inter_pos).abs() <= 0.05,
            "{ax} mirror symmetry at d=-{d}: diff(-d)={v_diff_neg:.6} vs inter(+d)={v_inter_pos:.6}"
        );
    }

    // Monotonicity and Lipschitz continuity across the sorted cut positions.
    for pair in positions.windows(2) {
        let d1 = pair[0];
        let d2 = pair[1];
        let v1 = measured[&d1].0;
        let v2 = measured[&d2].0;
        let step = (d2 - d1) as f64;
        assert!(
            v2 >= v1 - 1e-6,
            "axis {ax}: non-monotone volume from {d1} ({v1}) to {d2} ({v2})"
        );
        assert!(
            v2 - v1 <= max_area * step + VOLUME_TOL_MM3,
            "axis {ax}: jump {} exceeds max_area * step between {d1} and {d2}",
            v2 - v1
        );
    }
}

#[test]
fn spool_boolean_y_cuts_closed_forms_and_invariants() {
    let mut k = Kernel::default();
    let (spool, v_exact_full, v_orig, v_orig_zero) = import_verified_spool(&mut k);
    verify_axis_sweep(
        &mut k,
        &spool,
        v_exact_full,
        v_orig,
        v_orig_zero,
        1,
        &[
            -70, -60, -45, -30, -20, -10, -5, 0, 5, 10, 20, 30, 45, 60, 70,
        ],
        &[0, 5, 10, 20, 30, 45, 60, 70],
        18_200.0,
    );
}

#[test]
fn spool_boolean_x_cuts_closed_forms_and_invariants() {
    let mut k = Kernel::default();
    let (spool, v_exact_full, v_orig, v_orig_zero) = import_verified_spool(&mut k);
    verify_axis_sweep(
        &mut k,
        &spool,
        v_exact_full,
        v_orig,
        v_orig_zero,
        0,
        &[-70, -45, -26, -20, -10, 0, 10, 20, 45, 70],
        &[0, 10, 20, 45, 70],
        18_200.0,
    );
}

#[test]
fn spool_boolean_z_cuts_closed_forms_and_invariants() {
    let mut k = Kernel::default();
    let (spool, v_exact_full, v_orig, v_orig_zero) = import_verified_spool(&mut k);
    verify_axis_sweep(
        &mut k,
        &spool,
        v_exact_full,
        v_orig,
        v_orig_zero,
        2,
        &[-5, 10, 40, 70, 100, 130, 145],
        &[],
        10_300.0,
    );

    // Fine continuity check across the z = 130 inner-bore step plane (`129.5, 130.0, 130.5`).
    let mut prev_fine = None;
    for z_cut in [129.5, 130.0, 130.5] {
        let tool = halfspace_box(&mut k, 2, z_cut);
        let mut diff =
            boolean_difference(&mut k, spool.clone_standing_for_self().unwrap(), tool).unwrap();
        diff.refresh_cached_topology();
        assert_valid_brep(&mut k, &diff, &format!("fine z > {z_cut}"));
        let v = volume(&k, &diff);
        let v_exp = expected_cut_volume(2, z_cut);
        assert!(
            (v - v_exp).abs() <= VOLUME_TOL_MM3,
            "fine z > {z_cut}: {v} vs {v_exp}"
        );
        if let Some(prev_v) = prev_fine {
            let dv: f64 = v - prev_v;
            assert!(
                dv > 0.0 && dv <= 10_300.0 * 0.5 + 2.0 * VOLUME_TOL_MM3,
                "discontinuous jump {dv} across z = 130 step"
            );
        }
        prev_fine = Some(v);
    }
}

#[test]
fn primitive_cylinder_control_cuts_from_issue_2927() {
    // Control from #2927: solid cylinder of matching outer envelope (r = 65, h = 140)
    // cut by the same half-space boxes matches closed-form circular-segment volumes.
    let mut k = Kernel::default();
    let cyl = create_cylinder(&mut k, Point3::ZERO, 65.0, 140.0).unwrap();
    for (axis, at, expected) in [
        (1usize, 0.0, 0.5 * PI * 65.0 * 65.0 * 140.0),
        (1usize, -30.0, clipped_disk_area(65.0, -30.0) * 140.0),
        (2usize, 70.0, PI * 65.0 * 65.0 * 70.0),
    ] {
        let tool = halfspace_box(&mut k, axis, at);
        let mut diff =
            boolean_difference(&mut k, cyl.clone_standing_for_self().unwrap(), tool).unwrap();
        diff.refresh_cached_topology();
        assert_valid_brep(&mut k, &diff, &format!("cylinder axis={axis} at={at}"));
        let v = volume(&k, &diff);
        assert!(
            (v - expected).abs() < 1e-2,
            "cylinder axis={axis} at={at}: got {v}, expected {expected}"
        );
    }
}
