//! `cargo run --release -p solidcore-mesh --example dense_nurbs [--features nurbs-grid-mesh]
//! -- [<poles u> <poles v> <degree> <size> <chord>...]`
//!
//! Meshes one dense B-spline face (a wavy height field over a `size` x `size` square,
//! degree `degree` both ways, uniform clamped knots) made into a one-face sheet with
//! `CreateFaceFromSurface`, at each chord: time, triangles and measured deviation, with
//! the built-in grid and (feature `nurbs-grid-mesh`) with the nurbs-grid-mesh backend's
//! grid lines. Defaults: 318 x 281 poles, degree 7, size 490, chords 1.0, 0.312, 0.1.

use std::time::Instant;

use solidcore::geom::attribute::ObjectKey;
use solidcore::geom::extent2d::Extent2d;
use solidcore::geom::object::ObjectRef;
use solidcore::kernel::Kernel;
use solidcore::nurbs::b_spline_surface::BSplineSurface;
use solidcore::nurbs::math::{CPoint, NOW};
use solidcore::surface::surface::Surface;
use solidcore::topology::brep::p2::ArenaNotifyState;
use solidcore::topology::brep::p8::{create_face_from_surface, FaceBuildEnv};
use solidcore::topology::topology_traverser::Marks;
use solidcore::topology::Brep;
use solidcore::FaceId;
use solidcore_mesh::{tessellate, TessellateOptions};

/// Uniform clamped knots for `n` poles of degree `p` over `[0, 1]`.
fn knots(n: usize, p: usize) -> Vec<f64> {
    let spans = n - p;
    let mut k = vec![0.0; p + 1];
    for i in 1..spans {
        k.push(i as f64 / spans as f64);
    }
    k.extend(std::iter::repeat_n(1.0, p + 1));
    k
}

fn sheet(k: &mut Kernel, nu: usize, nv: usize, p: usize, size: f64) -> Result<Brep, String> {
    let pw: Vec<Vec<CPoint>> = (0..nu)
        .map(|i| {
            (0..nv)
                .map(|j| {
                    let (x, y) = (i as f64 / (nu - 1) as f64, j as f64 / (nv - 1) as f64);
                    // ridges a few knot spans wide, like a sculpted surface
                    let z = 0.04
                        * size
                        * ((x * 37.0).sin() * (y * 29.0).cos() + 0.5 * (x * y * 61.0).sin());
                    CPoint {
                        x: x * size,
                        y: y * size,
                        z,
                        w: NOW,
                    }
                })
                .collect()
        })
        .collect();
    let (u, v) = (knots(nu, p), knots(nv, p));
    let s = BSplineSurface::from_data(nu - 1, nv - 1, p, p, u.len() - 1, v.len() - 1, &pw, &u, &v)
        .map_err(|e| e.to_string())?;
    let mut b = Brep::new_with_infinite_region().map_err(|e| e.to_string())?;
    let mut state = ArenaNotifyState::new(&b);
    let mut env = state.env(None);
    let mut marks = Marks::default();
    let this_face = |f: FaceId| ObjectRef {
        type_code: 16002,
        key: ObjectKey((5 << 32) | u64::from(f.index())),
    };
    let mut fenv = FaceBuildEnv {
        env: &mut env,
        marks: &mut marks,
        this_face: &this_face,
    };
    create_face_from_surface(
        k,
        &mut b,
        Some(Surface::from(s)),
        &Extent2d::from_bounds(0.0, 0.0, 1.0, 1.0),
        &mut fenv,
    )
    .map_err(|e| e.to_string())?;
    Ok(b)
}

fn main() -> Result<(), String> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let num = |i: usize, d: f64| args.get(i).and_then(|a| a.parse().ok()).unwrap_or(d);
    let (nu, nv, p, size) = (
        num(0, 318.0) as usize,
        num(1, 281.0) as usize,
        num(2, 7.0) as usize,
        num(3, 490.0),
    );
    let chords: Vec<f64> = if args.len() > 4 {
        args.iter().skip(4).filter_map(|a| a.parse().ok()).collect()
    } else {
        vec![1.0, 0.312, 0.1]
    };
    let mut k = Kernel::default();
    let t = Instant::now();
    let b = sheet(&mut k, nu, nv, p, size)?;
    println!(
        "{nu} x {nv} poles, degree {p}, {size} wide: sheet made in {:.0} ms",
        t.elapsed().as_secs_f64() * 1e3
    );
    for chord in chords {
        #[allow(unused_mut)] // a second run with the feature
        let mut runs: Vec<(&str, TessellateOptions)> = vec![(
            "built-in grid",
            TessellateOptions::default().with_chord(chord),
        )];
        #[cfg(feature = "nurbs-grid-mesh")]
        {
            let mut o = TessellateOptions::default().with_chord(chord);
            o.backend = Some(std::sync::Arc::new(
                solidcore_mesh::backend::NurbsGridBackend::default(),
            ));
            runs.push(("nurbs-grid-mesh grid", o));
        }
        for (name, mut o) in runs {
            o.angle_tolerance = std::f64::consts::PI; // chord only, as nurbs-grid-mesh
            o.max_triangles_per_face = 5_000_000;
            if let Some(n) = std::env::var("PASSES").ok().and_then(|x| x.parse().ok()) {
                o.max_refine_passes = n;
            }
            let t = Instant::now();
            let m = tessellate(&k, &b, &o).map_err(|e| e.to_string())?;
            let ms = t.elapsed().as_secs_f64() * 1e3;
            if std::env::var("WORST").is_ok() {
                let d = solidcore_mesh::deviation(&k, &b, &m, chord).map_err(|e| e.to_string())?;
                if let Some(w) = d.worst {
                    println!(
                        "   worst {w}: {:?} uv {:?} ({} over)",
                        m.triangles.get(w).and_then(|t| m.triangle_points(*t)),
                        m.uvs.get(w),
                        d.over
                    );
                }
            }
            o.measure_deviation = true;
            let mm = tessellate(&k, &b, &o).map_err(|e| e.to_string())?;
            println!(
                "chord {chord}: {name}: {ms:.0} ms, {} triangles, deviation {:.3e} ({:.2}x chord), nurbs-grid-mesh deviation {}, unconverged {}, failures {}",
                m.triangles.len(),
                mm.stats.face_deviation.unwrap_or(f64::NAN),
                mm.stats.face_deviation.unwrap_or(f64::NAN) / chord,
                mm.stats.backend_deviation.map_or("-".into(), |d| format!("{d:.3e}")),
                m.unconverged_faces.len(),
                m.failures.len(),
            );
        }
    }
    Ok(())
}
