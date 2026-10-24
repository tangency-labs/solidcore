//! `cargo run --release -p solidcore-mesh --example tessellate -- [options] <model>...`
//!
//! Tessellates B-reps and writes or measures their meshes. A `<model>` is a B-rep file
//! (`.smb`, read with `SmBrep::ReadFromFile`), `prim:<name>` for one of the
//! `primitives` example's solids (`box`, `cylinder`, `cone`, `pyramid`, `sphere`,
//! `sphere-no-pole`, `partial-sphere`, `torus`, `partial-torus`) or `prims` for all
//! nine.
//!
//! Options:
//! - `--quality draft|normal|fine` (default normal), `--chord <abs>`, `--rel <fraction
//!   of the bounding-box diagonal>`, `--angle <degrees>`: the tolerances.
//! - `--stl <path>`: binary STL of the (last) model; `--json <path>`: the viewer JSON
//!   (schema below).
//! - `--edges-only`: only the edge polylines.
//! - `--grid-backend <min poles>` (feature `nurbs-grid-mesh`): grid lines of B-spline faces with at
//!   least that many control points come from the nurbs-grid-mesh backend.
//! - `--lods <r1,r2,...>`: nested levels at these relative chords ([`solidcore_mesh::Tessellator`]),
//!   each timed against an independent run, with its measured deviation.
//! - `--bench`: every model at relative chords 1e-2, 1e-3 and 1e-4, with time,
//!   counts, watertightness, measured deviation and volume and area against the
//!   B-rep's integrated ones (`--no-ref` skips the integration).
//!
//! JSON (schema 1): `{ "schema": 1, "source": "solidcore-mesh", "model", "chord",
//! "bbox": [x0, y0, z0, x1, y1, z1], "positions": [x, y, z, ...], "triangles": [i, j,
//! k, ...], "normals": [9 per triangle: per-corner unit normals], "triangle_faces": [face
//! index per triangle], "edges": [[[x, y, z], ...] per edge polyline] }`.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::time::Instant;

use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::geom::object::ObjectRef;
use solidcore::kernel::Kernel;
use solidcore::mass::integrator::properties::{brep_area, brep_volume, PropertyAttributes};
use solidcore::math::Point3;
use solidcore::primitive::api_primitives as api;
use solidcore::topology::brep::p10::TopologyEnv;
use solidcore::topology::brep::p7::read_from_file;
use solidcore::topology::database_io::FileType;
use solidcore::topology::Brep;
use solidcore::EntityRef;
use solidcore_mesh::{tessellate, tessellate_edges, Mesh, Quality, TessellateOptions, Tolerance};

const PRIMITIVES: [&str; 9] = [
    "box",
    "cylinder",
    "cone",
    "pyramid",
    "sphere",
    "sphere-no-pole",
    "partial-sphere",
    "torus",
    "partial-torus",
];

fn primitive(k: &mut Kernel, name: &str) -> solidcore::Result<Brep> {
    let o = Point3::new(0.0, 0.0, 0.0);
    let c = Point3::new(1.0, -2.0, 0.5);
    let t = Point3::new(0.0, 1.0, -1.0);
    match name {
        "box" => api::create_box(k, o, 5.0, 10.0, 20.0),
        "cylinder" => api::create_cylinder(k, o, 5.0, 20.0),
        "cone" => api::create_cone(k, o, 5.0, 2.0, 20.0),
        "pyramid" => api::create_pyramid(k, o, 5.0, 20.0),
        "sphere" => api::create_sphere(k, c, 3.0),
        "sphere-no-pole" => api::create_sphere_no_pole(k, c, 3.0),
        "partial-sphere" => api::create_partial_sphere(k, c, 3.0, 0.0, 90.0),
        "torus" => api::create_torus(k, t, 5.0, 1.5),
        "partial-torus" => api::create_partial_torus(k, t, 5.0, 1.5, 0.0, 90.0),
        _ => Err(solidcore::SolidError::invalid_input(
            "tessellate example",
            "unknown primitive",
        )),
    }
}

/// A distinct identity for every topology object the reader builds.
fn key_of(e: EntityRef) -> ObjectRef {
    let (type_code, kind, i) = match e {
        EntityRef::Brep(x) => (16011, 0u64, x.index()),
        EntityRef::Region(x) => (16012, 1, x.index()),
        EntityRef::Shell(x) => (16013, 2, x.index()),
        EntityRef::Faceuse(x) => (16022, 3, x.index()),
        EntityRef::Face(x) => (16002, 4, x.index()),
        EntityRef::Loopuse(x) => (16004, 5, x.index()),
        EntityRef::Loop(x) => (16003, 6, x.index()),
        EntityRef::Edgeuse(x) => (16005, 7, x.index()),
        EntityRef::Edge(x) => (16016, 8, x.index()),
        EntityRef::Vertexuse(x) => (16006, 9, x.index()),
        EntityRef::Vertex(x) => (16021, 10, x.index()),
    };
    ObjectRef {
        type_code,
        key: ObjectKey(((kind + 1) << 32) | u64::from(i)),
    }
}

fn read_model(k: &mut Kernel, path: &str) -> Result<Brep, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("{path}: {e}"))?;
    let mut b = Brep::new_with_infinite_region().map_err(|e| e.to_string())?;
    let mut attributes = AttributeStore::default();
    let mut records = BTreeMap::new();
    let mut env = TopologyEnv {
        all_attributes: &[],
        attributes: &mut attributes,
        key_of: &key_of,
        records: &mut records,
        next_key: 1 << 62,
    };
    read_from_file(k, &mut b, bytes, FileType::Ascii, false, &mut env)
        .map_err(|e| format!("{path}: {e}"))?;
    Ok(b)
}

fn load(k: &mut Kernel, model: &str) -> Result<Brep, String> {
    match model.strip_prefix("prim:") {
        Some(name) => primitive(k, name).map_err(|e| format!("{model}: {e}")),
        None => read_model(k, model),
    }
}

#[allow(clippy::struct_excessive_bools)] // command-line flags
struct Args {
    models: Vec<String>,
    opts: TessellateOptions,
    stl: Option<String>,
    json: Option<String>,
    edges_only: bool,
    lods: Option<Vec<f64>>,
    bench: bool,
    no_ref: bool,
    per_face: bool,
}

#[cfg(feature = "nurbs-grid-mesh")]
#[allow(clippy::unnecessary_wraps)] // the same signature as without the feature
fn grid_backend(
    min_poles: usize,
) -> Result<Option<std::sync::Arc<dyn solidcore_mesh::backend::FaceBackend>>, String> {
    Ok(Some(std::sync::Arc::new(
        solidcore_mesh::backend::NurbsGridBackend {
            min_poles,
            ..Default::default()
        },
    )))
}

#[cfg(not(feature = "nurbs-grid-mesh"))]
fn grid_backend(
    _: usize,
) -> Result<Option<std::sync::Arc<dyn solidcore_mesh::backend::FaceBackend>>, String> {
    Err("--grid-backend needs the nurbs-grid-mesh feature".into())
}

fn parse() -> Result<Args, String> {
    let mut a = Args {
        models: Vec::new(),
        opts: TessellateOptions::default(),
        stl: None,
        json: None,
        edges_only: false,
        lods: None,
        bench: false,
        no_ref: false,
        per_face: false,
    };
    let mut it = std::env::args().skip(1);
    let num = |s: Option<String>, what: &str| -> Result<f64, String> {
        s.ok_or(format!("{what} needs a value"))?
            .parse()
            .map_err(|e| format!("{what}: {e}"))
    };
    while let Some(x) = it.next() {
        match x.as_str() {
            "--quality" => {
                let q = match it.next().as_deref() {
                    Some("draft") => Quality::Draft,
                    Some("normal") => Quality::Normal,
                    Some("fine") => Quality::Fine,
                    _ => return Err("--quality draft|normal|fine".into()),
                };
                a.opts = TessellateOptions::preset(q);
            }
            "--chord" => a.opts.chord_tolerance = Tolerance::Absolute(num(it.next(), "--chord")?),
            "--rel" => a.opts.chord_tolerance = Tolerance::Relative(num(it.next(), "--rel")?),
            "--angle" => a.opts.angle_tolerance = num(it.next(), "--angle")?.to_radians(),
            "--stl" => a.stl = Some(it.next().ok_or("--stl needs a path")?),
            "--json" => a.json = Some(it.next().ok_or("--json needs a path")?),
            "--edges-only" => a.edges_only = true,
            "--max-tris" => a.opts.max_triangles_per_face = num(it.next(), "--max-tris")? as usize,
            "--passes" => a.opts.max_refine_passes = num(it.next(), "--passes")? as u32,
            "--per-face" => a.per_face = true,
            "--bench" => a.bench = true,
            "--no-ref" => a.no_ref = true,
            "--grid-backend" => {
                // nurbs-grid-mesh for B-spline faces of at least this many poles
                let min = num(it.next(), "--grid-backend")? as usize;
                a.opts.backend = grid_backend(min)?;
            }
            "--measure" => a.opts.measure_deviation = true,
            "--lods" => {
                let s = it.next().ok_or("--lods needs r1,r2,...")?;
                a.lods = Some(
                    s.split(',')
                        .map(|x| x.parse().map_err(|e| format!("--lods: {e}")))
                        .collect::<Result<_, _>>()?,
                );
            }
            "prims" => a
                .models
                .extend(PRIMITIVES.iter().map(|p| format!("prim:{p}"))),
            m if !m.starts_with("--") => a.models.push(m.to_string()),
            _ => return Err(format!("unknown option {x}")),
        }
    }
    if a.models.is_empty() {
        return Err("usage: tessellate [options] <model.smb | prim:<name> | prims>...".into());
    }
    Ok(a)
}

/// Integrated B-rep volume (solids) and area, for comparison.
fn reference(k: &Kernel, b: &Brep) -> (Option<f64>, Option<f64>) {
    let store = AttributeStore::new();
    let key = |_| ObjectKey(0);
    let attrs = PropertyAttributes {
        store: &store,
        key_of: &key,
    };
    let solid = solidcore::topology::brep::p2::is_manifold_solid(b).unwrap_or(false);
    let v = if solid {
        brep_volume(k, b, 1e-7, &attrs).ok()
    } else {
        None
    };
    (v, brep_area(k, b, 1e-7).ok())
}

/// Triangles with (almost) no area, and triangles a corner normal points against.
fn bad_triangles(m: &Mesh) -> (usize, usize) {
    let (mut degenerate, mut flipped) = (0, 0);
    let size2 = m.stats.model_size * m.stats.model_size;
    for (i, t) in m.triangles.iter().enumerate() {
        let Some([a, b, c]) = m.triangle_points(*t) else {
            continue;
        };
        let g = (b - a).cross(c - a);
        if g.length() <= 1e-14 * size2 {
            degenerate += 1;
            if std::env::var("WORST").is_ok() {
                println!(
                    "   degenerate {t:?}: {a:?} {b:?} {c:?} uv {:?}",
                    m.uvs.get(i)
                );
            }
            continue;
        }
        if m.normals
            .get(i)
            .is_some_and(|n| n.iter().any(|x| x.dot(g) <= 0.0))
        {
            flipped += 1;
        }
    }
    (degenerate, flipped)
}

fn summary(m: &Mesh) -> String {
    let r = m.edge_report();
    let (degenerate, flipped) = bad_triangles(m);
    let dev = m
        .stats
        .max_deviation
        .map_or("-".to_string(), |d| format!("{d:.3e}"));
    format!(
        "V {} T {} degenerate {degenerate} flipped {flipped} open {} nonmanifold {} inconsistent {} failures {} unconverged {} split {} ears {} inverted {} grid backend {} (failed {}) dev {dev}",
        m.used_vertices(),
        m.triangles.len(),
        r.open.len(),
        r.non_manifold.len(),
        r.inconsistent.len(),
        m.failures.len(),
        m.unconverged_faces.len(),
        m.split_boundary_faces.len(),
        m.ear_clipped_faces.len(),
        m.inverted_edgeuses.len(),
        m.backend_faces.len(),
        m.backend_failed_faces.len(),
    )
}

/// The B-rep's kind and the mesh's open edges split into those along a B-rep edge with
/// one face (a sheet's boundary, expected) and the others (holes).
fn topology(b: &Brep, m: &Mesh) -> String {
    let solid = solidcore::topology::brep::p2::is_manifold_solid(b).unwrap_or(false);
    let mut lamina = std::collections::BTreeSet::new();
    for e in &m.edges {
        let faces = solidcore::topology::edge::faces(b, e.edge, None).map_or(0, |f| f.len());
        if faces <= 1 {
            for (x, y) in e.vertices.iter().zip(e.vertices.iter().skip(1)) {
                lamina.insert((*x.min(y), *x.max(y)));
            }
        }
    }
    let r = m.edge_report();
    let holes = r
        .open
        .iter()
        .filter(|[x, y]| !lamina.contains(&(*x, *y)))
        .count();
    format!(
        "{} | open on sheet boundary {} | open elsewhere {holes}",
        if solid {
            "solid"
        } else {
            "not a manifold solid"
        },
        r.open.len() - holes
    )
}

fn rel(a: f64, b: Option<f64>) -> String {
    b.map_or("-".into(), |b| format!("{:+.2e}", (a - b) / b))
}

fn main() -> Result<(), String> {
    let a = parse()?;
    let mut last: Option<(String, Mesh)> = None;
    for model in &a.models {
        let mut k = Kernel::default();
        let b = match load(&mut k, model) {
            Ok(b) => b,
            Err(e) => {
                println!("{model}: not read: {e}");
                continue;
            }
        };
        if a.bench {
            let (vol, area) = if a.no_ref {
                (None, None)
            } else {
                reference(&k, &b)
            };
            for r in [1e-2, 1e-3, 1e-4] {
                let mut o = a.opts.clone();
                o.chord_tolerance = Tolerance::Relative(r);
                o.measure_deviation = false;
                let t = Instant::now();
                let m = tessellate(&k, &b, &o).map_err(|e| e.to_string())?;
                let ms = t.elapsed().as_secs_f64() * 1e3;
                let d = solidcore_mesh::deviation(&k, &b, &m, m.stats.chord_tolerance)
                    .map_err(|e| e.to_string())?;
                // the same mesh again, measured per face (over = beyond chord plus the
                // corners' own edge/face gap)
                let mut mo = o.clone();
                mo.measure_deviation = true;
                let mm = tessellate(&k, &b, &mo).map_err(|e| e.to_string())?;
                println!(
                    "{model} rel {r:e} chord {:.3e}: {ms:.1} ms | {} | {} | measured max {:.3e} ({:.2}x chord), B-rep gap {:.1e}, over chord+gap {}, flipped against surface {} | volume {} area {}",
                    m.stats.chord_tolerance,
                    summary(&m),
                    topology(&b, &m),
                    d.max,
                    d.max / m.stats.chord_tolerance,
                    mm.stats.boundary_gap.unwrap_or(0.0),
                    mm.stats.triangles_over,
                    d.flipped,
                    rel(m.volume(), vol),
                    rel(m.area(), area),
                );
                if d.over > 0 && std::env::var("WORST").is_ok() {
                    if let Some(w) = d.worst {
                        println!(
                            "   worst triangle {w}: face {:?} points {:?} uvs {:?}",
                            m.triangle_faces.get(w),
                            m.triangles.get(w).and_then(|t| m.triangle_points(*t)),
                            m.uvs.get(w)
                        );
                    }
                }
            }
            continue;
        }
        if let Some(levels) = &a.lods {
            let opts: Vec<TessellateOptions> = levels
                .iter()
                .map(|r| {
                    let mut o = a.opts.clone();
                    o.chord_tolerance = Tolerance::Relative(*r);
                    o
                })
                .collect();
            // one nested refinement, level by level (each level's time without and
            // with its deviation measurement), against independent runs
            let (vol, area) = if a.no_ref {
                (None, None)
            } else {
                reference(&k, &b)
            };
            let mut nested = solidcore_mesh::Tessellator::new(&k, &b).map_err(|e| e.to_string())?;
            let mut ms = Vec::new();
            let (mut total, mut total_alone) = (0.0, 0.0);
            for (r, o) in levels.iter().zip(&opts) {
                let t = Instant::now();
                let m = nested.level(o).map_err(|e| e.to_string())?;
                let level_ms = t.elapsed().as_secs_f64() * 1e3;
                let t = Instant::now();
                let alone = tessellate(&k, &b, o).map_err(|e| e.to_string())?;
                let alone_ms = t.elapsed().as_secs_f64() * 1e3;
                let d = solidcore_mesh::deviation(&k, &b, &m, m.stats.chord_tolerance)
                    .map_err(|e| e.to_string())?;
                total += level_ms;
                total_alone += alone_ms;
                println!(
                    "{model} level rel {r:e}: {level_ms:.1} ms (alone {alone_ms:.1} ms, {} triangles) | {} | measured deviation {:.3e} | volume {} area {}",
                    alone.triangles.len(),
                    summary(&m),
                    d.max,
                    rel(m.volume(), vol),
                    rel(m.area(), area),
                );
                ms.push(m);
            }
            println!(
                "{model}: {} nested levels in {total:.1} ms, independent runs {total_alone:.1} ms",
                ms.len()
            );
            last = ms.into_iter().last().map(|m| (model.clone(), m));
            continue;
        }
        let t = Instant::now();
        let m = if a.edges_only {
            tessellate_edges(&k, &b, &a.opts)
        } else {
            tessellate(&k, &b, &a.opts)
        }
        .map_err(|e| e.to_string())?;
        let ms = t.elapsed().as_secs_f64() * 1e3;
        println!(
            "{model}: {ms:.1} ms | {} | volume {:.9} area {:.9}",
            summary(&m),
            m.volume(),
            m.area()
        );
        for f in &m.failures {
            println!("   failure {:?}: {}", f.entity, f.error);
        }
        if a.per_face {
            let mut count: BTreeMap<u32, usize> = BTreeMap::new();
            for f in &m.triangle_faces {
                *count.entry(f.index()).or_insert(0) += 1;
            }
            println!("   triangles per face: {count:?}");
            println!("   unconverged: {:?}", m.unconverged_faces);
            if !a.edges_only {
                let d = solidcore_mesh::deviation(&k, &b, &m, m.stats.chord_tolerance)
                    .map_err(|e| e.to_string())?;
                let mut flips: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
                for t in &d.flipped_triangles {
                    if let Some(f) = m.triangle_faces.get(*t) {
                        flips.entry(f.index()).or_default().push(*t);
                    }
                }
                for (f, ts) in &flips {
                    let first = ts.first().and_then(|t| {
                        Some((m.triangle_points(*m.triangles.get(*t)?)?, *m.uvs.get(*t)?))
                    });
                    println!(
                        "   face {f}: {} flipped against the surface, first {first:?}",
                        ts.len()
                    );
                }
            }
        }
        last = Some((model.clone(), m));
    }
    if let Some((model, m)) = &last {
        if let Some(p) = &a.stl {
            std::fs::write(p, stl(m)).map_err(|e| format!("{p}: {e}"))?;
            println!("STL: {p}");
        }
        if let Some(p) = &a.json {
            std::fs::write(p, json(model, m)).map_err(|e| format!("{p}: {e}"))?;
            println!("JSON: {p}");
        }
    }
    Ok(())
}

/// Binary STL: 80-byte header, triangle count, then normal, three corners and a zero
/// attribute per triangle (single precision).
fn stl(m: &Mesh) -> Vec<u8> {
    let mut out = Vec::with_capacity(84 + 50 * m.triangles.len());
    let mut header = [0u8; 80];
    for (h, c) in header.iter_mut().zip(b"solidcore-mesh binary STL".iter()) {
        *h = *c;
    }
    out.extend_from_slice(&header);
    let n = u32::try_from(m.triangles.len()).unwrap_or(u32::MAX);
    out.extend_from_slice(&n.to_le_bytes());
    for t in m.triangles.iter().take(n as usize) {
        let Some([a, b, c]) = m.triangle_points(*t) else {
            continue;
        };
        let g = (b - a).cross(c - a).unitize_or_z();
        for p in [g, a, b, c] {
            for x in [p.x, p.y, p.z] {
                out.extend_from_slice(&(x as f32).to_le_bytes());
            }
        }
        out.extend_from_slice(&[0, 0]);
    }
    out
}

fn num(x: f64) -> String {
    if x.is_finite() {
        format!("{x:?}")
    } else {
        "null".into()
    }
}

/// The viewer JSON (schema in the module docs).
fn json(model: &str, m: &Mesh) -> String {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in &m.positions {
        for ((l, h), x) in lo.iter_mut().zip(hi.iter_mut()).zip([p.x, p.y, p.z]) {
            *l = l.min(x);
            *h = h.max(x);
        }
    }
    let mut s = String::new();
    let _ = write!(
        s,
        "{{\"schema\":1,\"source\":\"solidcore-mesh\",\"model\":{:?},\"chord\":{},\"bbox\":[{},{},{},{},{},{}],\"positions\":[",
        model,
        num(m.stats.chord_tolerance),
        num(lo[0]),
        num(lo[1]),
        num(lo[2]),
        num(hi[0]),
        num(hi[1]),
        num(hi[2])
    );
    let pos: Vec<String> = m
        .positions
        .iter()
        .map(|p| format!("{},{},{}", num(p.x), num(p.y), num(p.z)))
        .collect();
    s.push_str(&pos.join(","));
    s.push_str("],\"triangles\":[");
    let tris: Vec<String> = m
        .triangles
        .iter()
        .map(|t| format!("{},{},{}", t[0], t[1], t[2]))
        .collect();
    s.push_str(&tris.join(","));
    s.push_str("],\"normals\":[");
    let nrm: Vec<String> = m
        .normals
        .iter()
        .flat_map(|n| {
            n.iter()
                .map(|v| format!("{},{},{}", num(v.x), num(v.y), num(v.z)))
        })
        .collect();
    s.push_str(&nrm.join(","));
    s.push_str("],\"triangle_faces\":[");
    let faces: Vec<String> = m
        .triangle_faces
        .iter()
        .map(|f| f.index().to_string())
        .collect();
    s.push_str(&faces.join(","));
    s.push_str("],\"edges\":[");
    let edges: Vec<String> = m
        .edges
        .iter()
        .map(|e| {
            let pts: Vec<String> = e
                .vertices
                .iter()
                .filter_map(|v| m.positions.get(*v as usize))
                .map(|p| format!("[{},{},{}]", num(p.x), num(p.y), num(p.z)))
                .collect();
            format!("[{}]", pts.join(","))
        })
        .collect();
    s.push_str(&edges.join(","));
    s.push_str("]}\n");
    s
}
