//! Checks of `solidcore-mesh` on the primitives and on upstream's test models: watertight
//! meshes, measured deviation within the tolerance, convergence of volume and area,
//! boundary vertices on their edge curves, outward normals, determinism, the per-face
//! API, levels of detail and the options.
//!
//! A primitive or model that solidcore cannot build or read yet (an `Unported`
//! callee) is skipped with a note, so the checks run wherever the kernel gets that far.

// test helpers outside `#[test]` functions unwrap and index like the tests do
#![allow(
    clippy::unwrap_used,
    clippy::indexing_slicing,
    clippy::field_reassign_with_default
)]

use std::collections::BTreeMap;

use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::geom::object::ObjectRef;
use solidcore::kernel::Kernel;
use solidcore::mass::integrator::properties::{brep_area, brep_volume, PropertyAttributes};
use solidcore::math::Point3;
use solidcore::primitive::api_primitives as api;
use solidcore::topology::brep::p1 as bl;
use solidcore::topology::brep::p10::TopologyEnv;
use solidcore::topology::brep::p7::read_from_file;
use solidcore::topology::database_io::FileType;
use solidcore::topology::Brep;
use solidcore::EntityRef;
use solidcore_mesh::*;

type Build = fn(&mut Kernel) -> solidcore::Result<Brep>;

/// The `primitives` example's nine solids with their analytic volume and area.
fn primitives() -> Vec<(&'static str, Build, f64, f64)> {
    use std::f64::consts::PI;
    let (rb, rt, h) = (5.0_f64, 2.0_f64, 20.0_f64);
    let slant = ((rb - rt) * (rb - rt) + h * h).sqrt();
    let (r, big, small) = (3.0_f64, 5.0_f64, 1.5_f64);
    vec![
        (
            "box",
            |k| api::create_box(k, Point3::new(0.0, 0.0, 0.0), 5.0, 10.0, 20.0),
            1000.0,
            700.0,
        ),
        (
            "cylinder",
            |k| api::create_cylinder(k, Point3::new(0.0, 0.0, 0.0), 5.0, 20.0),
            PI * rb * rb * h,
            2.0 * PI * rb * rb + 2.0 * PI * rb * h,
        ),
        (
            "cone",
            |k| api::create_cone(k, Point3::new(0.0, 0.0, 0.0), 5.0, 2.0, 20.0),
            PI * h / 3.0 * (rb * rb + rb * rt + rt * rt),
            PI * (rb + rt) * slant + PI * rb * rb + PI * rt * rt,
        ),
        (
            "pyramid",
            |k| api::create_pyramid(k, Point3::new(0.0, 0.0, 0.0), 5.0, 20.0),
            25.0 * 20.0 / 3.0,
            25.0 + 2.0 * 5.0 * (2.5 * 2.5 + 400.0_f64).sqrt(),
        ),
        (
            "sphere",
            |k| api::create_sphere(k, Point3::new(1.0, -2.0, 0.5), 3.0),
            4.0 / 3.0 * PI * r * r * r,
            4.0 * PI * r * r,
        ),
        (
            "sphere-no-pole",
            |k| api::create_sphere_no_pole(k, Point3::new(1.0, -2.0, 0.5), 3.0),
            4.0 / 3.0 * PI * r * r * r,
            4.0 * PI * r * r,
        ),
        (
            "partial-sphere",
            |k| api::create_partial_sphere(k, Point3::new(1.0, -2.0, 0.5), 3.0, 0.0, 90.0),
            PI * r * r * r / 3.0,
            2.0 * PI * r * r,
        ),
        (
            "torus",
            |k| api::create_torus(k, Point3::new(0.0, 1.0, -1.0), 5.0, 1.5),
            2.0 * PI * PI * big * small * small,
            4.0 * PI * PI * big * small,
        ),
        (
            "partial-torus",
            |k| api::create_partial_torus(k, Point3::new(0.0, 1.0, -1.0), 5.0, 1.5, 0.0, 90.0),
            PI * PI * big * small * small / 2.0,
            PI * PI * big * small + 2.0 * PI * small * small,
        ),
    ]
}

fn build(name: &str, f: Build) -> Option<(Kernel, Brep)> {
    let mut k = Kernel::default();
    match f(&mut k) {
        Ok(b) => Some((k, b)),
        Err(e) => {
            eprintln!("skip {name}: does not build yet ({e})");
            None
        }
    }
}

/// The checks every closed solid's mesh must pass; returns the measured deviation.
fn check_closed(name: &str, k: &Kernel, b: &Brep, m: &Mesh) -> f64 {
    assert!(m.failures.is_empty(), "{name}: failures {:?}", m.failures);
    let r = m.edge_report();
    assert!(
        r.is_watertight(),
        "{name}: open {} non-manifold {} inconsistent {}",
        r.open.len(),
        r.non_manifold.len(),
        r.inconsistent.len()
    );
    // outward, non-degenerate triangles: each corner normal agrees with the facet
    let size = m.stats.model_size;
    for (t, n) in m.triangles.iter().zip(&m.normals) {
        let [a, p, c] = m.triangle_points(*t).unwrap();
        let g = (p - a).cross(c - a);
        assert!(
            g.length() > 1e-14 * size * size,
            "{name}: degenerate triangle {t:?}"
        );
        let g = g / g.length();
        for x in n {
            assert!(
                g.dot(*x) > 0.0,
                "{name}: corner normal against its triangle {t:?}"
            );
        }
    }
    // every edge polyline point lies on its curve: interior points within float noise, ends within
    // the B-rep's own vertex gap
    for e in &m.edges {
        if e.degenerate {
            continue;
        }
        let rec = b.edge(e.edge).unwrap();
        let curve = rec.curve.as_ref().unwrap();
        let last = e.params.len() - 1;
        for (i, (t, v)) in e.params.iter().zip(&e.vertices).enumerate() {
            let d = curve
                .evaluate_point(k, *t)
                .unwrap()
                .distance(m.positions[*v as usize]);
            let allowed = if i == 0 || i == last {
                rec.tolerance.max(1e-9 * size) * 10.0
            } else {
                1e-12 * size
            };
            assert!(d <= allowed, "{name}: edge point {d:e} off its curve");
        }
    }
    let d = deviation(k, b, m, m.stats.chord_tolerance).unwrap();
    assert_eq!(d.flipped, 0, "{name}: triangles wound against the surface");
    // within the chord tolerance, plus the B-rep's own edge/face gap (measured runs)
    let gap = m.stats.boundary_gap.unwrap_or(0.0);
    assert!(
        d.max <= 1.05 * (m.stats.chord_tolerance + gap),
        "{name}: measured deviation {} over chord {} + gap {gap}",
        d.max,
        m.stats.chord_tolerance
    );
    d.max
}

#[test]
fn primitives_watertight_within_tolerance_and_converging() {
    for (name, f, volume, area) in primitives() {
        let Some((k, b)) = build(name, f) else {
            continue;
        };
        let mut errors = Vec::new();
        for rel in [1e-2, 1e-3, 1e-4] {
            let opts = TessellateOptions::default().with_relative_chord(rel);
            let m = tessellate(&k, &b, &opts).unwrap();
            check_closed(name, &k, &b, &m);
            let (ev, ea) = ((m.volume() - volume) / volume, (m.area() - area) / area);
            eprintln!(
                "{name} rel {rel:e}: V {} T {} volume {ev:+.2e} area {ea:+.2e}",
                m.used_vertices(),
                m.triangles.len()
            );
            errors.push((ev.abs(), ea.abs()));
        }
        // a polyhedron is exact; curved solids converge as the chord shrinks
        if name == "box" || name == "pyramid" {
            assert!(
                errors.iter().all(|(v, a)| *v < 1e-12 && *a < 1e-12),
                "{name}: {errors:?}"
            );
        } else {
            assert!(
                errors[2].0 < errors[0].0 && errors[2].1 < errors[0].1,
                "{name}: {errors:?}"
            );
            assert!(
                errors[2].0 < 2e-3 && errors[2].1 < 2e-3,
                "{name}: {errors:?}"
            );
        }
    }
}

#[test]
fn deterministic_bit_for_bit() {
    for (name, f, _, _) in primitives() {
        let Some((k, b)) = build(name, f) else {
            continue;
        };
        let opts = TessellateOptions::default();
        let a = tessellate(&k, &b, &opts).unwrap();
        let c = tessellate(&k, &b, &opts).unwrap();
        assert_eq!(a, c, "{name}: two runs differ");
    }
}

#[test]
fn per_face_api_matches_tessellate() {
    for (name, f, _, _) in primitives() {
        let Some((k, b)) = build(name, f) else {
            continue;
        };
        let opts = TessellateOptions::default();
        let whole = tessellate(&k, &b, &opts).unwrap();
        // the protocol FaceTessellation documents
        let mut edges = tessellate_edge_set(&k, &b, &opts).unwrap();
        let ids = bl::faces(&b, None).unwrap();
        let mut faces: Vec<_> = ids
            .iter()
            .map(|f| {
                let mut t = FaceTessellation::new(&k, &b, *f, &edges).unwrap();
                t.refine(&k, &b, &edges, &opts).unwrap();
                t
            })
            .collect();
        for _ in 0..MAX_BOUNDARY_ROUNDS {
            let mut requests: Vec<SegmentSplit> = faces
                .iter()
                .flat_map(|t| t.boundary_requests().iter().copied())
                .collect();
            requests.sort_unstable();
            requests.dedup();
            if requests.is_empty() {
                break;
            }
            let changed = edges.split_segments(&k, &b, &requests).unwrap();
            if changed.is_empty() {
                break;
            }
            for t in &mut faces {
                if changed.iter().any(|e| t.uses_edge(*e)) {
                    t.refine(&k, &b, &edges, &opts).unwrap();
                }
            }
        }
        let faces: Vec<_> = ids
            .iter()
            .zip(&faces)
            .map(|(f, t)| (*f, t.mesh(&k, &b, &edges, &opts)))
            .collect();
        let joined = assemble(edges, faces, &opts, Some((&k, &b)));
        assert_eq!(whole, joined, "{name}: per-face assembly differs");
    }
}

#[test]
fn levels_of_detail_nest() {
    for (name, f, _, _) in primitives() {
        let Some((k, b)) = build(name, f) else {
            continue;
        };
        let levels: Vec<_> = [1e-2, 3e-3, 1e-3, 3e-4]
            .iter()
            .map(|r| TessellateOptions::default().with_relative_chord(*r))
            .collect();
        let t = std::time::Instant::now();
        let ms = tessellate_lods(&k, &b, &levels).unwrap();
        eprintln!(
            "{name}: {} levels in {:.1} ms",
            ms.len(),
            t.elapsed().as_secs_f64() * 1e3
        );
        let mut prev: Option<&Mesh> = None;
        for m in &ms {
            check_closed(name, &k, &b, m);
            let dev = m.stats.max_deviation.unwrap();
            if let Some(p) = prev {
                // the coarser level's vertices keep their ids and positions
                assert!(m.positions.len() >= p.positions.len());
                assert_eq!(
                    &m.positions[..p.positions.len()],
                    &p.positions[..],
                    "{name}: vertex ids moved"
                );
                // every edge polyline point of the coarser level is still there
                for (a, c) in p.edges.iter().zip(&m.edges) {
                    assert!(
                        a.vertices.iter().all(|v| c.vertices.contains(v)),
                        "{name}: edge point lost"
                    );
                }
                let pd = p.stats.max_deviation.unwrap();
                assert!(
                    dev <= pd || dev <= m.stats.chord_tolerance,
                    "{name}: deviation grew {pd} -> {dev}"
                );
            }
            eprintln!(
                "{name} level chord {:.3e}: V {} T {} deviation {dev:.3e}",
                m.stats.chord_tolerance,
                m.used_vertices(),
                m.triangles.len()
            );
            prev = Some(m);
        }
    }
}

#[test]
fn edges_only_mode() {
    for (name, f, _, _) in primitives() {
        let Some((k, b)) = build(name, f) else {
            continue;
        };
        let m = tessellate_edges(&k, &b, &TessellateOptions::default()).unwrap();
        assert!(m.triangles.is_empty());
        assert_eq!(m.edges.len(), bl::edges(&b, None).unwrap().len(), "{name}");
        let full = tessellate(&k, &b, &TessellateOptions::default()).unwrap();
        // the full mesh's polylines are these, with at most points added where faces
        // asked for them
        for (a, c) in m.edges.iter().zip(&full.edges) {
            assert!(a.params.iter().all(|t| c.params.contains(t)), "{name}");
        }
    }
}

#[test]
fn seams_and_poles() {
    // the sphere's pole vertices are fan centres; the torus has seams in both directions
    for (name, f, _, _) in primitives() {
        if !matches!(name, "sphere" | "torus" | "cylinder" | "partial-sphere") {
            continue;
        }
        let Some((k, b)) = build(name, f) else {
            continue;
        };
        for q in [Quality::Draft, Quality::Normal, Quality::Fine] {
            let m = tessellate(&k, &b, &TessellateOptions::preset(q)).unwrap();
            check_closed(name, &k, &b, &m);
            if name == "sphere" {
                // the poles are mesh vertices that fans of triangles meet at
                for pole in [Point3::new(1.0, -2.0, 3.5), Point3::new(1.0, -2.0, -2.5)] {
                    let (i, d) = m
                        .positions
                        .iter()
                        .enumerate()
                        .map(|(i, p)| (i, p.distance(pole)))
                        .min_by(|a, b| a.1.total_cmp(&b.1))
                        .unwrap();
                    assert!(d < 1e-9, "sphere: no vertex at the pole {pole:?} ({d})");
                    let fan = m
                        .triangles
                        .iter()
                        .filter(|t| t.contains(&(i as u32)))
                        .count();
                    assert!(fan >= 4, "sphere: pole fan of {fan}");
                }
            }
        }
    }
}

/// A backend that asks for a million lines each way on every face, unsorted, with
/// junk outside the box.
struct Hostile;

impl backend::FaceBackend for Hostile {
    fn name(&self) -> &'static str {
        "hostile"
    }

    fn grid_lines(
        &self,
        _k: &Kernel,
        face: &backend::BackendFace<'_>,
    ) -> Option<solidcore::Result<(Vec<f64>, Vec<f64>)>> {
        let b = face.bbox;
        let lines = |lo: f64, hi: f64| -> Vec<f64> {
            (0..1_000_000)
                .rev()
                .map(|i| lo + (hi - lo) * f64::from(i) / 1e6)
                .chain([f64::NAN, f64::INFINITY, lo - 1.0, hi + 1.0])
                .collect()
        };
        Some(Ok((
            lines(b.u_min(), b.u_max()),
            lines(b.v_min(), b.v_max()),
        )))
    }
}

#[test]
fn backend_lines_are_bounded() {
    let Some((k, b)) = build("box", primitives()[0].1) else {
        return;
    };
    let mut o = TessellateOptions::default();
    o.max_triangles_per_face = 5_000;
    o.backend = Some(std::sync::Arc::new(Hostile));
    let m = tessellate(&k, &b, &o).unwrap();
    assert!(m.failures.is_empty(), "{:?}", m.failures);
    assert!(m.backend_faces.len() == 6, "{:?}", m.backend_faces);
    let mut per_face: BTreeMap<u32, usize> = BTreeMap::new();
    for f in &m.triangle_faces {
        *per_face.entry(f.index()).or_insert(0) += 1;
    }
    // the budget, plus what the boundary and one last insertion add
    assert!(per_face.values().all(|n| *n <= 2 * 5_000), "{per_face:?}");
    let _ = check_closed("box with a hostile backend", &k, &b, &m);
}

#[test]
fn options_are_honoured() {
    let Some((k, b)) = build("cylinder", primitives()[1].1) else {
        return;
    };
    // invalid options are refused
    let mut bad = TessellateOptions::default();
    bad.angle_tolerance = 0.0;
    assert!(tessellate(&k, &b, &bad).is_err());
    let bad = TessellateOptions::default().with_chord(-1.0);
    assert!(tessellate(&k, &b, &bad).is_err());

    // presets get finer
    let counts: Vec<usize> = [Quality::Draft, Quality::Normal, Quality::Fine]
        .iter()
        .map(|q| {
            tessellate(&k, &b, &TessellateOptions::preset(*q))
                .unwrap()
                .triangles
                .len()
        })
        .collect();
    assert!(counts[0] < counts[1] && counts[1] < counts[2], "{counts:?}");

    // output toggles
    let mut o = TessellateOptions::default();
    o.output = OutputOptions {
        normals: false,
        uvs: false,
        face_ids: false,
        edge_polylines: false,
    };
    let m = tessellate(&k, &b, &o).unwrap();
    assert!(
        m.normals.is_empty()
            && m.uvs.is_empty()
            && m.triangle_faces.is_empty()
            && m.edges.is_empty()
    );
    assert!(!m.triangles.is_empty());

    // edge segment counts and lengths
    let mut o = TessellateOptions::default();
    o.edge.min_segments = 7;
    let m = tessellate(&k, &b, &o).unwrap();
    assert!(m
        .edges
        .iter()
        .filter(|e| !e.degenerate)
        .all(|e| e.params.len() >= 8));
    let mut o = TessellateOptions::default().with_relative_chord(1e-5);
    o.edge.max_segments = 10;
    let m = tessellate(&k, &b, &o).unwrap();
    assert!(m.edges.iter().all(|e| e.params.len() <= 11));
    let mut o = TessellateOptions::default();
    o.max_edge_length = Some(1.0);
    let m = tessellate(&k, &b, &o).unwrap();
    for t in &m.triangles {
        let [a, p, c] = m.triangle_points(*t).unwrap();
        assert!(a.distance(p).max(p.distance(c)).max(c.distance(a)) <= 1.0 + 1e-9);
    }

    // a per-face override and a per-edge hook
    let base = tessellate(&k, &b, &TessellateOptions::default()).unwrap();
    let wall = base.triangle_faces[0];
    let per_face = |m: &Mesh, f| m.triangle_faces.iter().filter(|x| **x == f).count();
    let mut o = TessellateOptions::default();
    o.face_overrides.insert(
        wall,
        FaceOverride {
            chord_tolerance: Some(base.stats.chord_tolerance / 10.0),
            ..FaceOverride::default()
        },
    );
    let finer = tessellate(&k, &b, &o).unwrap();
    assert!(per_face(&finer, wall) > per_face(&base, wall));
    check_closed("cylinder override", &k, &b, &finer);
    let mut o = TessellateOptions::default();
    o.edge_hook = Some(std::sync::Arc::new(|_| {
        Some(EdgeOverride {
            min_segments: Some(40),
            ..EdgeOverride::default()
        })
    }));
    let m = tessellate(&k, &b, &o).unwrap();
    assert!(m
        .edges
        .iter()
        .filter(|e| !e.degenerate)
        .all(|e| e.params.len() >= 41));
    check_closed("cylinder edge hook", &k, &b, &m);
}

// ---- upstream test models ------------------------------------------------------------

const CORPUS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../upstream/TestFiles/pt_TestFiles"
);

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

fn read(k: &mut Kernel, rel: &str) -> Option<Brep> {
    let dir = std::env::var("SOLIDCORE_CORPUS").unwrap_or_else(|_| CORPUS.to_string());
    let Ok(bytes) = std::fs::read(format!("{dir}/{rel}")) else {
        eprintln!("skip {rel}: corpus not found");
        return None;
    };
    let mut b = Brep::new_with_infinite_region().ok()?;
    let mut attributes = AttributeStore::default();
    let mut records = BTreeMap::new();
    let mut env = TopologyEnv {
        all_attributes: &[],
        attributes: &mut attributes,
        key_of: &key_of,
        records: &mut records,
        next_key: 1 << 62,
    };
    match read_from_file(k, &mut b, bytes, FileType::Ascii, false, &mut env) {
        Ok(()) => Some(b),
        Err(e) => {
            eprintln!("skip {rel}: not read yet ({e})");
            None
        }
    }
}

/// Manifold solids from upstream's test models.
const MODELS: [&str; 14] = [
    "Fillet/vase3.smb",
    "Fillet/crown.smb",
    "Defeature/Platte.smb",
    "Fillet/letter_g.smb",
    "Fillet/bolt.smb",
    "Fillet/cube_sphere2.smb",
    "Fillet/torii_fillet.smb",
    "Fillet/Heart.smb",
    "Fillet/box_diff_3cyls.smb",
    "Fillet/nurb_box1.smb",
    "Solids/revolution.smb",
    "Solids/sphere.smb",
    "Mouse.smb",
    "toruswvy.smb",
];

#[test]
fn corpus_models_watertight_within_tolerance() {
    for rel in MODELS {
        let mut k = Kernel::default();
        let Some(mut b) = read(&mut k, rel) else {
            continue;
        };
        let store = AttributeStore::new();
        let key = |_| ObjectKey(0);
        let attrs = PropertyAttributes {
            store: &store,
            key_of: &key,
        };
        let volume = brep_volume(&k, &mut b, 1e-6, &attrs).unwrap();
        let area = brep_area(&k, &mut b, 1e-6).unwrap();
        let mut errors = Vec::new();
        for r in [1e-2, 1e-3, 1e-4] {
            let mut o = TessellateOptions::default().with_relative_chord(r);
            o.measure_deviation = true;
            let m = tessellate(&k, &b, &o).unwrap();
            check_closed(rel, &k, &b, &m);
            errors.push(((m.volume() - volume) / volume, (m.area() - area) / area));
        }
        eprintln!("{rel}: relative volume/area errors at 1e-2, 1e-3, 1e-4: {errors:?}");
        assert!(
            errors[2].0.abs() < 1e-3 && errors[2].1.abs() < 1e-3,
            "{rel}: {errors:?}"
        );
    }
}
