//! STEP topology as data, one entity kind and one refusal at a time: small solids written here
//! entity by entity (boxes, a box with a void, a cone to its apex, sheets, loose faces), each
//! read and checked against what the file writes: counts, indices, flags and the geometry's
//! data. Hostile files are errors, never panics. Expected values are the file's own numbers.
//!
//! Origin: SolidCraft's `crates/step/src/topology_tests.rs` at `a5b71ca2` (its writer, the box,
//! the void, the surface model and the hostile edits), without the parts that test the schema-4
//! record (pcurves, intervals, seams, sides measured on the surfaces), which are not in the
//! record's reader.

use std::collections::BTreeMap;
use std::fmt::Write as _;

use super::*;
use crate::step::geom::{Placement, Vector};
use crate::step::p21;

type V = [f64; 3];

fn real(x: f64) -> String {
    format!("{x:.17E}")
}
fn sub3(a: V, b: V) -> V {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn neg(a: V) -> V {
    [-a[0], -a[1], -a[2]]
}
fn b(x: bool) -> &'static str {
    if x { ".T." } else { ".F." }
}
fn list(ids: &[u64]) -> String {
    ids.iter()
        .map(|i| format!("#{i}"))
        .collect::<Vec<_>>()
        .join(",")
}

/// Writes a STEP data section one entity at a time.
#[derive(Default)]
struct W {
    s: String,
    n: u64,
}

impl W {
    fn add(&mut self, body: impl AsRef<str>) -> u64 {
        self.n += 1;
        let _ = writeln!(self.s, "#{}={};", self.n, body.as_ref());
        self.n
    }
    fn pt(&mut self, p: V) -> u64 {
        self.add(format!(
            "CARTESIAN_POINT('',({},{},{}))",
            real(p[0]),
            real(p[1]),
            real(p[2])
        ))
    }
    fn dir(&mut self, d: V) -> u64 {
        self.add(format!(
            "DIRECTION('',({},{},{}))",
            real(d[0]),
            real(d[1]),
            real(d[2])
        ))
    }
    fn axis(&mut self, o: V, z: V, x: V) -> u64 {
        let (o, z, x) = (self.pt(o), self.dir(z), self.dir(x));
        self.add(format!("AXIS2_PLACEMENT_3D('',#{o},#{z},#{x})"))
    }
    fn vertex(&mut self, p: V) -> u64 {
        let p = self.pt(p);
        self.add(format!("VERTEX_POINT('',#{p})"))
    }
    /// A line edge from `a` to `b`; with `same_sense` false its LINE runs from b to a.
    fn line_edge(
        &mut self,
        (va, pa): (u64, V),
        (vb, pb): (u64, V),
        same_sense: bool,
        magnitude: f64,
    ) -> u64 {
        let (o, d) = if same_sense {
            (pa, sub3(pb, pa))
        } else {
            (pb, sub3(pa, pb))
        };
        let (o, d) = (self.pt(o), self.dir(d));
        let v = self.add(format!("VECTOR('',#{d},{})", real(magnitude)));
        let l = self.add(format!("LINE('',#{o},#{v})"));
        self.add(format!("EDGE_CURVE('',#{va},#{vb},#{l},{})", b(same_sense)))
    }
    fn circle(&mut self, c: V, z: V, x: V, r: f64) -> u64 {
        let a = self.axis(c, z, x);
        self.add(format!("CIRCLE('',#{a},{})", real(r)))
    }
    fn oe(&mut self, edge: u64, orientation: bool) -> u64 {
        self.add(format!("ORIENTED_EDGE('',*,*,#{edge},{})", b(orientation)))
    }
    fn bound(&mut self, oes: &[u64], outer: bool, orientation: bool) -> u64 {
        let l = self.add(format!("EDGE_LOOP('',({}))", list(oes)));
        self.add(format!(
            "{}('',#{l},{})",
            if outer {
                "FACE_OUTER_BOUND"
            } else {
                "FACE_BOUND"
            },
            b(orientation)
        ))
    }
    fn plane(&mut self, o: V, z: V, x: V) -> u64 {
        let a = self.axis(o, z, x);
        self.add(format!("PLANE('',#{a})"))
    }
    fn face(&mut self, bounds: &[u64], surface: u64, same_sense: bool) -> u64 {
        self.add(format!(
            "ADVANCED_FACE('',({}),#{surface},{})",
            list(bounds),
            b(same_sense)
        ))
    }
    /// The file: `unit` lengths, radians, an uncertainty of 1e-7 `unit`, every item in one
    /// representation.
    fn file_in(mut self, items: &[u64], unit: &str) -> Exchange {
        let l = self.add(unit);
        let a = self.add("(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.))");
        let s = self.add("(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT())");
        let u = self.add(format!(
            "UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07),#{l},'distance_accuracy_value','')"
        ));
        let c = self.add(format!(
            "(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#{u}))GLOBAL_UNIT_ASSIGNED_CONTEXT((#{l},#{a},#{s}))REPRESENTATION_CONTEXT('',''))"
        ));
        self.add(format!(
            "ADVANCED_BREP_SHAPE_REPRESENTATION('',({}),#{c})",
            list(items)
        ));
        parse(&self.s)
    }
    /// [`Self::file_in`] millimetres.
    fn file(self, items: &[u64]) -> Exchange {
        self.file_in(
            items,
            "(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.))",
        )
    }
}

fn parse(data: &str) -> Exchange {
    let t = format!(
        "ISO-10303-21;\nHEADER;FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));ENDSEC;\nDATA;\n{data}ENDSEC;\nEND-ISO-10303-21;\n"
    );
    match p21::parse(&t) {
        Ok(ex) => ex,
        Err(e) => panic!("{e}\n{t}"),
    }
}

/// How a box is written.
#[derive(Clone, Copy, Default)]
struct BoxStyle {
    /// Every other face on a plane facing into the box, with `same_sense` `.F.`
    flip_faces: bool,
    /// Every other edge with `same_sense` `.F.` (its line running from its end to its start).
    flip_edges: bool,
    /// Every face's bound written the other way round with orientation .F.
    flip_bounds: bool,
}

/// What [`box_faces`] wrote, to check the reader against.
#[derive(Default)]
struct Written {
    /// Each face's bound's loop as the face uses it: the corners, by index 0..8.
    loops: Vec<Vec<usize>>,
    /// Each vertex's point, by corner.
    corners: Vec<V>,
}

/// The box from `lo` to `hi`: its faces (outward normals), in the order a, side.
fn box_faces(w: &mut W, lo: V, hi: V, style: BoxStyle, out: &mut Written) -> Vec<u64> {
    let corner = |i: usize| {
        [
            if i & 1 == 0 { lo[0] } else { hi[0] },
            if i & 2 == 0 { lo[1] } else { hi[1] },
            if i & 4 == 0 { lo[2] } else { hi[2] },
        ]
    };
    let vs: Vec<(u64, V)> = (0..8).map(|i| (w.vertex(corner(i)), corner(i))).collect();
    out.corners = (0..8).map(corner).collect();
    let mut edges: BTreeMap<(usize, usize), u64> = BTreeMap::new();
    let mut faces = Vec::new();
    let mut count = 0;
    for a in 0..3 {
        for side in [false, true] {
            let (bx, cx) = ((a + 1) % 3, (a + 2) % 3);
            let bit = |axis: usize| 1usize << axis;
            let base = if side { bit(a) } else { 0 };
            let mut ring = [
                base,
                base | bit(bx),
                base | bit(bx) | bit(cx),
                base | bit(cx),
            ];
            if !side {
                ring.reverse();
            }
            out.loops.push(ring.to_vec());
            let mut oes = Vec::new();
            for k in 0..4 {
                let (p, q) = (ring[k], ring[(k + 1) % 4]);
                let key = (p.min(q), p.max(q));
                let e = if let Some(e) = edges.get(&key) {
                    *e
                } else {
                    let n = edges.len();
                    let flip = style.flip_edges && n % 2 == 1;
                    let magnitude = if n.is_multiple_of(3) { 1.0 } else { 2.5 };
                    let e = w.line_edge(vs[key.0], vs[key.1], !flip, magnitude);
                    edges.insert(key, e);
                    e
                };
                // Written the other way round: each edge the other way, from the last.
                oes.push((e, p == key.0));
            }
            let bounds = if style.flip_bounds {
                let rev: Vec<u64> = oes.iter().rev().map(|&(e, o)| w.oe(e, !o)).collect();
                vec![w.bound(&rev, true, false)]
            } else {
                let fwd: Vec<u64> = oes.iter().map(|&(e, o)| w.oe(e, o)).collect();
                vec![w.bound(&fwd, true, true)]
            };
            let mut normal = [0.0; 3];
            normal[a] = if side { 1.0 } else { -1.0 };
            let mut x = [0.0; 3];
            x[bx] = 1.0;
            let mut center = [
                0.5 * (lo[0] + hi[0]),
                0.5 * (lo[1] + hi[1]),
                0.5 * (lo[2] + hi[2]),
            ];
            center[a] = if side { hi[a] } else { lo[a] };
            let flip = style.flip_faces && count % 2 == 1;
            let plane = w.plane(center, if flip { neg(normal) } else { normal }, x);
            faces.push(w.face(&bounds, plane, !flip));
            count += 1;
        }
    }
    faces
}

fn closed_shell(w: &mut W, faces: &[u64]) -> u64 {
    w.add(format!("CLOSED_SHELL('',({}))", list(faces)))
}

fn one_solid(ex: &Exchange) -> StepSolid {
    let r = read_solids(ex);
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    assert_eq!(r.solids.len(), 1);
    r.solids.into_iter().next().unwrap()
}

/// Shells, faces, edges, vertices.
fn counts(s: &StepSolid) -> [usize; 4] {
    [
        s.shells.len() + s.voids.len(),
        s.faces.len(),
        s.edges.len(),
        s.vertices.len(),
    ]
}

fn cube(style: BoxStyle, out: &mut Written) -> (Exchange, Vec<u64>) {
    let mut w = W::default();
    let faces = box_faces(&mut w, [0.0; 3], [10.0, 20.0, 30.0], style, out);
    let shell = closed_shell(&mut w, &faces);
    let solid = w.add(format!("MANIFOLD_SOLID_BREP('cube',#{shell})"));
    (w.file(&[solid]), faces)
}

#[test]
fn a_cube_is_one_closed_shell_of_six_faces_twelve_edges_eight_vertices() {
    for style in [
        BoxStyle::default(),
        BoxStyle {
            flip_faces: true,
            ..Default::default()
        },
        BoxStyle {
            flip_edges: true,
            ..Default::default()
        },
        BoxStyle {
            flip_bounds: true,
            ..Default::default()
        },
        BoxStyle {
            flip_faces: true,
            flip_edges: true,
            flip_bounds: true,
        },
    ] {
        let mut written = Written::default();
        let (ex, face_ids) = cube(style, &mut written);
        let s = one_solid(&ex);
        assert_eq!(s.kind, SolidKind::ManifoldSolid);
        assert_eq!(s.name, "cube");
        assert_eq!(s.product, None);
        assert_eq!(s.units, Units { len: 1.0, ang: 1.0 });
        assert_eq!(s.tolerance, Some(1e-7 * 1.0));
        assert_eq!(counts(&s), [1, 6, 12, 8]);
        let shell = &s.shells[0];
        assert!(shell.declared_closed && shell.closed && shell.orientation);
        assert!(s.voids.is_empty());
        // Faces in the shell's order, each the file's face, not turned over.
        let listed: Vec<u64> = shell.faces.iter().map(|f| f.entity).collect();
        assert_eq!(listed, face_ids);
        for (k, of) in shell.faces.iter().enumerate() {
            assert_eq!(of.face as usize, k);
            assert!(of.orientation);
            let f = &s.faces[k];
            assert_eq!(f.entity, face_ids[k]);
            assert_eq!(f.kind, FaceKind::Advanced);
            // Same_sense as written: every other face .F. when its plane faces in.
            assert_eq!(f.same_sense, !(style.flip_faces && k % 2 == 1));
            assert_eq!(f.bounds.len(), 1);
            assert_eq!(f.loop_order(), [0]);
            let bound = &f.bounds[0];
            assert!(bound.outer);
            assert_eq!(bound.orientation, !style.flip_bounds);
            // The loop as the face uses it runs round the ring that was written.
            let corners: Vec<usize> = bound
                .face_loop()
                .iter()
                .map(|oe| {
                    let v = s.ends(*oe).unwrap()[0];
                    let p = s.vertices[v as usize].point;
                    written.corners.iter().position(|c| *c == p).unwrap()
                })
                .collect();
            assert_eq!(corners, written.loops[k]);
            let Surface::Plane { position } = &*f.surface else {
                panic!("{:?}", f.surface)
            };
            assert_eq!(
                position
                    .axis
                    .map(|z| z.iter().map(|x| x.abs()).sum::<f64>()),
                Some(1.0)
            );
        }
        // Edges: lines between their vertices as written, magnitude as written, each used
        // twice, once each way.
        let uses = s.edge_uses().unwrap();
        for (e, edge) in s.edges.iter().enumerate() {
            let (a, z) = (
                s.vertices[edge.start as usize].point,
                s.vertices[edge.end as usize].point,
            );
            let Curve::Line { origin, vector } = &*edge.curve else {
                panic!("{:?}", edge.curve)
            };
            let (o, d) = if edge.same_sense {
                (a, sub3(z, a))
            } else {
                (z, sub3(a, z))
            };
            assert_eq!(*origin, o);
            assert_eq!(vector.orientation, d);
            assert!(vector.magnitude == 1.0 || vector.magnitude == 2.5);
            assert_eq!(vector.scale, [1.0; 3]);
            assert_eq!(
                edge.curve_entity,
                ex.get(edge.entity).unwrap().params()[3]
                    .as_ref_id()
                    .unwrap()
            );
            assert_eq!(uses[e].len(), 2);
            assert_ne!(uses[e][0].orientation, uses[e][1].orientation);
        }
        let against = s.edges.iter().filter(|e| !e.same_sense).count();
        assert_eq!(against, if style.flip_edges { 6 } else { 0 });
    }
}

#[test]
fn lengths_and_the_uncertainty_are_in_millimetres() {
    let mut w = W::default();
    let mut written = Written::default();
    let faces = box_faces(
        &mut w,
        [0.0; 3],
        [1.0, 2.0, 3.0],
        BoxStyle::default(),
        &mut written,
    );
    let shell = closed_shell(&mut w, &faces);
    let solid = w.add(format!("MANIFOLD_SOLID_BREP('',#{shell})"));
    let inch = "(CONVERSION_BASED_UNIT('INCH',#9000)LENGTH_UNIT()NAMED_UNIT(*));\
#9000=LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4),#9001);\
#9001=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.))";
    let s = one_solid(&w.file_in(&[solid], inch));
    assert_eq!(s.units.len, 25.4);
    assert_eq!(s.tolerance, Some(1e-7 * 25.4));
    let mut points: Vec<V> = s.vertices.iter().map(|v| v.point).collect();
    points.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut want: Vec<V> = written
        .corners
        .iter()
        .map(|c| [c[0] * 25.4, c[1] * 25.4, c[2] * 25.4])
        .collect();
    want.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert_eq!(points, want);
}

#[test]
fn a_solid_with_a_void_keeps_its_oriented_shell() {
    let mut w = W::default();
    let mut written = Written::default();
    let outer = box_faces(
        &mut w,
        [0.0; 3],
        [10.0; 3],
        BoxStyle::default(),
        &mut written,
    );
    let outer = closed_shell(&mut w, &outer);
    let inner = box_faces(
        &mut w,
        [4.0; 3],
        [6.0; 3],
        BoxStyle::default(),
        &mut written,
    );
    let inner = closed_shell(&mut w, &inner);
    let void = w.add(format!("ORIENTED_CLOSED_SHELL('',*,#{inner},.F.)"));
    let solid = w.add(format!("BREP_WITH_VOIDS('',#{outer},(#{void}))"));
    let s = one_solid(&w.file(&[solid]));
    assert_eq!(s.kind, SolidKind::BrepWithVoids);
    assert_eq!(counts(&s), [2, 12, 24, 16]);
    assert_eq!(s.shells.len(), 1);
    assert_eq!(s.shells[0].entity, Some(outer));
    assert!(s.shells[0].orientation && s.shells[0].closed);
    let v = &s.voids[0];
    assert_eq!(v.entity, Some(void));
    assert!(v.declared_closed && v.closed);
    assert!(!v.orientation);
    let faces: Vec<u32> = v.faces.iter().map(|f| f.face).collect();
    assert_eq!(faces, (6..12).collect::<Vec<u32>>());
    assert!(v.faces.iter().all(|f| f.orientation));
}

#[test]
fn a_cone_reaches_its_apex_through_a_vertex_loop() {
    let mut w = W::default();
    let v = w.vertex([5.0, 0.0, 0.0]);
    let apex = w.vertex([0.0, 0.0, 10.0]);
    let k = w.circle([0.0; 3], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], 5.0);
    let e = w.add(format!("EDGE_CURVE('',#{v},#{v},#{k},.T.)"));
    // The cone widens down -z from its apex at z = 10 to radius 5 at z = 0.
    let a = w.axis([0.0; 3], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]);
    let cone = w.add(format!("CONICAL_SURFACE('',#{a},5.,0.4636476090008061)"));
    let vl = w.add(format!("VERTEX_LOOP('',#{apex})"));
    let point = w.add(format!("FACE_BOUND('',#{vl},.T.)"));
    let o = w.oe(e, true);
    let outer = w.bound(&[o], true, true);
    // The vertex loop written first: the outer bound is still the marked one.
    let side = w.face(&[point, outer], cone, true);
    let o = w.oe(e, false);
    let base_bound = w.bound(&[o], true, true);
    let plane = w.plane([0.0; 3], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
    let base = w.face(&[base_bound], plane, false);
    let shell = closed_shell(&mut w, &[side, base]);
    let solid = w.add(format!("MANIFOLD_SOLID_BREP('cone',#{shell})"));
    let s = one_solid(&w.file(&[solid]));
    assert_eq!(counts(&s), [1, 2, 1, 2]);
    assert!(s.shells[0].closed);
    let f = &s.faces[0];
    assert_eq!(f.bounds.len(), 2);
    assert_eq!(f.loop_order(), [1, 0]);
    // The vertex loop is read first: the apex is vertex 0.
    assert_eq!(
        f.bounds[0].r#loop,
        StepLoop::Vertex {
            entity: vl,
            vertex: 0
        }
    );
    assert!(!f.bounds[0].outer && f.bounds[0].face_loop().is_empty());
    assert_eq!(
        s.vertices[0],
        StepVertex {
            entity: apex,
            point: [0.0, 0.0, 10.0]
        }
    );
    assert_eq!(
        *f.surface,
        Surface::Cone {
            position: Placement {
                origin: [0.0; 3],
                axis: Some([0.0, 0.0, -1.0]),
                ref_direction: Some([1.0, 0.0, 0.0]),
            },
            radius: 5.0,
            semi_angle: 0.4636476090008061,
        }
    );
    // One closed edge: both ends are its one vertex, its curve the circle.
    let edge = &s.edges[0];
    assert_eq!((edge.start, edge.end, edge.curve_entity), (1, 1, k));
    assert!(matches!(&*edge.curve, Curve::Circle { radius, .. } if *radius == 5.0));
    assert!(!s.faces[1].same_sense);
    assert!(!s.faces[0].bounds[1].r#loop.seam_only());
}

/// A ball bounded by its meridian alone, used once each way.
#[test]
fn a_loop_of_seam_uses_alone() {
    let mut w = W::default();
    let (vn, vs) = (w.vertex([0.0, 0.0, 5.0]), w.vertex([0.0, 0.0, -5.0]));
    let k = w.circle([0.0; 3], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], 5.0);
    let e = w.add(format!("EDGE_CURVE('',#{vn},#{vs},#{k},.T.)"));
    let (o1, o2) = (w.oe(e, true), w.oe(e, false));
    let bound = w.bound(&[o1, o2], true, true);
    let a = w.axis([0.0; 3], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
    let sphere = w.add(format!("SPHERICAL_SURFACE('',#{a},5.)"));
    let face = w.face(&[bound], sphere, true);
    let shell = closed_shell(&mut w, &[face]);
    let solid = w.add(format!("MANIFOLD_SOLID_BREP('ball',#{shell})"));
    let s = one_solid(&w.file(&[solid]));
    assert_eq!(counts(&s), [1, 1, 1, 2]);
    assert!(s.shells[0].closed);
    assert!(s.faces[0].bounds[0].r#loop.seam_only());
    assert_eq!(
        s.edge_uses().unwrap(),
        [vec![
            StepEdgeUse {
                face: 0,
                bound: 0,
                position: 0,
                orientation: true
            },
            StepEdgeUse {
                face: 0,
                bound: 0,
                position: 1,
                orientation: false
            },
        ]]
    );
}

/// The cube with its first face listed through `wraps` `ORIENTED_FACE`s, the outer one `.F.`
fn cube_with_oriented_face(wraps: usize) -> (StepSolids, u64, u64) {
    let mut w = W::default();
    let mut written = Written::default();
    let mut faces = box_faces(
        &mut w,
        [0.0; 3],
        [1.0; 3],
        BoxStyle::default(),
        &mut written,
    );
    let base = faces[0];
    for k in 0..wraps {
        let o = if k + 1 == wraps { ".F." } else { ".T." };
        faces[0] = w.add(format!("ORIENTED_FACE('',*,#{},{o})", faces[0]));
    }
    let top = faces[0];
    let shell = closed_shell(&mut w, &faces);
    let solid = w.add(format!("MANIFOLD_SOLID_BREP('',#{shell})"));
    (read_solids(&w.file(&[solid])), base, top)
}

#[test]
fn oriented_faces_fold_into_the_shells_use() {
    for wraps in 1..=MAX_NESTING {
        let (r, base, top) = cube_with_oriented_face(wraps);
        assert!(r.failures.is_empty(), "{wraps}: {:?}", r.failures);
        let s = &r.solids[0];
        let of = s.shells[0].faces[0];
        assert_eq!((of.entity, of.face, of.orientation), (top, 0, false));
        assert_eq!(s.faces[0].entity, base);
    }
    // Two .F. wrappers turn the face back.
    let mut w = W::default();
    let mut written = Written::default();
    let mut faces = box_faces(
        &mut w,
        [0.0; 3],
        [1.0; 3],
        BoxStyle::default(),
        &mut written,
    );
    let inner = w.add(format!("ORIENTED_FACE('',*,#{},.F.)", faces[0]));
    faces[0] = w.add(format!("ORIENTED_FACE('',*,#{inner},.F.)"));
    let shell = closed_shell(&mut w, &faces);
    let solid = w.add(format!("MANIFOLD_SOLID_BREP('',#{shell})"));
    let s = one_solid(&w.file(&[solid]));
    assert!(s.shells[0].faces[0].orientation);
    // One more than the bound is refused.
    let (r, _, _) = cube_with_oriented_face(MAX_NESTING + 1);
    assert!(r.solids.is_empty());
    assert_eq!(r.failures[0].message, "oriented faces nested too deeply");
}

/// The cube with its first face's first oriented edge wrapped in `wraps` more `ORIENTED_EDGE`s,
/// each .F. (an odd number turn it over, so the edge is written the other way to start with).
fn cube_with_oriented_edge(wraps: usize) -> StepSolids {
    let mut w = W::default();
    let mut written = Written::default();
    box_faces(
        &mut w,
        [0.0; 3],
        [1.0; 3],
        BoxStyle::default(),
        &mut written,
    );
    // The first ORIENTED_EDGE written, rewritten under the wrappers.
    let text = w.s.clone();
    let first = text
        .lines()
        .find(|l| l.contains("=ORIENTED_EDGE("))
        .unwrap()
        .to_string();
    let (id, body) = first.split_once('=').unwrap();
    let id: u64 = id.trim_start_matches('#').parse().unwrap();
    let flipped = match (wraps % 2 == 1, body.ends_with(",.T.);")) {
        (false, _) => body.to_string(),
        (true, true) => body.replace(",.T.);", ",.F.);"),
        (true, false) => body.replace(",.F.);", ",.T.);"),
    };
    let mut inner = w.add(flipped.trim_end_matches(';'));
    for _ in 0..wraps {
        inner = w.add(format!("ORIENTED_EDGE('',*,*,#{inner},.F.)"));
    }
    w.s = w.s.replacen(
        &format!("EDGE_LOOP('',(#{id},"),
        &format!("EDGE_LOOP('',(#{inner},"),
        1,
    );
    let faces: Vec<u64> =
        w.s.lines()
            .filter(|l| l.contains("=ADVANCED_FACE("))
            .map(|l| {
                l.split_once('=')
                    .unwrap()
                    .0
                    .trim_start_matches('#')
                    .parse()
                    .unwrap()
            })
            .collect();
    let shell = closed_shell(&mut w, &faces);
    let solid = w.add(format!("MANIFOLD_SOLID_BREP('',#{shell})"));
    read_solids(&w.file(&[solid]))
}

#[test]
fn oriented_edges_in_oriented_edges_fold_into_one_orientation() {
    let plain = cube_with_oriented_edge(0);
    assert!(plain.failures.is_empty(), "{:?}", plain.failures);
    let want = plain.solids[0].faces[0].bounds[0].r#loop.clone();
    // ORIENTED_EDGE in ORIENTED_EDGE, up to the bound: the same edges, the same way round.
    for wraps in 1..MAX_NESTING {
        let r = cube_with_oriented_edge(wraps);
        assert!(r.failures.is_empty(), "{wraps}: {:?}", r.failures);
        let s = &r.solids[0];
        let (StepLoop::Edges { edges: got, .. }, StepLoop::Edges { edges: was, .. }) =
            (&s.faces[0].bounds[0].r#loop, &want)
        else {
            panic!()
        };
        let strip = |l: &[StepOrientedEdge]| {
            l.iter()
                .map(|o| (o.edge, o.orientation))
                .collect::<Vec<_>>()
        };
        assert_eq!(strip(got), strip(was), "{wraps}");
        assert_eq!(counts(s), [1, 6, 12, 8]);
    }
    let r = cube_with_oriented_edge(MAX_NESTING);
    assert_eq!(r.failures.len(), 1);
    assert!(
        r.failures[0]
            .message
            .ends_with("oriented edges nested too deeply"),
        "{:?}",
        r.failures
    );
}

#[test]
fn a_surface_model_of_open_shells() {
    let mut w = W::default();
    let mut written = Written::default();
    let mut faces = box_faces(
        &mut w,
        [0.0; 3],
        [1.0, 1.0, 1.0],
        BoxStyle::default(),
        &mut written,
    );
    let popped = faces.pop().unwrap();
    let open = w.add(format!("OPEN_SHELL('',({}))", list(&faces)));
    let single = box_faces(
        &mut w,
        [5.0; 3],
        [6.0; 3],
        BoxStyle::default(),
        &mut written,
    );
    let one = w.add(format!("OPEN_SHELL('',(#{}))", single[0]));
    let model = w.add(format!(
        "SHELL_BASED_SURFACE_MODEL('sheet',(#{open},#{one}))"
    ));
    // The face left out and the single box's other five faces are listed by no shell: loose
    // faces.
    let r = read_solids(&w.file(&[model]));
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    assert_eq!(r.solids.len(), 7);
    let s = &r.solids[0];
    assert_eq!(
        (s.kind, s.item, s.name.as_str()),
        (SolidKind::SurfaceModel, model, "sheet")
    );
    assert_eq!(counts(s), [2, 6, 12 + 4, 8 + 4]);
    assert!(s.shells.iter().all(|sh| !sh.declared_closed && !sh.closed));
    let once = s
        .edge_uses()
        .unwrap()
        .iter()
        .filter(|u| u.len() == 1)
        .count();
    assert_eq!(once, 4 + 4);
    let loose: Vec<u64> = r.solids[1..].iter().map(|s| s.item).collect();
    assert_eq!(loose, [&[popped][..], &single[1..]].concat());
    assert!(r.solids[1..].iter().all(|s| s.kind == SolidKind::LooseFace));

    // A model with a closed shell beside an open one: the closed shell first, a solid of its
    // own, then the model's open shell.
    let mut w = W::default();
    let faces = box_faces(
        &mut w,
        [0.0; 3],
        [1.0; 3],
        BoxStyle::default(),
        &mut written,
    );
    let shell = closed_shell(&mut w, &faces);
    let single = box_faces(
        &mut w,
        [5.0; 3],
        [6.0; 3],
        BoxStyle::default(),
        &mut written,
    );
    let one = w.add(format!("OPEN_SHELL('',({}))", list(&single)));
    let model = w.add(format!(
        "SHELL_BASED_SURFACE_MODEL('both',(#{one},#{shell}))"
    ));
    let ex = w.file(&[model]);
    assert_eq!(items(&ex, model), [shell, model]);
    let r = read_solids(&ex);
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let got: Vec<(SolidKind, u64, [usize; 4], bool)> = r
        .solids
        .iter()
        .map(|s| (s.kind, s.item, counts(s), s.shells[0].closed))
        .collect();
    assert_eq!(
        got,
        [
            (SolidKind::ModelClosedShell, shell, [1, 6, 12, 8], true),
            (SolidKind::SurfaceModel, model, [1, 6, 12, 8], true),
        ]
    );
    // An open shell the file calls open though every edge is used twice.
    assert!(!r.solids[1].shells[0].declared_closed);
}

#[test]
fn a_face_no_shell_lists_is_a_sheet_of_its_own() {
    let mut w = W::default();
    let mut written = Written::default();
    let faces = box_faces(
        &mut w,
        [0.0; 3],
        [1.0; 3],
        BoxStyle::default(),
        &mut written,
    );
    let shell = closed_shell(&mut w, &faces[1..]);
    let solid = w.add(format!("MANIFOLD_SOLID_BREP('',#{shell})"));
    // A loose ORIENTED_FACE over a loose face: the wrapper is the sheet, not its face.
    let other = box_faces(
        &mut w,
        [5.0; 3],
        [6.0; 3],
        BoxStyle::default(),
        &mut written,
    );
    let wrapped = w.add(format!("ORIENTED_FACE('',*,#{},.F.)", other[0]));
    // A loose ORIENTED_FACE over a face a shell lists: not a sheet.
    w.add(format!("ORIENTED_FACE('',*,#{},.F.)", faces[1]));
    closed_shell(&mut w, &other[1..]);
    let ex = w.file(&[solid]);
    assert_eq!(loose_faces(&ex), [faces[0], wrapped]);
    let r = read_solids(&ex);
    assert!(r.failures.is_empty(), "{:?}", r.failures);
    let got: Vec<(SolidKind, u64)> = r.solids.iter().map(|s| (s.kind, s.item)).collect();
    assert_eq!(
        got,
        [
            (SolidKind::ManifoldSolid, solid),
            (SolidKind::LooseFace, faces[0]),
            (SolidKind::LooseFace, wrapped),
        ]
    );
    // The solid misses a face: read, not closed.
    assert!(!r.solids[0].shells[0].closed && r.solids[0].shells[0].declared_closed);
    let sheet = &r.solids[1];
    assert_eq!(counts(sheet), [1, 1, 4, 4]);
    assert_eq!(
        sheet.shells[0],
        StepShell {
            entity: None,
            declared_closed: false,
            orientation: true,
            faces: vec![StepOrientedFace {
                entity: faces[0],
                face: 0,
                orientation: true
            }],
            closed: false,
        }
    );
    // A loose face in no representation is read in millimetres, with no uncertainty.
    assert_eq!((sheet.units, sheet.tolerance), (Units::default(), None));
    let turned = &r.solids[2];
    assert_eq!(
        turned.shells[0].faces[0],
        StepOrientedFace {
            entity: wrapped,
            face: 0,
            orientation: false
        }
    );
    assert_eq!(turned.faces[0].entity, other[0]);
}

#[test]
fn a_solid_without_a_name_takes_its_products() {
    let mut w = W::default();
    let mut written = Written::default();
    let faces = box_faces(
        &mut w,
        [0.0; 3],
        [1.0; 3],
        BoxStyle::default(),
        &mut written,
    );
    let shell = closed_shell(&mut w, &faces);
    let solid = w.add(format!("MANIFOLD_SOLID_BREP('',#{shell})"));
    let product = w.add("PRODUCT('p1','bracket','',())");
    let pdf = w.add(format!("PRODUCT_DEFINITION_FORMATION('','',#{product})"));
    let pd = w.add(format!("PRODUCT_DEFINITION('design','',#{pdf},$)"));
    let pds = w.add(format!("PRODUCT_DEFINITION_SHAPE('','',#{pd})"));
    // The representation `file` writes after the relationship and its five entities.
    let rep = w.n + 7;
    w.add(format!("SHAPE_DEFINITION_REPRESENTATION(#{pds},#{rep})"));
    let ex = w.file(&[solid]);
    assert_eq!(
        ex.get(rep).unwrap().name(),
        "ADVANCED_BREP_SHAPE_REPRESENTATION"
    );
    let s = one_solid(&ex);
    assert_eq!(
        (s.name.as_str(), s.product.as_deref()),
        ("bracket", Some("bracket"))
    );
}

/// The cube with one solid written twice over the same shell.
fn two_solids_on_one_shell() -> (Exchange, [u64; 2]) {
    let mut w = W::default();
    let mut written = Written::default();
    let faces = box_faces(
        &mut w,
        [0.0; 3],
        [1.0; 3],
        BoxStyle::default(),
        &mut written,
    );
    let shell = closed_shell(&mut w, &faces);
    let a = w.add(format!("MANIFOLD_SOLID_BREP('a',#{shell})"));
    let b = w.add(format!("MANIFOLD_SOLID_BREP('b',#{shell})"));
    (w.file(&[a, b]), [a, b])
}

#[test]
fn a_shared_shell_counts_for_each_solid_and_solids_are_bounded() {
    let (ex, [a, b]) = two_solids_on_one_shell();
    assert_eq!(read_solids(&ex).solids.len(), 2);
    // Eleven faces' work: the first solid's six, then the second runs out at its sixth.
    let r = read_solids_in(&ex, 11, MAX_SOLIDS);
    assert_eq!(r.solids.iter().map(|s| s.item).collect::<Vec<_>>(), [a]);
    assert_eq!(
        r.failures,
        [StepFailure {
            item: b,
            message: format!("more than {MAX_FILE_FACES} faces in one file")
        }]
    );
    assert_eq!(read_solids_in(&ex, 12, MAX_SOLIDS).solids.len(), 2);
    let r = read_solids_in(&ex, MAX_FILE_FACES, 1);
    assert_eq!(r.solids.len(), 1);
    assert_eq!(
        r.failures,
        [StepFailure {
            item: b,
            message: "more than 1 solids in one file".into()
        }]
    );
}

/// The cube, edited, read. Its first face is #48, that face's loop #41, its first edge #21 and
/// its first vertex #2.
fn edited(edit: impl Fn(&str) -> String) -> StepSolids {
    let mut w = W::default();
    let mut written = Written::default();
    let faces = box_faces(
        &mut w,
        [0.0; 3],
        [1.0; 3],
        BoxStyle::default(),
        &mut written,
    );
    let shell = closed_shell(&mut w, &faces);
    w.add(format!("MANIFOLD_SOLID_BREP('',#{shell})"));
    let text = edit(&w.s);
    read_solids(&parse(&text))
}

/// The id of the first entity whose line contains `what`.
fn first_id(t: &str, what: &str) -> u64 {
    let line = t.lines().find(|l| l.contains(what)).unwrap();
    line.split_once('=')
        .unwrap()
        .0
        .trim_start_matches('#')
        .parse()
        .unwrap()
}

/// The first reference after `after` in `t`.
fn first_ref(t: &str, after: &str) -> String {
    let (_, rest) = t.split_once(after).unwrap();
    rest.chars().take_while(char::is_ascii_digit).collect()
}

#[test]
fn hostile_topology_is_an_error_never_a_panic() {
    type Edit = Box<dyn Fn(&str) -> String>;
    assert_eq!(edited(str::to_string).solids.len(), 1);
    let cases: Vec<(String, Edit)> = vec![
        (
            "missing entity #9999149".into(),
            Box::new(|t: &str| {
                t.replace("MANIFOLD_SOLID_BREP('',#", "MANIFOLD_SOLID_BREP('',#9999")
            }),
        ),
        (
            "missing entity #9999".into(),
            Box::new(|t: &str| t.replacen("CLOSED_SHELL('',(#", "CLOSED_SHELL('',(#9999,#", 1)),
        ),
        (
            "face #48 is used twice".into(),
            Box::new(|t: &str| {
                let f = first_ref(t, "CLOSED_SHELL('',(#");
                t.replacen("CLOSED_SHELL('',(#", &format!("CLOSED_SHELL('',(#{f},#"), 1)
            }),
        ),
        (
            "oriented faces nested too deeply".into(),
            Box::new(|t: &str| {
                format!("{t}#9000=ORIENTED_FACE('',*,#9000,.T.);").replacen(
                    "CLOSED_SHELL('',(#",
                    "CLOSED_SHELL('',(#9000,#",
                    1,
                )
            }),
        ),
        (
            "face #48: oriented edges nested too deeply".into(),
            Box::new(|t: &str| {
                format!("{t}#9000=ORIENTED_EDGE('',*,*,#9000,.T.);").replacen(
                    "EDGE_LOOP('',(#",
                    "EDGE_LOOP('',(#9000,#",
                    1,
                )
            }),
        ),
        (
            "face #48: edge loop #41 is not closed".into(),
            Box::new(|t: &str| {
                let (head, rest) = t.split_once("EDGE_LOOP('',(#").unwrap();
                let (_, rest) = rest.split_once(',').unwrap();
                format!("{head}EDGE_LOOP('',({rest}")
            }),
        ),
        (
            "face #48: edge loop #41 is not connected".into(),
            Box::new(|t: &str| {
                // The loop's second and third oriented edges swapped.
                let (head, rest) = t.split_once("EDGE_LOOP('',(").unwrap();
                let (list, rest) = rest.split_once(')').unwrap();
                let mut l: Vec<&str> = list.split(',').collect();
                l.swap(1, 2);
                format!("{head}EDGE_LOOP('',({}){rest}", l.join(","))
            }),
        ),
        (
            "face #48: #2: not supported yet: vertex VERTEX_X".into(),
            Box::new(|t: &str| t.replacen("=VERTEX_POINT(", "=VERTEX_X(", 1)),
        ),
        (
            "face #48: #21: not supported yet: edge SUBEDGE".into(),
            Box::new(|t: &str| t.replacen("=EDGE_CURVE(", "=SUBEDGE(", 1)),
        ),
        (
            "face #48: vertex #2: #18: expected CARTESIAN_POINT, found DIRECTION".into(),
            Box::new(|t: &str| {
                let (head, rest) = t.split_once("VERTEX_POINT('',#").unwrap();
                let (_, rest) = rest.split_once(')').unwrap();
                let d = first_id(t, "=DIRECTION");
                format!("{head}VERTEX_POINT('',#{d}){rest}")
            }),
        ),
        (
            "face #48: not supported yet: surface HYPERBOLA".into(),
            Box::new(|t: &str| t.replacen("PLANE('',", "HYPERBOLA('',", 1)),
        ),
        (
            format!("face #48: edge loop #41 with {} edges", MAX_LOOP + 5),
            Box::new(|t: &str| {
                let first = first_ref(t, "EDGE_LOOP('',(#");
                let many = format!("#{first},").repeat(MAX_LOOP + 1);
                t.replacen("EDGE_LOOP('',(#", &format!("EDGE_LOOP('',({many}#"), 1)
            }),
        ),
        (
            "face same_sense: expected a boolean".into(),
            Box::new(|t: &str| {
                let (head, rest) = t.split_once("ADVANCED_FACE(").unwrap();
                format!(
                    "{head}ADVANCED_FACE({}",
                    rest.replacen(",.T.);", ",.U.);", 1)
                )
            }),
        ),
        (
            "face #48: #41: not supported yet: loop POLY_LOOP".into(),
            Box::new(|t: &str| t.replacen("EDGE_LOOP(", "POLY_LOOP(", 1)),
        ),
        (
            "not supported yet: FACETED_BREP".into(),
            Box::new(|t: &str| t.replace("MANIFOLD_SOLID_BREP", "FACETED_BREP")),
        ),
        (
            "#48: not supported yet: face FACE".into(),
            Box::new(|t: &str| t.replacen("=ADVANCED_FACE(", "=FACE(", 1)),
        ),
        (
            "face #48 has no bounds".into(),
            Box::new(|t: &str| {
                let (head, rest) = t.split_once("=ADVANCED_FACE('',(").unwrap();
                let (_, rest) = rest.split_once(')').unwrap();
                format!("{head}=ADVANCED_FACE('',(){rest}")
            }),
        ),
        (
            "face #48: edge #21: edge geometry nested too deeply".into(),
            Box::new(|t: &str| {
                // The first edge's curve a surface curve that is its own 3D curve.
                let (head, rest) = t.split_once("EDGE_CURVE(").unwrap();
                let mut parts: Vec<String> = rest.splitn(4, ',').map(str::to_string).collect();
                let (_, tail) = parts[3].split_once(',').unwrap();
                parts[3] = format!("#9000,{tail}");
                format!(
                    "{head}EDGE_CURVE({}#9000=SURFACE_CURVE('',#9000,(),.CURVE_3D.);",
                    parts.join(",")
                )
            }),
        ),
    ];
    for (what, edit) in cases {
        let r = edited(edit);
        assert!(
            r.solids.is_empty() && r.failures.len() == 1,
            "{what}: {:?} {}",
            r.failures,
            r.solids.len()
        );
        assert_eq!(r.failures[0].message, what);
    }
}

#[test]
fn a_shell_of_too_many_faces_is_refused_before_it_is_read() {
    let r = edited(|t: &str| {
        let first = first_ref(t, "CLOSED_SHELL('',(#");
        let many = format!("#{first},").repeat(MAX_FACES);
        t.replacen(
            "CLOSED_SHELL('',(#",
            &format!("CLOSED_SHELL('',({many}#"),
            1,
        )
    });
    assert!(r.solids.is_empty());
    // The shell is written just before the solid.
    assert_eq!(
        r.failures[0].message,
        format!(
            "shell #{} with {} faces",
            r.failures[0].item - 1,
            MAX_FACES + 6
        )
    );
}

#[test]
fn an_edge_on_a_surface_curve_or_trimmed_curve_is_its_3d_curve() {
    let mut w = W::default();
    let mut written = Written::default();
    let faces = box_faces(
        &mut w,
        [0.0; 3],
        [1.0; 3],
        BoxStyle::default(),
        &mut written,
    );
    // The first edge's line under a trimmed curve under a surface curve whose pcurve's
    // representation has no items: the pcurve is not read.
    let text = w.s.clone();
    let (id, body) = text
        .lines()
        .find(|l| l.contains("=EDGE_CURVE("))
        .unwrap()
        .split_once('=')
        .unwrap();
    let edge: u64 = id.trim_start_matches('#').parse().unwrap();
    let line: u64 = body
        .split(',')
        .nth(3)
        .unwrap()
        .trim_start_matches('#')
        .parse()
        .unwrap();
    let ctx = w.add("(GEOMETRIC_REPRESENTATION_CONTEXT(2)PARAMETRIC_REPRESENTATION_CONTEXT()REPRESENTATION_CONTEXT('',''))");
    let rep = w.add(format!("DEFINITIONAL_REPRESENTATION('',(),#{ctx})"));
    let plane = w.plane([0.0; 3], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
    let pc = w.add(format!("PCURVE('',#{plane},#{rep})"));
    let trimmed = w.add(format!(
        "TRIMMED_CURVE('',#{line},(PARAMETER_VALUE(0.)),(PARAMETER_VALUE(1.)),.T.,.PARAMETER.)"
    ));
    let sc = w.add(format!("SURFACE_CURVE('',#{trimmed},(#{pc}),.CURVE_3D.)"));
    let new_body = body.replacen(&format!("#{line},"), &format!("#{sc},"), 1);
    w.s = w.s.replacen(
        &format!("#{edge}={body}"),
        &format!("#{edge}={new_body}"),
        1,
    );
    let shell = closed_shell(&mut w, &faces);
    let solid = w.add(format!("MANIFOLD_SOLID_BREP('',#{shell})"));
    let s = one_solid(&w.file(&[solid]));
    let e = s.edges.iter().find(|e| e.entity == edge).unwrap();
    assert_eq!(e.curve_entity, line);
    assert!(
        matches!(&*e.curve, Curve::Line { vector: Vector { magnitude, .. }, .. } if *magnitude == 1.0)
    );
}

/// A whole sphere bounded by one vertex loop, as `jack_in_the_box.stp` writes its balls: data,
/// with no edge.
#[test]
fn a_face_bounded_by_a_vertex_loop_alone() {
    let mut w = W::default();
    let v = w.vertex([0.0, 0.0, 5.0]);
    let vl = w.add(format!("VERTEX_LOOP('',#{v})"));
    let bound = w.add(format!("FACE_BOUND('',#{vl},.T.)"));
    let a = w.axis([0.0; 3], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
    let sphere = w.add(format!("SPHERICAL_SURFACE('',#{a},5.)"));
    let face = w.face(&[bound], sphere, true);
    let shell = closed_shell(&mut w, &[face]);
    let solid = w.add(format!("MANIFOLD_SOLID_BREP('ball',#{shell})"));
    let s = one_solid(&w.file(&[solid]));
    assert_eq!(counts(&s), [1, 1, 0, 1]);
    let f = &s.faces[0];
    assert_eq!(f.bounds.len(), 1);
    assert_eq!(
        f.bounds[0],
        StepBound {
            entity: bound,
            outer: false,
            orientation: true,
            r#loop: StepLoop::Vertex {
                entity: vl,
                vertex: 0
            },
        }
    );
    assert_eq!(f.loop_order(), [0]);
    // No edge is used fewer than twice.
    assert!(s.shells[0].closed);
    assert_eq!(s.edge_uses().unwrap(), Vec::<Vec<StepEdgeUse>>::new());
}
