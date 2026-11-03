//! Display tessellation plans: the parameter-space triangulation of every face of a body,
//! with the surfaces to evaluate it on, kept apart from the evaluation itself.
//!
//! The planner ([`crate::plan`]) decides *where* the vertices go (in each face's (u, v) parameters, with the
//! shared edge nodes fixed in 3D so neighbouring faces meet exactly) and *how* they are
//! connected. Evaluating the surfaces at those parameters is the bulk of the work and has no
//! branching on topology, so it can run anywhere: [`TessPlan::evaluate`] does it on the CPU,
//! and the GPU evaluator runs the same maths in a compute shader. Both evaluate the kernel's
//! [`EvalSurface`] forms, so the two cannot drift apart.

use crate::{Mesh, RefusedFace};
use solidcore::eval::EvalSurface;
use solidcore_convert::{kvec3, vec3};
use solidcore_record::Vec3;

/// One face's share of a [`TessPlan`].
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FacePlan {
    /// The B-rep face (index in the body's face order).
    pub face: u32,
    /// Index into [`TessPlan::surfaces`]; `None` when every position and normal is given.
    pub surface: Option<u32>,
    /// Parameters of every vertex.
    pub uv: Vec<[f64; 2]>,
    /// Positions of the first `fixed.len()` vertices: the shared edge nodes, used exactly as
    /// they are so the neighbouring face gets bit-identical boundary points.
    pub fixed: Vec<Vec3>,
    /// Normals of every vertex, only for faces without an evaluable surface.
    pub normals: Vec<Vec3>,
    /// The face is reversed against its surface: normals point the other way.
    pub flip: bool,
    /// Counter-clockwise seen from outside (orientation already applied).
    pub triangles: Vec<[u32; 3]>,
}

impl FacePlan {
    /// Triangle `t` is on this face's own nodes (the evaluated mesh keeps only those).
    fn has_nodes(&self, t: &[u32; 3]) -> bool {
        t.iter().all(|&i| (i as usize) < self.uv.len())
    }
}

/// A body's display tessellation, ready to evaluate.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TessPlan {
    /// The surfaces the faces are evaluated on, in the kernel evaluator's form.
    pub surfaces: Vec<EvalSurface>,
    /// Each face's share, in the body's face order.
    pub faces: Vec<FacePlan>,
    /// B-rep edges as polylines (the shared edge nodes), in the body's edge order.
    pub edges: Vec<Vec<Vec3>>,
    /// Per edge: the faces it bounds.
    pub edge_faces: Vec<Vec<u32>>,
}

impl TessPlan {
    /// Evaluate on the CPU into an indexed mesh (vertices are shared within a face, not
    /// across faces, so creases keep their own normals), with the faces left out of it.
    ///
    /// A face a vertex of which the kernel's evaluator refuses (no value there, or a pole of
    /// higher order with no limit normal) is left out whole, nothing substituted, and returned
    /// with why: the caller reports it as unmeshed.
    #[must_use]
    pub fn evaluate(&self) -> (Mesh, Vec<RefusedFace>) {
        let mut out = Mesh::default();
        let mut refused = Vec::new();
        for f in &self.faces {
            let surface = f.surface.and_then(|i| self.surfaces.get(i as usize));
            let mut positions = Vec::with_capacity(f.uv.len());
            let mut normals = Vec::with_capacity(f.uv.len());
            let mut why = None;
            for (k, uv) in f.uv.iter().enumerate() {
                let (p, n) = match surface {
                    Some(s) => match s.eval(uv[0], uv[1]) {
                        Ok(pn) => solidcore_convert::record_pair(pn),
                        Err(e) => {
                            why = Some(format!("vertex {k} at {uv:?}: {}", crate::eval_message(e)));
                            break;
                        }
                    },
                    None => (Vec3::ZERO, f.normals.get(k).copied().unwrap_or(Vec3::Z)),
                };
                positions.push(f.fixed.get(k).copied().unwrap_or(p));
                normals.push(if f.flip { -n } else { n });
            }
            if let Some(reason) = why {
                refused.push(RefusedFace { face: f.face as usize, reason });
                continue;
            }
            let base = out.positions.len();
            out.positions.extend(positions);
            out.normals.extend(normals);
            for t in &f.triangles {
                if f.has_nodes(t) {
                    let at = |i: u32| u32::try_from(base + i as usize).unwrap_or(u32::MAX);
                    out.triangles.push([at(t[0]), at(t[1]), at(t[2])]);
                    out.tri_face.push(f.face);
                }
            }
        }
        out.edges = self.edges.clone();
        out.edge_faces = self.edge_faces.clone();
        (out, refused)
    }

    /// Which of the body's `faces` faces have triangles in [`TessPlan::evaluate`]'s mesh,
    /// without evaluating it.
    pub fn faces_with_triangles(&self, faces: usize) -> Vec<bool> {
        let mut has = vec![false; faces];
        for f in &self.faces {
            if f.triangles.iter().any(|t| f.has_nodes(t))
                && let Some(h) = has.get_mut(f.face as usize)
            {
                *h = true;
            }
        }
        has
    }

    /// Vertices across all faces.
    pub fn vertex_count(&self) -> usize {
        self.faces.iter().map(|f| f.uv.len()).sum()
    }

    /// This plan for the body moved by the rotation and translation `m` (by columns, the
    /// translation in `m[3]`): the same parameters and triangles on the moved surfaces, so a
    /// part placed many times is planned once. Only for proper rigid motions (a reflection
    /// would turn the triangles inside out; the caller checks).
    pub fn transformed(&self, m: &[[f64; 4]; 4]) -> TessPlan {
        let vector = |v: Vec3| {
            Vec3::new(
                m[0][0] * v.x + m[1][0] * v.y + m[2][0] * v.z,
                m[0][1] * v.x + m[1][1] * v.y + m[2][1] * v.z,
                m[0][2] * v.x + m[1][2] * v.y + m[2][2] * v.z,
            )
        };
        let point = |p: Vec3| vector(p) + Vec3::new(m[3][0], m[3][1], m[3][2]);
        let surfaces = self
            .surfaces
            .iter()
            .map(|s| match s {
                EvalSurface::Plane { origin, du, dv } => {
                    EvalSurface::Plane { origin: kvec3(point(vec3(*origin))), du: kvec3(vector(vec3(*du))), dv: kvec3(vector(vec3(*dv))) }
                }
                // An affine map keeps a surface closed or open.
                EvalSurface::Nurbs(n) => EvalSurface::nurbs_closed(
                    n.degree(),
                    n.knots().clone(),
                    n.count(),
                    // Homogeneous: the translation scales with the weight.
                    n.ctrl()
                        .iter()
                        .map(|c| {
                            let q = vector(Vec3::new(c[0], c[1], c[2])) + Vec3::new(m[3][0], m[3][1], m[3][2]) * c[3];
                            [q.x, q.y, q.z, c[3]]
                        })
                        .collect(),
                    n.closed(),
                ),
                EvalSurface::Revolved { profile, origin, axis, transform, swap } => {
                    // `m` (as rows) after the surface's own affine map.
                    let r = |i: usize| [m[0][i], m[1][i], m[2][i], m[3][i]];
                    let rows = [r(0), r(1), r(2)];
                    let t = transform;
                    let composed = rows.map(|row| {
                        let mut out = [0.0; 4];
                        for (j, slot) in out.iter_mut().enumerate() {
                            *slot = (0..3).map(|k| row[k] * t[k][j]).sum::<f64>() + if j == 3 { row[3] } else { 0.0 };
                        }
                        out
                    });
                    EvalSurface::Revolved { profile: profile.clone(), origin: *origin, axis: *axis, transform: composed, swap: *swap }
                }
            })
            .collect();
        let faces = self
            .faces
            .iter()
            .map(|f| FacePlan {
                fixed: f.fixed.iter().map(|p| point(*p)).collect(),
                normals: f.normals.iter().map(|n| vector(*n)).collect(),
                ..f.clone()
            })
            .collect();
        TessPlan {
            surfaces,
            faces,
            edges: self.edges.iter().map(|e| e.iter().map(|p| point(*p)).collect()).collect(),
            edge_faces: self.edge_faces.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vec::VecExt as _;
    use solidcore::eval::EvalCurve;

    fn close(a: Vec3, b: Vec3, tol: f64) -> bool {
        a.dist(b) <= tol
    }

    /// A moved plan evaluates to the plan's points moved, on a rational spline surface (the
    /// translation scales with each weight) and a revolved one (composed with its own map).
    #[test]
    fn a_transformed_plan_evaluates_to_the_moved_points() {
        let (c, s) = (0.6f64.cos(), 0.6f64.sin());
        // Turned about z by 0.6 rad, then moved: by columns.
        let m = [[c, s, 0.0, 0.0], [-s, c, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [5.0, -2.0, 7.0, 1.0]];
        let apply = |p: Vec3| Vec3::new(c * p.x - s * p.y + 5.0, s * p.x + c * p.y - 2.0, p.z + 7.0);
        let w = std::f64::consts::FRAC_1_SQRT_2;
        let nurbs = EvalSurface::nurbs(
            [2, 1],
            [vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], vec![0.0, 0.0, 1.0, 1.0]],
            [3, 2],
            vec![[1.0, 0.0, 0.0, 1.0], [w, w, 0.0, w], [0.0, 1.0, 0.0, 1.0], [1.0, 0.0, 2.0, 1.0], [w, w, 2.0 * w, w], [0.0, 1.0, 2.0, 1.0]],
        );
        let revolved = EvalSurface::Revolved {
            profile: line(Vec3::new(2.0, 0.0, 0.0), Vec3::new(3.0, 0.0, 1.0)),
            origin: solidcore::Vec3::ZERO,
            axis: solidcore::Vec3::Z,
            transform: [[1.0, 0.0, 0.0, 1.0], [0.0, 0.0, -1.0, 2.0], [0.0, 1.0, 0.0, 3.0]],
            swap: false,
        };
        let face = |i: u32| FacePlan {
            face: i,
            surface: Some(i),
            uv: vec![[0.1, 0.2], [0.7, 0.4], [0.3, 0.9]],
            fixed: vec![Vec3::new(1.0, 2.0, 3.0)],
            triangles: vec![[0, 1, 2]],
            ..Default::default()
        };
        let plan = TessPlan {
            surfaces: vec![nurbs, revolved],
            faces: vec![face(0), face(1)],
            edges: vec![vec![Vec3::new(1.0, 1.0, 1.0)]],
            edge_faces: vec![vec![0, 1]],
        };
        let (a, b) = (plan.evaluate().0, plan.transformed(&m).evaluate().0);
        assert_eq!(a.triangles, b.triangles);
        for (p, q) in a.positions.iter().zip(&b.positions) {
            assert!(close(apply(*p), *q, 1e-12), "{p:?} -> {q:?}");
        }
        for (n, k) in a.normals.iter().zip(&b.normals) {
            assert!(close(apply(*n) - apply(Vec3::ZERO), *k, 1e-12));
        }
        assert!(close(apply(a.edges[0][0]), b.edges[0][0], 1e-12));
    }

    fn line(a: Vec3, b: Vec3) -> EvalCurve {
        EvalCurve::new(1, vec![0.0, 0.0, 1.0, 1.0], vec![[a.x, a.y, a.z, 1.0], [b.x, b.y, b.z, 1.0]])
    }

    #[test]
    fn evaluate_uses_fixed_points_and_drops_bad_indices() {
        let plan = TessPlan {
            surfaces: vec![EvalSurface::Plane { origin: solidcore::Vec3::ZERO, du: solidcore::Vec3::X, dv: solidcore::Vec3::Y }],
            faces: vec![FacePlan {
                face: 3,
                surface: Some(0),
                uv: vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
                fixed: vec![Vec3::new(0.0, 0.0, 1e-9)],
                normals: vec![],
                flip: true,
                triangles: vec![[0, 1, 2], [0, 1, 7]],
            }],
            edges: vec![],
            edge_faces: vec![],
        };
        let m = plan.evaluate().0;
        assert_eq!(m.positions.len(), 3);
        assert_eq!(m.positions.first().copied(), Some(Vec3::new(0.0, 0.0, 1e-9)));
        assert_eq!(m.triangles, vec![[0, 1, 2]]);
        assert_eq!(m.tri_face, vec![3]);
        assert!(m.normals.iter().all(|n| close(*n, -Vec3::Z, 1e-12)));
        // Faces with triangles, counted without evaluating, are the faces the mesh has: not a
        // face whose only triangle is off its nodes, nor one with none, nor one out of range.
        let mut plan = plan;
        let only = |face: u32, triangles: Vec<[u32; 3]>| FacePlan { face, uv: vec![[0.0, 0.0]; 3], triangles, ..Default::default() };
        plan.faces.extend([only(0, vec![[0, 1, 9]]), only(1, vec![]), only(2, vec![[2, 1, 0]]), only(7, vec![[0, 1, 2]])]);
        let has = plan.faces_with_triangles(5);
        let mut evaluated = vec![false; 5];
        for f in plan.evaluate().0.tri_face {
            if let Some(h) = evaluated.get_mut(f as usize) {
                *h = true;
            }
        }
        assert_eq!(has, evaluated);
        assert_eq!(has, vec![false, false, true, true, false]);
    }
}
