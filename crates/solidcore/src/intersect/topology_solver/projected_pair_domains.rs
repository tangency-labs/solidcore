//! Original-domain accounting for native projected topology pair maxima.
//!
//! Not upstream. Complete tree object groups are recorded before scalar pruning;
//! a retained point's enclosure never stands for an unvisited domain's proof.
//! Deviations: at most 4,096 recorded operands and 16,384 domain pairs. Exceeding
//! either budget refuses an activated native pair query; no Cartesian allocation.

use super::p1::{self, LocalSolver, LocalSolverKind, TopologyObjects};
use super::{projected_affine_face as native, projected_bound, projected_residual};
use crate::blend::fillet_section::interval::Ivl;
use crate::curve::curve::{p1::affine_curve_endpoints, p9, Curve};
use crate::error::{Result, SolidError};
use crate::geom::extent1d::Extent1d;
use crate::geom::global_solver::{GlobalSolver, SolutionRequested, SolverOperation};
use crate::geom::solution_array::{ObjectKind, Solution, SolutionObject};
use crate::geom::tree::Tree;
use crate::math::{Vec3, BIG_DOUBLE};
use crate::nurbs::b_spline_curve::BSplineCurve;
use crate::Kernel;

#[derive(Clone, Copy)]
enum Domain<'a, O> {
    Point(Vec3, Option<O>),
    Curve {
        curve: &'a Curve,
        interval: Extent1d,
        owner: O,
    },
    Face(O),
    Unknown,
}
struct Group<'a, O> {
    first: Vec<Domain<'a, O>>,
    second: Vec<Domain<'a, O>>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum CaptureState {
    Fresh,
    Recorded,
}

/// Solver-local proof history; references retain the original native operands.
pub(super) struct Ledger<'a, O> {
    groups: Vec<Group<'a, O>>,
    operands: usize,
    uncovered: bool,
    capture: CaptureState,
    target: Option<f64>,
    normal: Option<Vec3>,
    pub active: bool,
    pub in_tree: bool,
    best: Option<(Solution<O>, Ivl)>,
}
impl<O> Default for Ledger<'_, O> {
    fn default() -> Self {
        Self {
            groups: Vec::new(),
            operands: 0,
            uncovered: false,
            capture: CaptureState::Fresh,
            target: None,
            normal: None,
            active: false,
            in_tree: false,
            best: None,
        }
    }
}
#[derive(Clone, Copy)]
enum Role {
    Point,
    Vertex,
    Curve,
    Curve2,
    Edge,
    Face,
    Unknown,
}
fn roles(kind: LocalSolverKind) -> [Role; 2] {
    use LocalSolverKind as K;
    match kind {
        K::PointVertex => [Role::Point, Role::Vertex],
        K::VertexVertex => [Role::Vertex, Role::Vertex],
        K::PointCurve => [Role::Point, Role::Curve],
        K::PointEdge => [Role::Point, Role::Edge],
        K::VertexEdge => [Role::Vertex, Role::Edge],
        K::PointFace => [Role::Point, Role::Face],
        K::VertexFace => [Role::Vertex, Role::Face],
        K::CurveVertex => [Role::Curve, Role::Vertex],
        K::EdgeVertex => [Role::Edge, Role::Vertex],
        K::CurveCurve => [Role::Curve, Role::Curve2],
        K::CurveEdge => [Role::Curve, Role::Edge],
        K::EdgeEdge => [Role::Edge, Role::Edge],
        K::CurveFace => [Role::Curve, Role::Face],
        K::EdgeFace => [Role::Edge, Role::Face],
        K::FaceFace => [Role::Face, Role::Face],
        _ => [Role::Unknown, Role::Unknown],
    }
}
fn domain<'a, O: SolutionObject>(
    objects: &'a dyn TopologyObjects<O>,
    local: &LocalSolver<'a, O>,
    role: Role,
    object: Option<O>,
) -> Result<Option<Domain<'a, O>>> {
    let missing = || native::refused("original pair domain has an unnamed operand");
    Ok(Some(match role {
        Role::Point => Domain::Point(local.point, local.point_vertex),
        Role::Vertex => {
            let owner = object.ok_or_else(missing)?;
            Domain::Point(objects.vertex_point(owner)?, Some(owner))
        }
        Role::Curve | Role::Curve2 => {
            let second = matches!(role, Role::Curve2);
            let supplied = if second { local.curve2 } else { local.curve };
            let face = if second {
                local.curve_face2
            } else {
                local.curve_face
            };
            if face.is_some() {
                return Ok(Some(Domain::Unknown));
            }
            if let Some(supplied) = supplied {
                Domain::Curve {
                    curve: supplied.curve,
                    interval: if second {
                        local.interval2
                    } else {
                        local.interval
                    },
                    owner: supplied.handle,
                }
            } else {
                let handle = object.ok_or_else(missing)?;
                let curve = objects.curve(handle)?;
                Domain::Curve {
                    curve,
                    interval: curve.natural_interval(),
                    owner: handle,
                }
            }
        }
        Role::Edge => {
            let handle = object.ok_or_else(missing)?;
            let edge = if handle.kind() == ObjectKind::Edge {
                handle
            } else {
                objects.curve_edge(handle)?.ok_or_else(missing)?
            };
            if !local.in_shape(objects, Some(edge)) {
                return Ok(None);
            }
            let curve = objects.curve(objects.edge_curve(edge)?.ok_or_else(missing)?)?;
            Domain::Curve {
                curve,
                interval: objects.edge_interval(edge)?,
                owner: edge,
            }
        }
        Role::Face => {
            let handle = object.ok_or_else(missing)?;
            let face = if handle.kind() == ObjectKind::Face {
                handle
            } else {
                objects.surface_face(handle)?.ok_or_else(missing)?
            };
            if !local.in_shape(objects, Some(face)) {
                return Ok(None);
            }
            Domain::Face(face)
        }
        Role::Unknown => Domain::Unknown,
    }))
}
impl<'a, O: SolutionObject> Ledger<'a, O> {
    pub(super) fn eligible(solver: &GlobalSolver<'_, O>) -> bool {
        solver.operation == SolverOperation::ProjectedMaximize
            && solver.solution_requested == SolutionRequested::Single
    }
    fn begin(&mut self, solver: &GlobalSolver<'_, O>) {
        if self.capture == CaptureState::Fresh {
            self.capture = CaptureState::Recorded;
            self.normal = solver.opt_vectors.first().copied();
            self.target = if solver.best_answer_so_far_sq == BIG_DOUBLE {
                None
            } else {
                p1::best_answer_so_far(solver)
            };
            if !solver.solutions.is_empty() {
                self.uncovered = true;
            }
        }
    }
    fn append(&mut self, first: Vec<Domain<'a, O>>, second: Vec<Domain<'a, O>>) -> Result<()> {
        if first.is_empty() || second.is_empty() {
            return Ok(());
        }
        for (curves, faces) in [(&first, &second), (&second, &first)] {
            if faces.iter().any(|d| matches!(d, Domain::Face(_))) {
                for d in curves {
                    if let Domain::Curve {
                        curve, interval, ..
                    } = d
                    {
                        if affine_curve_endpoints(curve, *interval)?.is_some() {
                            self.active = true;
                        }
                    }
                }
            }
        }
        self.operands = self
            .operands
            .saturating_add(first.len())
            .saturating_add(second.len());
        if self.operands > 4096 {
            self.uncovered = true;
        } else {
            self.groups.push(Group { first, second });
        }
        if self.active && self.uncovered {
            return Err(native::refused(
                "native pair ledger cannot cover the original query history",
            ));
        }
        Ok(())
    }
    pub fn record_local(
        &mut self,
        solver: &GlobalSolver<'_, O>,
        objects: &'a dyn TopologyObjects<O>,
        local: &LocalSolver<'a, O>,
        first: Option<O>,
        second: Option<O>,
    ) -> Result<()> {
        if !Self::eligible(solver) || self.in_tree {
            return Ok(());
        }
        self.begin(solver);
        let [a, b] = roles(local.kind);
        let (a, b) = if local.swap_order { (b, a) } else { (a, b) };
        let first = domain(objects, local, a, first)?.into_iter().collect();
        let second = domain(objects, local, b, second)?.into_iter().collect();
        self.append(first, second)
    }
    pub fn record_trees(
        &mut self,
        solver: &GlobalSolver<'_, O>,
        objects: &'a dyn TopologyObjects<O>,
        local: &LocalSolver<'a, O>,
        first: &Tree<O>,
        second: &Tree<O>,
    ) -> Result<()> {
        if !Self::eligible(solver) {
            return Ok(());
        }
        self.begin(solver);
        let [a, b] = roles(local.kind);
        let (a, b) = if local.swap_order { (b, a) } else { (a, b) };
        let mut groups = [Vec::new(), Vec::new()];
        for ((output, tree), role) in groups.iter_mut().zip([first, second]).zip([a, b]) {
            for id in tree.object_list()? {
                if output.len() >= 4096 {
                    self.uncovered = true;
                    break;
                }
                if let Some(d) = domain(objects, local, role, tree.list(id)?.object)? {
                    output.push(d);
                }
            }
        }
        let [first, second] = groups;
        self.append(first, second)
    }
    pub fn consider(&mut self, solution: &Solution<O>, enclosure: Ivl) {
        self.active = true;
        if self
            .best
            .as_ref()
            .is_none_or(|(_, old)| enclosure.hi > old.hi)
        {
            self.best = Some((*solution, enclosure));
        }
    }
    fn query_curve(&self, solution: &Solution<O>) -> Option<(O, &'a Curve)> {
        self.groups
            .iter()
            .flat_map(|g| g.first.iter().chain(&g.second))
            .find_map(|d| match d {
                Domain::Curve { curve, owner, .. } if solution.objects.contains(&Some(*owner)) => {
                    Some((*owner, *curve))
                }
                _ => None,
            })
    }
    fn contains_solution(
        &self,
        objects: &dyn TopologyObjects<O>,
        solution: &Solution<O>,
    ) -> Result<bool> {
        if solution.num_objects != 2 {
            return Ok(false);
        }
        let Some(a_owner) = solution.objects.first().copied().flatten() else {
            return Ok(false);
        };
        let Some(b_owner) = solution.objects.get(1).copied().flatten() else {
            return Ok(false);
        };
        let mut checks = 0usize;
        for group in &self.groups {
            for a in &group.first {
                for b in &group.second {
                    checks += 1;
                    if checks > 16384 {
                        return Err(native::refused(
                            "candidate membership exceeds the original domain-pair budget",
                        ));
                    }

                    let mut parameter = 0;
                    if contains_point(objects, *a, a_owner, solution, &mut parameter)?
                        && contains_point(objects, *b, b_owner, solution, &mut parameter)?
                        && parameter == solution.num_variables
                    {
                        return Ok(true);
                    }
                }
            }
        }
        Ok(false)
    }
    pub fn finish(
        &mut self,
        k: &Kernel,
        objects: &'a dyn TopologyObjects<O>,
        solver: &mut GlobalSolver<'_, O>,
    ) -> Result<()> {
        if !self.active {
            return Ok(());
        }
        if self.uncovered {
            return Err(native::refused(
                "native pair ledger has uncovered original domains",
            ));
        }
        let normal = *solver
            .opt_vectors
            .first()
            .ok_or_else(|| native::refused("missing ledger projection normal"))?;
        if !Self::eligible(solver) || self.normal != Some(normal) {
            return Err(native::refused(
                "native pair metric changed during its domain ledger",
            ));
        }
        if !solver.tolerance_3d.is_finite()
            || solver.tolerance_3d < 0.0
            || self.target.is_some_and(|target| !target.is_finite())
        {
            return Err(native::refused(
                "native pair has invalid retention or external target",
            ));
        }
        let mut chosen = match self.best {
            Some((solution, bound)) if self.contains_solution(objects, &solution)? => {
                Some((solution, bound))
            }
            _ => None,
        };
        let current = solver.solutions.get(0);
        let current_admissible = current
            .map(|current| self.contains_solution(objects, current))
            .transpose()?
            .unwrap_or(false);
        if let Some(current) = current.filter(|_| current_admissible) {
            if let Some((owner, curve)) = self.query_curve(current) {
                match native::retained_bound(k, objects, current, owner, curve, normal) {
                    Ok(bound) if chosen.as_ref().is_none_or(|(_, old)| bound.hi > old.hi) => {
                        chosen = Some((*current, bound));
                    }
                    Ok(_) | Err(SolidError::Unsupported { .. }) => {}
                    Err(error) => return Err(error),
                }
            }
        }
        let (solution, bound) = chosen.ok_or_else(|| {
            native::refused("original pair domains have no certified feasible candidate")
        })?;
        let mut pairs = 0usize;
        for group in &self.groups {
            for a in &group.first {
                for b in &group.second {
                    pairs += 1;
                    if pairs > 16384 {
                        return Err(native::refused(
                            "native pair ledger exceeds its domain-pair budget",
                        ));
                    }
                    k.check_cancelled("intersect::topology_solver::projected_pair_domains")?;
                    certify_domain_pair(k, objects, *a, *b, normal, bound.hi)?;
                }
            }
        }
        if let Some(target) = self.target {
            if !solver.tolerance_3d.is_finite() || solver.tolerance_3d < 0.0 {
                return Err(native::refused("invalid ledger retention slack"));
            }
            let cutoff = Ivl::from(target.abs()) - Ivl::from(solver.tolerance_3d);
            if bound.hi < cutoff.lo {
                solver.solutions.set_size(0);
                return Ok(());
            }
            if bound.lo < cutoff.hi {
                return Err(native::refused(
                    "native pair enclosure straddles external cutoff",
                ));
            }
        }
        solver.solutions.set_size(0);
        solver.solutions.add(&solution);
        solver.best_answer_so_far_sq =
            solution.start.solution_value * solution.start.solution_value;
        self.best = Some((solution, bound));
        Ok(())
    }
}
fn certify_domain_pair<O: SolutionObject>(
    k: &Kernel,
    objects: &dyn TopologyObjects<O>,
    a: Domain<'_, O>,
    b: Domain<'_, O>,
    normal: Vec3,
    upper: f64,
) -> Result<()> {
    match (a, b) {
        (
            Domain::Curve {
                curve, interval, ..
            },
            Domain::Face(face),
        )
        | (
            Domain::Face(face),
            Domain::Curve {
                curve, interval, ..
            },
        ) => {
            let endpoints = affine_curve_endpoints(curve, interval)?.ok_or_else(|| {
                native::refused("mixed non-affine curve/face domain has no certificate")
            })?;
            let net = projected_bound::face_net(k, objects, face)?;
            native::certify(
                k,
                &net,
                &objects.face_uv_domain(face)?,
                &endpoints,
                normal,
                upper,
            )
        }
        (Domain::Point(point, _), Domain::Face(face))
        | (Domain::Face(face), Domain::Point(point, _)) => {
            let net = projected_bound::face_net(k, objects, face)?;
            if projected_residual::maximum_is_bounded(
                k,
                &net,
                &objects.face_uv_domain(face)?,
                point,
                normal,
                upper,
            )? {
                Ok(())
            } else {
                Err(native::refused(
                    "original point/face domain exceeds retained enclosure",
                ))
            }
        }
        (
            Domain::Curve {
                curve: a,
                interval: ia,
                ..
            },
            Domain::Curve {
                curve: b,
                interval: ib,
                ..
            },
        ) => {
            let bound = p9::projected_pair_domain_upper(k, [a, b], [ia, ib], normal)?;
            if bound.is_none_or(|bound| bound <= upper) {
                Ok(())
            } else {
                Err(native::refused(
                    "original curve pair domain exceeds retained enclosure",
                ))
            }
        }
        (
            Domain::Curve {
                curve, interval, ..
            },
            Domain::Point(point, _),
        )
        | (
            Domain::Point(point, _),
            Domain::Curve {
                curve, interval, ..
            },
        ) => {
            let control =
                crate::nurbs::math::CPoint::new(point.x, point.y, point.z, crate::nurbs::math::NOW);
            let point_curve = Curve::from(BSplineCurve::from_data(
                3,
                1,
                3,
                1,
                &[control, control],
                &[0.0, 0.0, 1.0, 1.0],
            )?);
            let bound = p9::projected_pair_domain_upper(
                k,
                [curve, &point_curve],
                [interval, point_curve.natural_interval()],
                normal,
            )?;
            if bound.is_none_or(|bound| bound <= upper) {
                Ok(())
            } else {
                Err(native::refused(
                    "original curve/point domain exceeds retained enclosure",
                ))
            }
        }
        (Domain::Point(a, _), Domain::Point(b, _)) => {
            let bound = projected_bound::pair_metric(
                [a.x, a.y, a.z].map(Ivl::from),
                [b.x, b.y, b.z].map(Ivl::from),
                normal,
            );
            if bound.hi.is_finite() && bound.hi <= upper {
                Ok(())
            } else {
                Err(native::refused(
                    "original point pair exceeds retained enclosure",
                ))
            }
        }
        _ => Err(native::refused(
            "mixed topology pair domain has no native certificate",
        )),
    }
}

fn contains_point<O: SolutionObject>(
    objects: &dyn TopologyObjects<O>,
    domain: Domain<'_, O>,
    owner: O,
    solution: &Solution<O>,
    parameter: &mut usize,
) -> Result<bool> {
    Ok(match domain {
        Domain::Curve {
            curve,
            interval,
            owner: expected,
        } => {
            let t = solution.start.param(*parameter)?;
            *parameter += 1;
            expected == owner
                && t.is_finite()
                && interval.min() <= t
                && t <= interval.max()
                && curve.natural_interval().min() <= t
                && t <= curve.natural_interval().max()
        }
        Domain::Face(expected) => {
            let uv = crate::math::Vec2::new(
                solution.start.param(*parameter)?,
                solution.start.param(*parameter + 1)?,
            );
            *parameter += 2;
            expected == owner
                && uv.x.is_finite()
                && uv.y.is_finite()
                && objects.face_uv_domain(expected)?.contains_point(uv, 0.0)
        }
        Domain::Point(point, expected) => {
            expected == Some(owner) && objects.vertex_point(owner)? == point
        }
        Domain::Unknown => false,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn empty_cartesian_groups_do_not_accumulate_or_consume_domain_budget() {
        let mut ledger: Ledger<'_, crate::topology::brep_cache::CacheObject> = Ledger::default();
        for _ in 0..5000 {
            ledger
                .append(Vec::new(), vec![Domain::Point(Vec3::ZERO, None)])
                .unwrap();
        }
        assert!(ledger.groups.is_empty());
        assert_eq!(ledger.operands, 0);
    }
}
