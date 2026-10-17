//! Rail classification and rollover processing for a fillet solver.
//!
//! Upstream: `source/SMLib/src/SmFilletSolver.cpp` (5714–5822, 7161–7899).
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`FilletSolver::test_roll_over`] | `SmFilletSolver::TestRollOver` |
//! | [`FilletSolver::process_roll_over`] | `SmFilletSolver::ProcessRollOver` |
//!
//! Deviations: pseudo-B-rep objects use the arenas from `blend.fillet_geom_p1`;
//! solver geometry handles are resolved explicitly into that arena. Resolvers distinguish
//! target geometry from pseudo-B-rep geometry and resolve corner types. Outputs
//! stay in the caller's lists, including on failure. The full-rollover recovery
//! stops at `ProcessFullRollOver` (p3); the multi-interval path stops at its first
//! missing binormal or `sm_FindRolloverSideEdge` (p2). Its subsequent topology
//! splitting, cliff creation, corner adjustment and surface rebuilding require
//! that helper's outputs and are deferred with it. `TestRollOver` stops at the real
//! `CurveOnClassify`'s `ReSet` dependency; after classification, short-end squeezing
//! and refinement stop at `ConnectIntervals` (`curve_class_p2`) and
//! `sm_RefineRailCurveClassification` (`fillet_solver_p3`). No recursion.

use super::p1::{FilletGeomRef, FilletSolver};
use crate::blend::fillet_geom::p1::{
    FilletBrep, FilletCornerId, FilletCornerType, FilletGeomId, FilletGeomType,
};
use crate::error::{ignored, EntityRef, Result, SolidError};
use crate::intersect::curve_class::p1::CurveClassification;
use crate::intersect::point_class::PointClassType;
use crate::kernel::Kernel;
use crate::query::face::p2::{curve_on_classify, CurveOnClassifyObjects};
use crate::surface::surf_of_revolution::SurfOfRevolution;
use crate::topology::{edge, edgeuse::p1 as eu, edgeuse::p2, Brep, EdgeId, EdgeuseId};

const OP: &str = "blend::fillet_solver::p5";

/// Resolve handles used by `TestRollOver`, without assuming the two B-reps' arena
/// indices identify the same objects. Curves returned here must be resolvable by
/// the inherited classification geometry resolver.
pub trait RolloverObjects: CurveOnClassifyObjects<EntityRef> {
    /// Identity of a rail's 3D curve in the pseudo-B-rep.
    fn rail_curve(&self, rail: EdgeId) -> EntityRef;
    /// Identity of the original face's UV curve on the rail edgeuse's mate.
    fn rail_mate_uv_curve(&self, mate: EdgeuseId) -> Option<EntityRef>;
    /// The actual virtual `GetCornerType` result for the specified corner.
    ///
    /// # Errors
    /// A resolver must report an unresolved corner, rather than invent its type.
    fn corner_type(&self, corner: FilletCornerId) -> Result<FilletCornerType>;
}

impl FilletSolver {
    /// Classify both rails against their original faces, squeeze tiny end
    /// intervals, refine their classifications and process rollover. Iterates a
    /// snapshot of the solver's geometry list, as upstream does. `geometry_id`
    /// resolves each solver-local handle into the caller's pseudo-B-rep arena.
    ///
    /// # Errors
    /// Missing rails/curves/mates are invalid inputs (upstream dereferences them).
    /// Classification and unported recovery errors propagate; completed earlier
    /// changes and the caller's adjusted-corner list are retained.
    pub fn test_roll_over(
        &mut self,
        k: &Kernel,
        target: &Brep,
        fillets: &mut FilletBrep,
        objects: &impl RolloverObjects,
        geometry_id: &impl Fn(FilletGeomRef) -> Result<FilletGeomId>,
        adjusted_corners: &mut Vec<FilletCornerId>,
    ) -> Result<()> {
        let geoms = self.fillet_geoms();
        for reference in geoms {
            let id = geometry_id(reference)?;
            for ii in 0..2 {
                let rail = fillets
                    .geom(id)?
                    .rail(ii)
                    .ok_or_else(|| SolidError::invalid_input(OP, "TestRollOver: missing rail"))?;
                let curve = edge::curve(fillets.brep(), rail)?.ok_or_else(|| {
                    SolidError::invalid_input(OP, "TestRollOver: missing rail curve")
                })?;
                let interval = curve.natural_interval();
                let primary = eu::primary_edgeuse(fillets.brep(), rail)?.ok_or_else(|| {
                    SolidError::invalid_input(OP, "TestRollOver: missing primary edgeuse")
                })?;
                let mate = eu::mate(fillets.brep(), primary)?
                    .ok_or_else(|| SolidError::invalid_input(OP, "TestRollOver: missing mate"))?;
                let mut cc = CurveClassification::with_curve_interval(
                    k,
                    objects.rail_curve(rail),
                    interval,
                    self.approx_tol_3d,
                );
                cc.set_uv_curve1(objects.rail_mate_uv_curve(mate));
                fillets.fillet_edge_mut(rail)?.set_curve_class(Some(cc));
                let face = fillets.fillet_edge(rail)?.original_face().ok_or_else(|| {
                    SolidError::invalid_input(OP, "TestRollOver: missing original face")
                })?;
                // Install the classification before calling the fallible operation:
                // upstream leaves that object attached even when classification fails.
                let cc = fillets
                    .fillet_edge_mut(rail)?
                    .curve_class
                    .as_mut()
                    .ok_or_else(|| SolidError::null_pointer(OP))?;
                curve_on_classify(k, target, face, objects, true, cc, true, false, true)?;
                let min_interval_size = 50.0 * self.approx_tol_3d;
                let count = cc.size();
                if count > 1 {
                    for jj in (0..count).rev() {
                        let interval = fillets
                            .fillet_edge(rail)?
                            .curve_class()
                            .ok_or_else(|| SolidError::null_pointer(OP))?
                            .sub_interval(jj)?;
                        let curve = edge::curve(fillets.brep(), rail)?.ok_or_else(|| {
                            SolidError::invalid_input(OP, "TestRollOver: missing rail curve")
                        })?;
                        let p1 = curve.evaluate_point(interval.min())?;
                        let p2 = curve.evaluate_point(interval.max())?;
                        if p1.distance(p2) > min_interval_size {
                            continue;
                        }
                        if jj == 0 || jj == count - 1 {
                            // UNPORTED(intersect.curve_class_p2): SmCurveClassification::ConnectIntervals
                            return Err(SolidError::unported(
                                "SmCurveClassification::ConnectIntervals",
                                "intersect.curve_class_p2",
                            ));
                        }
                    }
                }
                // This warning-only call cannot be silently discarded while its
                // implementation is missing: refinement changes the classification.
                // UNPORTED(blend.fillet_solver_p3): sm_RefineRailCurveClassification
                ignored(refine_rail_classification())?;
            }
            self.process_roll_over(
                k,
                target,
                fillets,
                id,
                &|corner| objects.corner_type(corner),
                adjusted_corners,
            )?;
        }
        Ok(())
    }

    /// Process an existing pair of rail classifications. Blend patches and
    /// single-interval rails that both remain classified need no repair.
    /// `corner_type` resolves the corner on rail 0's start vertex only when full
    /// rollover is needed; `adjusted_corners` is append-only upstream output.
    ///
    /// # Errors
    /// A missing classification or corner returns `SM_ERR_NULL_POINTER`. Unsupported
    /// full-rollover corner types return `SM_ERR`. Recovery through unavailable
    /// units returns `Unported`; invalid handles are reported by the arena.
    pub fn process_roll_over(
        &mut self,
        k: &Kernel,
        target: &Brep,
        fillets: &mut FilletBrep,
        id: FilletGeomId,
        corner_type: &impl Fn(FilletCornerId) -> Result<FilletCornerType>,
        adjusted_corners: &mut Vec<FilletCornerId>,
    ) -> Result<()> {
        let geom = fillets.geom(id)?;
        if geom.fillet_geom_type() == FilletGeomType::Blends {
            return Ok(());
        }
        let rail0 = geom
            .rail(0)
            .ok_or_else(|| SolidError::invalid_input(OP, "missing rail 0"))?;
        let cc0 = fillets
            .fillet_edge(rail0)?
            .curve_class()
            .ok_or_else(|| SolidError::null_pointer(OP))?;
        let rail1 = geom
            .rail(1)
            .ok_or_else(|| SolidError::invalid_input(OP, "missing rail 1"))?;
        let cc1 = fillets
            .fillet_edge(rail1)?
            .curve_class()
            .ok_or_else(|| SolidError::null_pointer(OP))?;
        let counts = [cc0.size(), cc1.size()];
        if counts == [1, 1] {
            let mid0 = cc0.at(0)?.mid.point_class();
            let mid1 = cc1.at(0)?.mid.point_class();
            if mid0 != PointClassType::Unknown && mid1 != PointClassType::Unknown {
                return Ok(());
            }
            let vertex = edge::vertex(fillets.brep(), rail0)?
                .ok_or_else(|| SolidError::invalid_input(OP, "rail without vertex"))?;
            let corner = fillets
                .fillet_vertex(vertex)?
                .fillet_corner()
                .ok_or_else(|| SolidError::null_pointer(OP))?;
            if !matches!(
                corner_type(corner)?,
                FilletCornerType::NByOneClosed | FilletCornerType::OneByOne
            ) {
                return Err(SolidError::sm_err(OP));
            }
            // Upstream tries rail 0 first, then rail 1 (including when both
            // are unknown). No mutations are fabricated by the unavailable call.
            // UNPORTED(blend.fillet_solver_p3): SmFilletSolver::ProcessFullRollOver
            return Err(SolidError::unported(
                "SmFilletSolver::ProcessFullRollOver",
                "blend.fillet_solver_p3",
            ));
        }
        let solver_eu = self
            .edgeuse(0)?
            .ok_or_else(|| SolidError::invalid_input(OP, "missing solver edgeuse"))?;
        let convex = p2::is_convex_radial_sector(k, target, solver_eu, 10, true, 0.0)?;
        let surface = fillets
            .geom(id)?
            .fillet_surface()
            .ok_or_else(|| SolidError::null_pointer(OP))?;
        let analytic = if surface.is_kind_of(crate::surface::surface::p2::B_SPLINE_SURFACE_TYPE)? {
            SurfOfRevolution::is_nurb_surface_surf_of_revolution(k, surface, 1.0)?.is_some()
        } else {
            false
        };
        // Copies are made before the segment walk, including for empty classes.
        for rail in [rail0, rail1] {
            let curve = edge::curve(fillets.brep(), rail)?
                .ok_or_else(|| SolidError::invalid_input(OP, "missing rail curve"))?;
            let _curve_copy = curve.copy(k)?;
            let primary = eu::primary_edgeuse(fillets.brep(), rail)?
                .ok_or_else(|| SolidError::invalid_input(OP, "missing primary edgeuse"))?;
            let (uv, _) = eu::uv_trim_curve(fillets.brep(), primary)?;
            let _uv_copy = uv.ok_or_else(|| SolidError::null_pointer(OP))?.copy(k)?;
        }
        let _surface_copy = surface.copy(k)?;
        let _ = (convex, analytic, adjusted_corners);
        if counts.contains(&0) {
            return Ok(());
        }
        // With positive counts outside [1,1], at least one rail reaches this
        // helper in the first segment. Its outputs govern all remaining repair.
        // UNPORTED(blend.fillet_solver_p2): sm_FindRolloverSideEdge
        Err(SolidError::unported(
            "sm_FindRolloverSideEdge",
            "blend.fillet_solver_p2",
        ))
    }
}

fn refine_rail_classification() -> Result<()> {
    // UNPORTED(blend.fillet_solver_p3): sm_RefineRailCurveClassification
    Err(SolidError::unported(
        "sm_RefineRailCurveClassification",
        "blend.fillet_solver_p3",
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blend::fillet_geom::p1::{FilletEdge, FilletGeom, FilletVertex};
    use crate::intersect::point_class::PointClassType;
    use crate::math::Vec3;

    fn classified_rails(k: &Kernel, fb: &mut FilletBrep) -> FilletGeomId {
        let mut geom = FilletGeom::new(k, None, FilletGeomType::Default).unwrap();
        for i in 0..2 {
            let a = fb.brep_mut().add_vertex(Vec3::ZERO).unwrap();
            let b = fb.brep_mut().add_vertex(Vec3::X).unwrap();
            let edge = fb.brep_mut().add_edge(a, b).unwrap();
            let vu = fb
                .brep_mut()
                .add_vertexuse_record(crate::topology::Vertexuse {
                    vertex: Some(a),
                    ..crate::topology::Vertexuse::default()
                })
                .unwrap();
            let eu = fb
                .brep_mut()
                .add_edgeuse_record(crate::topology::Edgeuse {
                    edge: Some(edge),
                    vertexuse: Some(vu),
                    ..crate::topology::Edgeuse::default()
                })
                .unwrap();
            fb.brep_mut().edge_mut(edge).unwrap().edgeuses.push(eu);

            fb.insert_fillet_vertex(a, FilletVertex::new(k, Some(FilletCornerId(42 + i as u32))));
            let mut record = FilletEdge::with_corner(None);
            let mut cc = CurveClassification::new(k);
            cc.set_size(k, 1).unwrap();
            cc.at_mut(0)
                .unwrap()
                .mid
                .set_class_object(PointClassType::Edge, Some(EntityRef::Edge(edge)));
            record.set_curve_class(Some(cc));
            fb.insert_fillet_edge(edge, record);
            geom.set_rail(i, Some(edge));
        }
        fb.add_geom(geom)
    }

    #[test]
    fn resolver_is_lazy_and_failure_preserves_caller_output() {
        let k = Kernel::default();
        let mut fb = FilletBrep::new();
        let id = classified_rails(&k, &mut fb);
        let mut solver = FilletSolver::new(1e-5, 1e-3, 1e-3);
        let target = Brep::new();
        let mut adjusted = vec![FilletCornerId(9), FilletCornerId(9)];
        solver
            .process_roll_over(
                &k,
                &target,
                &mut fb,
                id,
                &|_| panic!("known rails must not resolve a corner"),
                &mut adjusted,
            )
            .unwrap();
        let rail = fb.geom(id).unwrap().rail(1).unwrap();
        fb.fillet_edge_mut(rail)
            .unwrap()
            .curve_class
            .as_mut()
            .unwrap()
            .at_mut(0)
            .unwrap()
            .mid
            .set_class_object(PointClassType::Unknown, None);
        let refusal = SolidError::invalid_input(OP, "unresolved corner");
        let result = solver.process_roll_over(
            &k,
            &target,
            &mut fb,
            id,
            &|corner| {
                assert_eq!(corner, FilletCornerId(42));
                Err(refusal.clone())
            },
            &mut adjusted,
        );
        assert_eq!(result, Err(refusal));
        assert_eq!(adjusted, [FilletCornerId(9), FilletCornerId(9)]);
        let result = solver.process_roll_over(
            &k,
            &target,
            &mut fb,
            id,
            &|_| Ok(FilletCornerType::OneByOne),
            &mut adjusted,
        );
        assert!(matches!(
            result,
            Err(SolidError::Unported {
                unit: "blend.fillet_solver_p3",
                ..
            })
        ));
    }
}
