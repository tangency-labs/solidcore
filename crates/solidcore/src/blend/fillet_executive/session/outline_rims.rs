//! A lone round ends on the body's outline (not upstream; solidcore #2541, WR-729,
//! `fixes.fillet_round_ends_on_outline`): before its trace, a lone round's solver rolls its
//! ball on the surfaces tangent to it along the rims it touches
//! ([`crate::blend::fillet_outline::outline_supports`]).
//!
//! | Rust                                  | Upstream |
//! |---------------------------------------|----------|
//! | `FilletSession::end_round_on_outline` | none     |
//!
//! Deviations: none of upstream's; there is no recursion.

use super::{invalid, FilletObject, FilletSession, SolverObject};
use crate::blend::fillet_solver::p1::{
    FilletSurfaceGenerator, FilletSurfaceGeneratorType, OffsetSurfaceSlot,
};
use crate::blend::fillet_custom_solver::OffsetSurface;
use crate::blend::fillet_outline::outline_supports;
use crate::blend::fillet_solver::p5::is_concave_edge;
use crate::error::{Result, SolidError};
use crate::kernel::Kernel;
use crate::surface::offset_surface::SS_UNKNOWN;
use crate::topology::edgeuse::p1 as eu1;

const OP: &str = "blend::fillet_executive::session::outline_rims";

impl FilletSession<'_> {
    /// With `fixes.fillet_round_ends_on_outline` (and the cliff overflow, the refusal of
    /// rails leaving their faces and the exact round classes it builds on), a lone constant
    /// round of circular section
    /// whose ball would take material past a rim of its face rolls on the surfaces that make
    /// it end on the body's outline: each solver offset surface whose face it leaves is
    /// replaced by one of the same offset on the surface tangent to the ball along the rim
    /// ([`outline_supports`]), in the solver and in its fillet geometries. A round built
    /// together with others meets them in corners, which this construction does not make,
    /// so only a round on its own ends on the outline.
    ///
    /// # Errors
    /// Those of the topology queries and of [`outline_supports`].
    pub(super) fn end_round_on_outline(&mut self, k: &Kernel) -> Result<()> {
        let lone = k.fixes.fillet_round_ends_on_outline
            && k.fixes.fillet_cliff_overflow
            && k.fixes.fillet_refuse_rails_leaving_faces
            && k.fixes.fillet_exact_round_classes
            && self.executive.solvers.len() == 1;
        if !lone {
            return Ok(());
        }
        let Some(&FilletObject::Solver(id)) = self.executive.solvers.first() else {
            return Err(invalid("object is not a solver"));
        };
        // a constant round of circular section: its ball is one ball all along
        let object = self.objects.solver(id)?;
        let base = object.base();
        let circular = base.fsg.as_ref().map(FilletSurfaceGenerator::fillet_surface_generator_type)
            == Some(FilletSurfaceGeneratorType::Circular);
        if !matches!(object, SolverObject::ConstantRadius(_)) || !circular {
            return Ok(());
        }
        let tangency = base.tangency_tolerance;
        let uses = base.edgeuses;
        let mut offsets = [0.0; 2];
        for (offset, rail) in offsets.iter_mut().zip(0..2) {
            *offset = base
                .surface(rail)?
                .ok_or_else(|| SolidError::null_pointer(OP))?
                .surface
                .offset_distance();
        }
        let geoms = base.fillet_geoms();
        let target = &mut *self.objects.target;
        let [Some(u0), Some(u1)] = uses else {
            return Ok(());
        };
        let Some(edge) = eu1::edge(target, u0)? else {
            return Ok(());
        };
        // a convex round, which takes material (a concave one past its wall stays refused)
        if !crate::topology::brep::p2::is_manifold_solid(target)? || is_concave_edge(k, target, edge)? {
            return Ok(());
        }
        let mut faces = Vec::with_capacity(2);
        for u in [u0, u1] {
            let face = eu1::faceuse(target, u)?
                .map(|fu| eu1::faceuse_face(target, fu))
                .transpose()?
                .flatten()
                .ok_or_else(|| SolidError::null_pointer(OP))?;
            faces.push(face);
        }
        let faces = <[_; 2]>::try_from(faces).map_err(|_| SolidError::null_pointer(OP))?;
        let Some(supports) = outline_supports(k, target, edge, faces, offsets, tangency)? else {
            return Ok(());
        };
        for (rail, support) in supports.into_iter().enumerate() {
            let Some(support) = support else {
                continue;
            };
            let solver = self.objects.solver_mut(id)?.base_mut();
            let old = solver
                .surface(rail)?
                .ok_or_else(|| SolidError::null_pointer(OP))?;
            let (key, distance) = (old.key, old.surface.offset_distance());
            let base_id = old.surface.base_id;
            let surface = OffsetSurface::new(k, distance, support, base_id, SS_UNKNOWN)?;
            solver.set_surface(
                rail,
                Some(OffsetSurfaceSlot {
                    key,
                    surface: surface.clone(),
                }),
            )?;
            for &g in &geoms {
                if let Some(slot) = self.objects.pseudo.geom_mut(g)?.off_surfs.get_mut(rail) {
                    if slot.is_some() {
                        *slot = Some(surface.clone());
                    }
                }
            }
        }
        Ok(())
    }
}
