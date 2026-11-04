//! A round turning at its radius law's corners, split there on the live session (not
//! upstream): the fork's `SmFilletExecutive::NoteRadiusCornerSurfaces`,
//! `SplitFilletFacesAtRadiusCorners` and `sm_ReplaceIsoCrvOnSurf`
//! (`source/SMLib/src/SmFilletExecutive.cpp`, tangency-labs/smlib@6eecf04, with 58d03ef's
//! copy of the pieces; WR-685's law-corner addendum, #2273).
//!
//! The circular section skins a round whose radius law has a corner piece by piece and
//! joins the pieces C0 there ([`crate::blend::fillet_standard_solver::p1::CircularCrossSectionFsg::piecewise_surface`]);
//! SMLib's check refuses a face (and the edges round it) with such a kink inside. The
//! merge gives the target copies of the fillet faces and drops the fillet B-rep, so the
//! surfaces are noted before it (with the pieces kept for them, where a corner dips or
//! peaks) and found again after it by their equality: each face on one is split at its
//! surface's C0 knots, its surface cut back to the side it lies on, and the new edges
//! along the turn given the iso curve of the surface they lie on, exactly.
//!
//! | Rust | Fork |
//! |------|------|
//! | `FilletSession::note_radius_corner_surfaces` (crate) | `SmFilletExecutive::NoteRadiusCornerSurfaces` |
//! | `FilletSession::split_fillet_faces_at_radius_corners` (crate) | `SmFilletExecutive::SplitFilletFacesAtRadiusCorners` |
//! | `replace_iso_crv_on_surf` (private) | `sm_ReplaceIsoCrvOnSurf` |
//!
//! Deviations: the noted copies and their pieces are one list of pairs
//! (`m_vRadiusCornerSurfaces`, `m_vRadiusCornerPieces`), owned by the session and emptied
//! by the split. `BuildLawCorners`, which the split calls with the pieces where a corner
//! dips or peaks, is not ported yet: the faces keep the joined surface there. No
//! recursion.

use super::trim_at_creases::with_target_env_of;
use super::{FilletSession, SolverObject};
use crate::blend::law_corners::{find_law_corner_pieces, LawCornerPieces};
use crate::curve::curve::{Curve, CurveKind};
use crate::error::{EntityRef, Result, SolidError};
use crate::geom::extent1d::Extent1d;
use crate::geom::extent2d::Extent2d;
use crate::geom::geom_utility::p2::ContinuityType;
use crate::geom::object::ObjectRef;
use crate::kernel::Kernel;
use crate::math::EFF_ZERO_SQRT;
use crate::nurbs::b_spline_curve_n_lib::TrimMode;
use crate::nurbs::b_spline_surface::SurfParamType;
use crate::surface::surface::Surface;
use crate::topology::brep::p3;
use crate::topology::edge::{self, NotificationStores};
use crate::topology::face::p1::{calculate_uv_domain_from_uv_trim_curves, shrink_geometry};
use crate::topology::{Brep, EdgeId};

const OP: &str = "blend::fillet_executive::session::law_corners";

/// A fillet surface noted before the merge (a copy), with its kept pieces where a law
/// corner dips or peaks.
pub(super) type RadiusCornerSurface = (Surface, Option<LawCornerPieces>);

impl FilletSession<'_> {
    /// `NoteRadiusCornerSurfaces` (the fork's 6eecf04, 58d03ef): copies of the fillet
    /// surfaces of the solvers whose radius has a corner, where those surfaces turn (a C0
    /// knot), each with a copy of the pieces kept for it (matched by its corners,
    /// [`find_law_corner_pieces`]) where a corner dips or peaks.
    ///
    /// # Errors
    ///
    /// Those of the solvers' corners, the geometries' surfaces, the knot search and the
    /// pieces' match.
    pub(super) fn note_radius_corner_surfaces(&mut self, k: &Kernel) -> Result<()> {
        let mut noted = Vec::new();
        for index in 0..self.objects.solvers.len() {
            // `GetRadiusCornerParameters` (virtual): the variable radius solver's
            let geoms = match self.objects.solvers.get(index).and_then(Option::as_ref) {
                Some(SolverObject::VariableRadius(vr)) => {
                    if vr.radius_corner_parameters(k)?.is_empty() {
                        continue;
                    }
                    vr.curve_based.solver.fillet_geoms()
                }
                _ => continue,
            };
            for geom in geoms {
                // (the surface is the intermediate body's face's by now)
                let Some(surface) = self.objects.fillet_surface(geom)? else {
                    continue;
                };
                let (u, v) = surface.find_knots_by_max_continuity(k, ContinuityType::C0, false)?;
                if u.is_empty() && v.is_empty() {
                    continue;
                }
                let copy = surface.clone();
                // the exact pieces where a law corner dips or peaks: those kept with the
                // same corners
                let base = self
                    .objects
                    .solvers
                    .get_mut(index)
                    .and_then(Option::as_mut)
                    .ok_or_else(|| SolidError::null_pointer(OP))?
                    .base_mut();
                let pieces = find_law_corner_pieces(k, &mut base.law_corner_pieces, &copy)?
                    .and_then(|i| base.law_corner_pieces.get(i))
                    .filter(|kept| kept.has_dip_or_peak())
                    .cloned();
                noted.push((copy, pieces));
            }
        }
        self.objects.radius_corner_surfaces.extend(noted);
        Ok(())
    }

    /// `SplitFilletFacesAtRadiusCorners` (the fork's 6eecf04): each face of the target on
    /// a surface noted by [`Self::note_radius_corner_surfaces`] split at its surface's C0
    /// knots (`SmBrep::SubdivideFaceAtDiscontinuities`), each resulting face's surface cut
    /// at the turns to the side the face lies on and the face shrunk back to it, and the
    /// split's edges given the iso curves they lie along. Empties the noted list.
    ///
    /// # Errors
    ///
    /// Those of the surface equality, the split, the knot search, the trims, the shrink and
    /// the curve replacement.
    pub(super) fn split_fillet_faces_at_radius_corners(&mut self, k: &mut Kernel) -> Result<()> {
        let noted = std::mem::take(&mut self.objects.radius_corner_surfaces);
        if noted.is_empty() {
            return Ok(());
        }
        // the keys a caller-named target recorded before the split name its entities
        // afterwards; the split's new edge curves are named by the session's identity, as
        // the merge's are, and rebound to the target's naming after it (as `merge_fillet`
        // does, tangency-labs/solidcore#423), also after an error
        let recorded = self.objects.target.geometry_owner_keys().clone();
        let result = self.split_noted_faces(k, &noted);
        let naming = self.objects.target.naming();
        let serial = self.objects.target_serial;
        let session_owner = |e| p3::topology_ref(serial, e);
        let restored =
            super::piecewise::normalize_session_owners(self.objects.target, session_owner, |e| {
                super::piecewise::naming_owner(naming, recorded.get(&e).copied(), e)
                    .unwrap_or_else(|| session_owner(e))
            });
        result.and(restored)
    }

    /// The body of [`Self::split_fillet_faces_at_radius_corners`] on the noted surfaces.
    fn split_noted_faces(&mut self, k: &mut Kernel, noted: &[RadiusCornerSurface]) -> Result<()> {
        let serial = self.objects.target_serial;
        let key_of = move |entity| p3::topology_ref(serial, entity);
        for face in crate::topology::brep::p1::faces(self.objects.target, None)? {
            // the merge gives the target an exact copy of each fillet surface
            let Some(surface) = crate::topology::face::p2::surface(self.objects.target, face)?
            else {
                continue;
            };
            let mut fillet = false;
            for (copy, _) in noted {
                if copy.approx_eq(k, surface)? {
                    fillet = true;
                    break;
                }
            }
            if !fillet {
                continue;
            }
            let mut pieces = with_target_env_of(&mut self.objects, k, |k, b, env, _| {
                let key_of = env.key_of;
                let curve_ref = |e: EdgeId, kind: i64| ObjectRef {
                    key: key_of(EntityRef::Edge(e)).key,
                    type_code: kind,
                };
                crate::topology::brep::p10::subdivide_face_at_discontinuities(
                    k,
                    b,
                    env,
                    &curve_ref,
                    face,
                    ContinuityType::C0,
                )
            })?;
            if !pieces.contains(&face) {
                pieces.push(face);
            }
            // Each face's surface is cut at the turns to the side the face lies on, so the
            // turn is its boundary, then back to the face (whose trims end within tolerance
            // of the turn).
            let target = &mut *self.objects.target;
            for &piece in &pieces {
                let bx =
                    calculate_uv_domain_from_uv_trim_curves(k, target, piece, true, false, false)?
                        .domain;
                let mid = bx.mid();
                let surface = crate::topology::face::p2::surface(target, piece)?
                    .ok_or_else(|| SolidError::null_pointer(OP))?;
                let domain = surface.natural_uv_domain_with_nurb(k)?;
                let (mut lo, mut hi) = (domain.min(), domain.max());
                let (u, v) = surface.find_knots_by_max_continuity(k, ContinuityType::C0, false)?;
                for knot in u {
                    if knot <= lo.x || knot >= hi.x {
                        continue;
                    }
                    if knot < mid.x {
                        lo.x = knot;
                    } else {
                        hi.x = knot;
                    }
                }
                for knot in v {
                    if knot <= lo.y || knot >= hi.y {
                        continue;
                    }
                    if knot < mid.y {
                        lo.y = knot;
                    } else {
                        hi.y = knot;
                    }
                }
                let mut side = Extent2d::from_bounds(lo.x, lo.y, hi.x, hi.y);
                if side != domain {
                    target
                        .face_mut(piece)?
                        .surface
                        .as_mut()
                        .ok_or_else(|| SolidError::null_pointer(OP))?
                        .trim_with_domain(k, &mut side)?;
                }
                shrink_geometry(k, target, piece, &key_of)?;
            }

            // The split's edges lie along iso lines of the surface as it was, as curves on
            // that surface (with its turn inside); each becomes that iso curve itself.
            let mut edges = Vec::new();
            for &piece in &pieces {
                for e in crate::topology::face::p4::edges(target, piece, None, false, None)? {
                    if !edges.contains(&e) {
                        edges.push(e);
                    }
                }
            }
            let super::SessionObjects {
                target,
                notification_records,
                attributes,
                ..
            } = &mut self.objects;
            let mut stores = NotificationStores {
                records: notification_records,
                attributes,
                key_of: &key_of,
            };
            for e in edges {
                replace_iso_crv_on_surf(k, target, &mut stores, e)?;
            }
            // (the fork's BuildLawCorners(sPieces, pieces) where a law corner dips or peaks,
            // 7babe17, is not built yet: the faces keep the joined surface there)
        }
        Ok(())
    }
}

/// `sm_ReplaceIsoCrvOnSurf` (the fork's 6eecf04): an edge whose curve is a curve on a
/// surface along an iso line of that surface gets that iso curve, exactly, as its curve
/// (its part between the edge's ends, its way round), where it runs between the edge's
/// own ends within the edge's tolerance; other edges are left as they are.
///
/// # Errors
///
/// Those of the evaluations, the iso curve, its trim and reversal, and
/// [`crate::topology::brep::p2::replace_curve`].
fn replace_iso_crv_on_surf(
    k: &mut Kernel,
    b: &mut Brep,
    stores: &mut NotificationStores<'_>,
    e: EdgeId,
) -> Result<()> {
    let Some(curve) = edge::curve(b, e)? else {
        return Ok(());
    };
    let CurveKind::CrvOnSurf(cos) = curve.kind() else {
        return Ok(());
    };
    // the edge's ends in the surface's parameters
    let ivl = edge::interval(b, e)?;
    let uv0 = cos.uv_curve().evaluate_point(k, ivl.min())?;
    let uv1 = cos.uv_curve().evaluate_point(k, ivl.max())?;
    let surface = cos.surface();
    let domain = surface.natural_uv_domain_with_nurb(k)?;
    let tol = EFF_ZERO_SQRT * (1.0 + domain.size().length());
    let (param, iso, from, to) = if (uv0.x - uv1.x).abs() <= tol {
        (SurfParamType::U, uv0.x, uv0.y, uv1.y)
    } else if (uv0.y - uv1.y).abs() <= tol {
        (SurfParamType::V, uv0.y, uv0.x, uv1.x)
    } else {
        return Ok(()); // not along an iso line
    };
    if (to - from).abs() <= tol {
        return Ok(());
    }
    let mut iso_curve: Curve = surface
        .create_iso_parametric_curve(k, param, iso, 0.0, None, false)?
        .curve;

    // the iso curve runs along the other parameter: keep the edge's part, its way round
    let natural = iso_curve.natural_interval();
    let mut part = Extent1d::from_min_max(
        natural.min().max(from.min(to)),
        natural.max().min(from.max(to)),
    );
    if part.length() <= tol {
        return Ok(());
    }
    iso_curve.trim(k, &mut part, TrimMode::Notify)?;
    if to < from {
        let whole = iso_curve.natural_interval();
        iso_curve.reverse_parameterization(k, &whole)?;
    }

    // it must run between the edge's own ends
    let end0 = curve.evaluate_point(k, ivl.min())?;
    let end1 = curve.evaluate_point(k, ivl.max())?;
    let iso_ivl = iso_curve.natural_interval();
    let iso0 = iso_curve.evaluate_point(k, iso_ivl.min())?;
    let iso1 = iso_curve.evaluate_point(k, iso_ivl.max())?;
    let edge_tol = edge::tolerance(b, e)?;
    if !iso0.closer_than(edge_tol, end0) || !iso1.closer_than(edge_tol, end1) {
        return Ok(());
    }
    crate::topology::brep::p2::replace_curve(k, b, e, iso_curve, stores)
}
