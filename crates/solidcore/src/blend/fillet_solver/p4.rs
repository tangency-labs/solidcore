//! Fillet geometry of a fillet solver: the ordered intersection points along the
//! filleted edge, the analytic attempt, the offset-surface trace of the fillet centre
//! line (with cusp recovery), the split and self-intersection post-processing, and
//! the point evaluation the fillet intersector steps through.
//!
//! Upstream: `source/SMLib/src/SmFilletSolver.cpp` (1495–2673, 3544–3768),
//! `source/SMLib/inc/SmFilletSolver.h`.
//!
//! | Rust                                        | Upstream                                          |
//! |---------------------------------------------|---------------------------------------------------|
//! | [`FilletSolver::calc_fillet_geom`]          | `SmFilletSolver::CalcFilletGeom`                  |
//! | [`FilletSolver::compute_point_values`]      | `SmFilletSolver::ComputePointValues`              |
//! | [`FilletGeomResolver`]                      | the `SmFilletGeom*` entries of `m_vFilletGeoms` (pointer identity) |
//! | `TsRef`, `TsectStore` (private)             | `SmTArray<SmTsectPnt*> sTsectPnts` / `sDeletePnts` (pointers into `SmFilletVertexuse::m_vTsectPnt` or to `new SmTsectPnt`) |
//! | `get_or_create_uv_trim_curve` (private)     | `SmEdgeuse::GetOrCreateUVTrimCurve` (existing-curve path only; unit `topology.edgeuse_p3`) |
//!
//! Deviations:
//! - **Object model.** The solver's virtuals are the [`FilletSolverOps`] trait object
//!   passed in (`this`); `SmContext` arguments are the [`Kernel`]. The filleted edge,
//!   its edge uses and their UV trim curves belong to the *target* [`Brep`]; rails,
//!   fillet vertex uses and fillet geometries to the pseudo B-rep [`FilletBrep`]. The
//!   solver's opaque [`FilletGeomRef`] handles are resolved (and new geometries
//!   registered) through a [`FilletGeomResolver`].
//! - **Intersection point list.** `sTsectPnts` mixes pointers into fillet vertex uses
//!   (`SmFilletVertexuse::GetTsectPnt()`, the rail-end points) and points allocated by
//!   this function (G1-knot points, the cusp-recovery point). Here the list holds
//!   `TsRef` handles: a vertex use id, or an index into a local append-only store
//!   (`sDeletePnts`, which lives as long as the call). Every in-place write upstream
//!   makes through those pointers (the `smos_CleanUpNoise` of the rail-end UVs, the
//!   `m_dCurveParameter = -SM_BIG_DOUBLE` sort mark, everything `DoPointIntersection`
//!   writes) goes back to the storage the handle names; `DoPointIntersection` works on
//!   a snapshot that is written back point by point.
//! - **Function static.** `ComputePointValues`' `static double sdStaticStep` (set once
//!   per process by the first finite-difference evaluation) is the
//!   [`Kernel::fillet_point_step`] latch: per kernel instead of per process. Difftest
//!   runners seed it from the C++ process's value.
//! - **Intersector.** `SmFilletIntersector sFI(*pOffsetSrf1, ..)` traces copies of the
//!   solver's offset surfaces seen as [`Surface`]s (the intersector borrows its
//!   surfaces while the solver stays mutable).
//! - **Unported callees** (each a private stub at the end of this file returning
//!   [`SolidError::Unported`]): `CalcCylinderFilletGeom`, `CalcTorusFilletGeom`,
//!   `RecordFilletError` (`blend.fillet_solver_p2`); `sm_FindFilletTsectPnt`,
//!   `StepSolve` (`blend.fillet_solver_p3`); `SmFilletGeom::TopologySplit`
//!   (`blend.fillet_geom_p3`); `SmFilletGeom::CalcBlendingGeom`
//!   (`blend.fillet_geom_p4`); `SmFilletEdge::CalcCrossSection`
//!   (`blend.fillet_geom_p5`); `SmEdgeuse::EvaluateBinormal` and the creating path of
//!   `SmEdgeuse::GetOrCreateUVTrimCurve` (`topology.edgeuse_p3`). DISPATCH GAP: the
//!   base `SmSurfaceIntersector::DoPointIntersection` does not dispatch to the
//!   `SmFilletIntersector` overrides (`intersect.surface_intersector_p4`), so the trace
//!   returns [`SolidError::Unported`] there.
//! - **Iteration.** No recursion. Loops are bounded: two surfaces, two rails of two
//!   edge uses, the G1 knots and points, ten trace passes, twenty cusp skips, three
//!   drop tolerances.

// The skeleton's shared state and helpers are used as the chunks C1-C6 land; the
// integrator of the last chunk removes this line (it then warns as unfulfilled).
#![expect(dead_code, reason = "blend.fillet_solver_p4 chunks C1-C6 pending")]

use super::p1::{FilletGeomRef, FilletSolver, FilletSolverOps};
use crate::blend::fillet_custom_solver::OffsetSurface;
use crate::blend::fillet_geom::p1::{FilletBrep, FilletErrorType, FilletGeomId, FilletGeomType};
use crate::blend::fillet_intersector::{FilletIntersector, StepSolve as StepSolveFlags};
use crate::curve::curve::{Curve, CurveKind};
use crate::error::{Result, SolidError};
use crate::geom::extent1d::Extent1d;
use crate::geom::extent2d::Extent2d;
use crate::geom::vector2d::Point2;
use crate::intersect::surface_drop_lib::p1::SM_SUCCESS;
use crate::intersect::surface_intersector::p1::TsectPnt;
use crate::kernel::Kernel;
use crate::math::{Point3, Vec3};
use crate::nurbs::b_spline_curve::BSplineCurve;
use crate::nurbs::knots::Side;
use crate::surface::surface::p4::NoObject;
use crate::surface::surface::p7::drop_vectors_at;
use crate::surface::surface::{Surface, SurfaceKind};
use crate::topology::edgeuse::p1 as eu;
use crate::topology::{Brep, EdgeId, EdgeuseId, VertexuseId};

const OP: &str = "blend::fillet_solver::p4";

fn invalid(reason: &'static str) -> SolidError {
    SolidError::invalid_input(OP, reason)
}

// =====================================================================================
// SHARED: data model (sTsectPnts), state passed between the chunks, public resolver
// =====================================================================================

/// Resolves the solver's [`FilletGeomRef`] handles (`m_vFilletGeoms`) into the pseudo
/// B-rep's [`FilletGeomId`]s, and registers the geometries `CalcFilletGeom` adds
/// (`m_vFilletGeoms.Add`, after `TopologySplit` or a self-intersection trace).
pub trait FilletGeomResolver {
    /// The pseudo-B-rep geometry of a solver handle.
    ///
    /// # Errors
    ///
    /// A resolver reports an unknown handle rather than invent one.
    fn geometry_id(&self, reference: FilletGeomRef) -> Result<FilletGeomId>;

    /// Records that the solver handle `reference` (new in `m_vFilletGeoms`) is the
    /// pseudo-B-rep geometry `id`.
    ///
    /// # Errors
    ///
    /// A resolver reports a handle it cannot record.
    fn register_geometry(&mut self, reference: FilletGeomRef, id: FilletGeomId) -> Result<()>;
}

/// The map `reference.index -> id`.
impl FilletGeomResolver for Vec<FilletGeomId> {
    fn geometry_id(&self, reference: FilletGeomRef) -> Result<FilletGeomId> {
        self.get(reference.index)
            .copied()
            .ok_or_else(|| invalid("unknown fillet geometry handle"))
    }

    fn register_geometry(&mut self, reference: FilletGeomRef, id: FilletGeomId) -> Result<()> {
        if reference.index == self.len() {
            self.push(id);
            return Ok(());
        }
        let slot = self
            .get_mut(reference.index)
            .ok_or_else(|| invalid("fillet geometry handle past the end"))?;
        *slot = id;
        Ok(())
    }
}

/// An entry of `sTsectPnts` (an `SmTsectPnt*`): which storage the point lives in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TsRef {
    /// `&pVU->GetTsectPnt()`: the point stored in a pseudo-B-rep fillet vertex use.
    Vertexuse(VertexuseId),
    /// A point this call allocated (`new SmTsectPnt`, owned by `sDeletePnts`): an index
    /// into [`TsectStore`].
    Local(usize),
}

/// `sDeletePnts`: the points this call allocates. Append-only, so a [`TsRef::Local`]
/// stays valid for the whole call (as upstream's pointers until `sObjs` deletes them).
#[derive(Debug, Clone, Default)]
struct TsectStore {
    local: Vec<TsectPnt>,
}

impl TsectStore {
    /// `new SmTsectPnt` + `sDeletePnts.Add`: stores `point` and returns its handle.
    fn alloc(&mut self, point: TsectPnt) -> TsRef {
        self.local.push(point);
        TsRef::Local(self.local.len() - 1)
    }

    /// `*pTsectPnt`.
    fn get<'a>(&'a self, fillets: &'a FilletBrep, r: TsRef) -> Result<&'a TsectPnt> {
        match r {
            TsRef::Vertexuse(vu) => Ok(fillets.fillet_vertexuse(vu)?.tsect_pnt()),
            TsRef::Local(i) => self
                .local
                .get(i)
                .ok_or_else(|| invalid("local intersection point out of range")),
        }
    }

    /// `*pTsectPnt`, for an in-place write (it lands in the vertex use for a rail-end
    /// point, as upstream's write through the pointer does).
    fn get_mut<'a>(
        &'a mut self,
        fillets: &'a mut FilletBrep,
        r: TsRef,
    ) -> Result<&'a mut TsectPnt> {
        match r {
            TsRef::Vertexuse(vu) => Ok(fillets.fillet_vertexuse_mut(vu)?.tsect_pnt_mut()),
            TsRef::Local(i) => self
                .local
                .get_mut(i)
                .ok_or_else(|| invalid("local intersection point out of range")),
        }
    }

    /// Copies of the points `refs` names, in order (for a callee that takes the array
    /// of points; write them back with [`TsectStore::write_back`]).
    fn snapshot(&self, fillets: &FilletBrep, refs: &[TsRef]) -> Result<Vec<TsectPnt>> {
        refs.iter()
            .map(|r| self.get(fillets, *r).cloned())
            .collect()
    }

    /// Writes `points` (a [`TsectStore::snapshot`] of `refs`, possibly changed) back
    /// to the storage each handle names, in order (a handle listed twice gets its last
    /// value; upstream's two pointers would alias one point).
    fn write_back(
        &mut self,
        fillets: &mut FilletBrep,
        refs: &[TsRef],
        points: &[TsectPnt],
    ) -> Result<()> {
        if refs.len() != points.len() {
            return Err(invalid("intersection point snapshot length changed"));
        }
        for (r, p) in refs.iter().zip(points) {
            *self.get_mut(fillets, *r)? = p.clone();
        }
        Ok(())
    }
}

/// The objects `CalcFilletGeom` reaches besides the solver.
struct Cx<'a> {
    /// The kernel (`m_crContext`).
    k: &'a Kernel,
    /// The target B-rep: the filleted edge, its edge uses and UV trim curves.
    target: &'a mut Brep,
    /// The pseudo B-rep: rails, fillet vertex uses (rail-end points), geometries.
    fillets: &'a mut FilletBrep,
    /// `m_vFilletGeoms` handle resolution.
    geoms: &'a mut dyn FilletGeomResolver,
}

/// The locals of `CalcFilletGeom` that live across the chunks.
#[derive(Debug, Clone)]
struct CalcState {
    // ---- C1: setup (:1499-1553) ------------------------------------------------------
    /// `pEdge = GetEdgeuse(0)->GetEdge()` (target B-rep).
    edge: EdgeId,
    /// `pEdgeCurve`: the edge's curve (a kind derived from `SmBSplineCurve`, checked
    /// with [`bspline_curve_cast`]); virtual calls go through [`Curve`].
    edge_curve: Curve,
    /// `sIvl = pEdge->GetInterval()`.
    ivl: Extent1d,
    /// `sPVSt`: position and first derivative at `sIvl.GetMin()`.
    pv_start: [Vec3; 2],
    /// `sPVEn`: position and first derivative at `sIvl.GetMax()`.
    pv_end: [Vec3; 2],
    /// `pBaseSrf1`, `pBaseSrf2`: copies of `GetSurface(i)->GetBaseSurface()` taken at
    /// the start (upstream holds the pointers; nothing in the call replaces a base).
    base_surfaces: [Surface; 2],
    /// `pFilletGeom0 = m_vFilletGeoms[0]`.
    geom0: FilletGeomId,
    /// `bSplitOriginalFillet = pFilletGeom0->GetSplitFlag()`.
    split_original_fillet: bool,
    /// `pSplitTsectPnt` (`NULL` until a G1-knot point at the split parameter is made).
    split_tsect_pnt: Option<TsRef>,
    // ---- C1/C2: the point list (:1582-1920) ------------------------------------------
    /// `sDeletePnts` (the points this call allocates).
    store: TsectStore,
    /// `sTsectPnts`.
    tsect_pnts: Vec<TsRef>,
    /// `bDoSorting`.
    do_sorting: bool,
    // ---- C3: trace setup (:2028-2090) ------------------------------------------------
    /// The trace locals C3 sets and C4/C5 use.
    trace: TraceLocals,
}

/// The locals of the trace (:2028-2090), set by C3.
#[derive(Debug, Clone)]
struct TraceLocals {
    /// `sDomain1`, `sDomain2` after `SetupOffsetExtension`.
    domains: [Extent2d; 2],
    /// `dThisApproxTol3d = GetThisApproxTol3d()`.
    approx_tol_3d: f64,
    /// `dAngleTol = GetThisAngTolRad()`.
    angle_tol: f64,
    /// `sDir`: the trace start direction.
    dir: Vec3,
    /// `bReverseTracing`.
    reverse_tracing: bool,
    /// `bExtendBefore` (the local; `m_bExtendBefore` is the solver's).
    extend_before: bool,
    /// `bExtendAfter` (the local).
    extend_after: bool,
}

impl CalcState {
    /// The state after setup (:1499-1553): an empty point list, sorting on, and the
    /// trace locals at neutral values until C3 sets them.
    fn new(
        edge: EdgeId,
        edge_curve: Curve,
        ivl: Extent1d,
        pv_start: [Vec3; 2],
        pv_end: [Vec3; 2],
        base_surfaces: [Surface; 2],
        geom0: FilletGeomId,
        split_original_fillet: bool,
    ) -> Self {
        let domain = base_surfaces[0].natural_uv_domain();
        Self {
            edge,
            edge_curve,
            ivl,
            pv_start,
            pv_end,
            base_surfaces,
            geom0,
            split_original_fillet,
            split_tsect_pnt: None,
            store: TsectStore::default(),
            tsect_pnts: Vec::new(),
            do_sorting: true,
            trace: TraceLocals {
                domains: [domain, domain],
                approx_tol_3d: 0.0,
                angle_tol: 0.0,
                dir: pv_start[1],
                reverse_tracing: false,
                extend_before: false,
                extend_after: false,
            },
        }
    }

    /// `*sTsectPnts[i]`.
    fn point<'a>(&'a self, fillets: &'a FilletBrep, i: usize) -> Result<&'a TsectPnt> {
        let r = *self
            .tsect_pnts
            .get(i)
            .ok_or_else(|| invalid("intersection point index out of range"))?;
        self.store.get(fillets, r)
    }

    /// `*sTsectPnts[i]`, for an in-place write.
    fn point_mut<'a>(
        &'a mut self,
        fillets: &'a mut FilletBrep,
        i: usize,
    ) -> Result<&'a mut TsectPnt> {
        let r = *self
            .tsect_pnts
            .get(i)
            .ok_or_else(|| invalid("intersection point index out of range"))?;
        self.store.get_mut(fillets, r)
    }
}

/// One pass of the `lWhichSurf` loop after collection (:1593-1752): the UV trim curve
/// of `GetEdgeuse(lWhichSurf)` and `sAllTsectPnts`.
#[derive(Debug, Clone)]
struct SurfPass {
    /// `pUVCurveFilletedEdge` (a copy of the edge use's trim curve).
    uv_curve: Curve,
    /// `sAllTsectPnts`.
    all: Vec<TsRef>,
}

/// The fillet intersector of the trace: over copies of the offset surfaces, solution
/// objects unused, the solver's geometry handles.
type TraceIntersector<'s> = FilletIntersector<'s, 's, NoObject, FilletGeomRef>;

// =====================================================================================
// PUBLIC API: orchestration (keep the upstream control flow here; chunk bodies below)
// =====================================================================================

impl FilletSolver {
    /// `CalcFilletGeom()`: the fillet geometry of the solver's sector. Collects the
    /// rail-end points of the first fillet geometry, the points at the filleted edge's
    /// G1 knots (and at the split parameter), orders them along a UV trim curve of
    /// the edge, tries an analytic (cylinder / torus) fillet, and otherwise traces the
    /// offset-surface intersection through the points (recovering from cusps), then
    /// splits the geometry when asked, finishes self-intersection patches and undoes a
    /// reversed trace. `ops` is the solver (its virtuals); `target` holds the filleted
    /// edge; `fillets` the pseudo B-rep; `geoms` resolves `m_vFilletGeoms`.
    ///
    /// # Errors
    ///
    /// Upstream's statuses (`SM_ERR`, `SM_ERR_NULL_POINTER`, the callees');
    /// [`SolidError::InvalidInput`] where upstream dereferences `NULL`;
    /// [`SolidError::Unported`] from unported callees.
    pub fn calc_fillet_geom(
        ops: &mut dyn FilletSolverOps,
        k: &Kernel,
        target: &mut Brep,
        fillets: &mut FilletBrep,
        geoms: &mut dyn FilletGeomResolver,
    ) -> Result<()> {
        let mut cx = Cx {
            k,
            target,
            fillets,
            geoms,
        };

        // :1495-1553 SetExtensionFlags, edge, curve end points, surfaces, geometry 0
        let mut st = c1_setup(ops, &mut cx)?;

        // :1582-1867 build the ordered sTsectPnts on one of the two edge uses
        for which_surf in 0..2 {
            st.tsect_pnts.clear();
            // :1587-1752 trim curve, closed test, split parameter, G1 points, rail ends
            let Some(pass) = c1_collect_points(ops, &mut cx, &mut st, which_surf)? else {
                // closed UV trim curve on the first edge use: try the other one
                continue;
            };
            if st.do_sorting {
                // :1754-1858 drop to the trim curve and insertion-sort
                c2_sort_points(ops, &mut cx, &mut st, which_surf, &pass)?;
                if st.tsect_pnts.len() > 1 {
                    break;
                }
            } else {
                // :1860-1863 both trim curves closed: sTsectPnts.Append(sAllTsectPnts)
                st.tsect_pnts.extend(pass.all.iter().copied());
            }
        }

        // :1872-1920 no points; drop near-duplicate end points
        c2_clean_ends(ops, &mut cx, &mut st)?;

        // :1955-2011 analytic fillet
        if c3_try_analytic(ops, &mut cx, &mut st)? {
            return Ok(());
        }

        // :2028-2034 domains (and the intersector's surfaces)
        let surfaces = c3_trace_setup(ops, &mut cx, &mut st)?;
        {
            // :2037-2052 SmFilletIntersector sFI(..), extension distance
            let mut fi = c3_make_intersector(ops, &mut cx, &st, &surfaces)?;
            // :2054-2090 trace locals, reversal at a tangent last point
            c3_trace_start(ops, &mut cx, &mut st)?;
            // :2093-2467 while (!bFoundAllPoints)
            c4_trace(ops, &mut cx, &mut st, &mut fi)?;
        }

        // :2470-2671 split, self-intersection patches, undo the reversal
        c5_post_process(ops, &mut cx, &mut st)
    }

    /// `ComputePointValues(aUVValues, crUVDomain1, crUVDomain2, dCurveTraceDirection,
    /// rTsectPnt, pOptPreviousPnt)`: the surface values of `point` at the rail UV
    /// points, the centre-line tangent (cross product of the normals, or finite
    /// differences of solver steps for variable offsets and tangent points), the point
    /// type and the UV tangents. `ops` is the solver (its virtuals); the finite
    /// difference step is the kernel's [`Kernel::fillet_point_step`] latch.
    ///
    /// # Errors
    ///
    /// Those of `ComputeSurfaceValues`, `StepSolve` and the evaluations; `SM_ERR` when
    /// neither step succeeds or the stepped normals are parallel; a zero tangent's
    /// `Unitize` status; [`SolidError::Unported`] from unported callees.
    pub fn compute_point_values(
        ops: &mut dyn FilletSolverOps,
        k: &Kernel,
        uv_values: [Point2; 2],
        uv_domain1: &Extent2d,
        uv_domain2: &Extent2d,
        curve_trace_direction: f64,
        point: &mut TsectPnt,
        previous: Option<&TsectPnt>,
    ) -> Result<()> {
        c6_compute_point_values(
            ops,
            k,
            uv_values,
            uv_domain1,
            uv_domain2,
            curve_trace_direction,
            point,
            previous,
        )
    }
}

// =====================================================================================
// SHARED: helpers used by more than one chunk
// =====================================================================================

/// `GetSurface(rail)` (an `SmOffsetSurface*`; upstream dereferences it unchecked).
fn offset_surface(ops: &dyn FilletSolverOps, rail: usize) -> Result<&OffsetSurface> {
    ops.solver()
        .surface(rail)?
        .map(|slot| &slot.surface)
        .ok_or_else(|| invalid("no offset surface"))
}

/// `GetSurface(rail)`, mutable (offset distance changes, stored extensions).
fn offset_surface_mut(ops: &mut dyn FilletSolverOps, rail: usize) -> Result<&mut OffsetSurface> {
    ops.solver_mut()
        .surface_mut(rail)?
        .map(|slot| &mut slot.surface)
        .ok_or_else(|| invalid("no offset surface"))
}

/// `GetSurface(rail)->EvaluatePoint(uv, pnt)` on the solver's own surface (an
/// extended base it builds is stored, as upstream's `const` cast does).
fn offset_evaluate_point(ops: &mut dyn FilletSolverOps, rail: usize, uv: Point2) -> Result<Point3> {
    offset_surface_mut(ops, rail)?
        .surface
        .evaluate_point_extending(uv)
}

/// `GetSurface(rail)->EvaluateNormal(uv, TRUE, TRUE, n)` on the solver's own surface.
fn offset_evaluate_normal(ops: &mut dyn FilletSolverOps, rail: usize, uv: Point2) -> Result<Vec3> {
    offset_surface_mut(ops, rail)?
        .surface
        .evaluate_normal_extending(uv, Side::Left, Side::Left)
}

/// `SM_CAST_PTR(SmBSplineCurve, pCurve)`: the B-spline part of the kinds derived from
/// `SmBSplineCurve` (`SmLine`, `SmEllipse`/`SmCircle`, `SmParabola`, `SmHyperbola`),
/// `None` (upstream `NULL`) for the others.
const fn bspline_curve_cast(curve: &Curve) -> Option<&BSplineCurve> {
    match curve.kind() {
        CurveKind::BSpline(b) => Some(b),
        CurveKind::Line(l) => Some(l.base()),
        CurveKind::Ellipse(e) => Some(e.base()),
        CurveKind::Hyperbola(h) => Some(h.base()),
        CurveKind::Parabola(p) => Some(p.base()),
        CurveKind::Hermite(_)
        | CurveKind::Projected(_)
        | CurveKind::Iso(_)
        | CurveKind::Offset(_)
        | CurveKind::Conic(_)
        | CurveKind::CrvInVolume(_)
        | CurveKind::CrvOnSurf(_)
        | CurveKind::Composite(_) => None,
    }
}

/// `SM_CAST_PTR(SmBSplineSurface, pSurface) != NULL`: the kinds derived from
/// `SmBSplineSurface` (`SmPlane`, `SmSurfOfRevolution`).
const fn is_bspline_surface(surface: &Surface) -> bool {
    matches!(
        surface.kind(),
        SurfaceKind::BSpline(_) | SurfaceKind::Plane(_) | SurfaceKind::SurfOfRevolution(_)
    )
}

/// `m_vFilletGeoms`, resolved.
fn fillet_geom_ids(
    ops: &dyn FilletSolverOps,
    geoms: &dyn FilletGeomResolver,
) -> Result<Vec<FilletGeomId>> {
    ops.solver()
        .fillet_geoms
        .iter()
        .map(|r| geoms.geometry_id(*r))
        .collect()
}

/// `m_vFilletGeoms.Add(pNewFilletGeom)`: a new solver handle (one past the largest
/// index in use) for `id`, registered with the resolver.
fn add_fillet_geom(
    ops: &mut dyn FilletSolverOps,
    geoms: &mut dyn FilletGeomResolver,
    id: FilletGeomId,
) -> Result<FilletGeomRef> {
    let solver = ops.solver_mut();
    let index = solver
        .fillet_geoms
        .iter()
        .map(|r| r.index + 1)
        .max()
        .unwrap_or(0);
    let reference = FilletGeomRef { index };
    geoms.register_geometry(reference, id)?;
    solver.fillet_geoms.push(reference);
    Ok(reference)
}

/// The guess point of a centre-line point at parameter `t` of the filleted edge
/// (B2 at :1641-1657, F2 at :2395-2425): `GetEdgeuse(0)->EvaluateBinormal(t, FALSE,
/// ..)`, the binormal scaled to `|GetSurface(1)->GetOffsetDistance()|` (then halved
/// when `half`: `sBinVec*dOffset*0.5`), dropped into the parameter space of
/// `base0` (`pBaseSrf1->DropVectors(uv, TRUE, TRUE, ..)`) from the binormal's UV
/// point, evaluated on `GetSurface(0)`. Returns the guess point and the edge-use
/// tangent (`sTangent`).
fn binormal_guess_point(
    ops: &mut dyn FilletSolverOps,
    cx: &mut Cx<'_>,
    base0: &Surface,
    t: f64,
    half: bool,
) -> Result<(Point3, Vec3)> {
    let eu0 = ops
        .solver()
        .edgeuse(0)?
        .ok_or_else(|| invalid("no edge use 0"))?;
    let b = evaluate_binormal(cx.k, cx.target, eu0, t)?;
    let offset = offset_surface(ops, 1)?.offset_distance().abs();
    let mut bin_vec = b.binormal * offset;
    if half {
        bin_vec *= 0.5;
    }
    let uv_dir = drop_vectors_at(base0, b.uv, Side::Left, Side::Left, &[bin_vec])?
        .first()
        .copied()
        .ok_or_else(|| invalid("short drop result"))?;
    let uv = b.uv + uv_dir;
    let pnt = offset_evaluate_point(ops, 0, uv)?;
    Ok((pnt, b.tangent))
}

// =====================================================================================
// P4-CHUNK(C1) BEGIN: setup and point collection, SmFilletSolver.cpp:1495-1752
// (chunk-local imports go inside the functions: never edit the shared `use` block)
// =====================================================================================

/// :1495-1553: `SetExtensionFlags()`, the filleted edge, its curve (B-spline kinds
/// only, `NER`), interval and end-point derivatives, the base surfaces, geometry 0 and
/// its split flag.
fn c1_setup(ops: &mut dyn FilletSolverOps, cx: &mut Cx<'_>) -> Result<CalcState> {
    let _ = (ops, cx);
    // P4-CHUNK(C1): SmFilletSolver.cpp:1495-1553
    Err(SolidError::Unported {
        upstream: "SmFilletSolver::CalcFilletGeom (p4 chunk C1 pending)",
        unit: "blend.fillet_solver_p4",
    })
}

/// :1587-1752, one pass of the `lWhichSurf` loop up to the sort: the UV trim curve
/// of `GetEdgeuse(which_surf)` (`GetOrCreateUVTrimCurve`, `SER`); `None` when it is
/// closed on pass 0 (`continue`), `st.do_sorting = false` when closed on pass 1; the
/// split parameter (temporary `m_vG1Knots` push/pop); a point per G1 knot
/// (`binormal_guess_point(.., half = true)`, `sm_FindFilletTsectPnt`, split capture);
/// the rail-end points (`SM_FV_MATE` / `SM_FV_RAIL_X_EXTENDED_EDGEUSE` skipped,
/// `SM_OT_SAME` in front, `smos_CleanUpNoise` written back into the vertex use).
fn c1_collect_points(
    ops: &mut dyn FilletSolverOps,
    cx: &mut Cx<'_>,
    st: &mut CalcState,
    which_surf: usize,
) -> Result<Option<SurfPass>> {
    let _ = (ops, cx, st, which_surf);
    // P4-CHUNK(C1): SmFilletSolver.cpp:1587-1752
    Err(SolidError::Unported {
        upstream: "SmFilletSolver::CalcFilletGeom (p4 chunk C1 pending)",
        unit: "blend.fillet_solver_p4",
    })
}

// P4-CHUNK(C1) END

// =====================================================================================
// P4-CHUNK(C2) BEGIN: sorting and end cleanup, SmFilletSolver.cpp:1754-1920
// (chunk-local imports go inside the functions: never edit the shared `use` block)
// =====================================================================================

/// :1754-1858: `dScale`; for each of `pass.all` (the `-SM_BIG_DOUBLE` UV skip kept),
/// `GlobalPointSolve` on `pass.uv_curve` (exactly one solution, else `SM_ERR`), the
/// boundary adjustment along the tangent line, and the deduplicating insertion sort
/// into `st.tsect_pnts` that marks each inserted point's `curve_parameter` with
/// `-SM_BIG_DOUBLE` (written back through its handle).
fn c2_sort_points(
    ops: &mut dyn FilletSolverOps,
    cx: &mut Cx<'_>,
    st: &mut CalcState,
    which_surf: usize,
    pass: &SurfPass,
) -> Result<()> {
    let _ = (ops, cx, st, which_surf, pass);
    // P4-CHUNK(C2): SmFilletSolver.cpp:1754-1858
    Err(SolidError::Unported {
        upstream: "SmFilletSolver::CalcFilletGeom (p4 chunk C2 pending)",
        unit: "blend.fillet_solver_p4",
    })
}

/// :1872-1920: `SM_ERR` without points; with more than two, drop the penultimate
/// point when the last two are within `m_dThisApproxTol3d` on `GetSurface(0)` (UV of
/// surface 0, statuses ignored), then the second when the first two are.
fn c2_clean_ends(ops: &mut dyn FilletSolverOps, cx: &mut Cx<'_>, st: &mut CalcState) -> Result<()> {
    let _ = (ops, cx, st);
    // P4-CHUNK(C2): SmFilletSolver.cpp:1872-1920
    Err(SolidError::Unported {
        upstream: "SmFilletSolver::CalcFilletGeom (p4 chunk C2 pending)",
        unit: "blend.fillet_solver_p4",
    })
}

// P4-CHUNK(C2) END

// =====================================================================================
// P4-CHUNK(C3) BEGIN: analytic attempt and trace setup, SmFilletSolver.cpp:1955-2090
// (chunk-local imports go inside the functions: never edit the shared `use` block)
// =====================================================================================

/// :1955-2011: under `CanMakeAnalytic()`, the base surfaces as B-splines (one offset
/// level unwrapped, `NER`), `sm_TestAnalyticSurface` on both, then
/// `CalcCylinderFilletGeom` and `CalcTorusFilletGeom`. `true` when one succeeded
/// (`return SM_SUCCESS`); their failures fall through to the trace.
fn c3_try_analytic(
    ops: &mut dyn FilletSolverOps,
    cx: &mut Cx<'_>,
    st: &mut CalcState,
) -> Result<bool> {
    let _ = (ops, cx, st);
    // P4-CHUNK(C3): SmFilletSolver.cpp:1955-2011
    Err(SolidError::Unported {
        upstream: "SmFilletSolver::CalcFilletGeom (p4 chunk C3 pending)",
        unit: "blend.fillet_solver_p4",
    })
}

/// :2028-2034: `st.trace.domains` = the base surfaces' natural domains after
/// `SetupOffsetExtension`; returns the solver's offset surfaces as [`Surface`]s for
/// the intersector.
fn c3_trace_setup(
    ops: &mut dyn FilletSolverOps,
    cx: &mut Cx<'_>,
    st: &mut CalcState,
) -> Result<[Surface; 2]> {
    let _ = (ops, cx, st);
    // P4-CHUNK(C3): SmFilletSolver.cpp:2028-2034
    Err(SolidError::Unported {
        upstream: "SmFilletSolver::CalcFilletGeom (p4 chunk C3 pending)",
        unit: "blend.fillet_solver_p4",
    })
}

/// :2037-2052: `SmFilletIntersector sFI(*pOffsetSrf1, sDomain1, *pOffsetSrf2,
/// sDomain2, *this)`; with an extension flag set and an arc edge curve
/// (`IsArc(5, m_dThisApproxTol3d, ..)`), `SetExtensionDistance(|5 * offset 0|)`.
fn c3_make_intersector<'s>(
    ops: &mut dyn FilletSolverOps,
    cx: &mut Cx<'_>,
    st: &CalcState,
    surfaces: &'s [Surface; 2],
) -> Result<TraceIntersector<'s>> {
    let _ = (ops, cx, st, surfaces);
    // P4-CHUNK(C3): SmFilletSolver.cpp:2037-2052
    Err(SolidError::Unported {
        upstream: "SmFilletSolver::CalcFilletGeom (p4 chunk C3 pending)",
        unit: "blend.fillet_solver_p4",
    })
}

/// :2054-2090: the trace locals (tolerances, `sDir = sPVSt[1]`); a TANGENT last point
/// reverses `st.tsect_pnts`, sets `sDir = -sPVEn[1]` and `bReverseTracing`; the local
/// extension flags from the solver's.
fn c3_trace_start(
    ops: &mut dyn FilletSolverOps,
    cx: &mut Cx<'_>,
    st: &mut CalcState,
) -> Result<()> {
    let _ = (ops, cx, st);
    // P4-CHUNK(C3): SmFilletSolver.cpp:2054-2090
    Err(SolidError::Unported {
        upstream: "SmFilletSolver::CalcFilletGeom (p4 chunk C3 pending)",
        unit: "blend.fillet_solver_p4",
    })
}

// P4-CHUNK(C3) END

// =====================================================================================
// P4-CHUNK(C4) BEGIN: the trace loop, SmFilletSolver.cpp:2093-2467
// (chunk-local imports go inside the functions: never edit the shared `use` block)
// =====================================================================================

/// :2093-2467, `while (!bFoundAllPoints)`: the ten-pass cap (`SetStatus`,
/// `RecordFilletError`, `SM_ERR`); the reversed extension flags (every pass);
/// `DoPointIntersection` on a snapshot of the points (written back; a failure only
/// sets the status); missed points (`curve_parameter < 0`); the boundary branch; the
/// cusp branch (centre-line end, `GlobalPointSolve`, twenty skips with
/// `binormal_guess_point(.., half = false)`, `sm_FindFilletTsectPnt`,
/// `IsPointOnCurve`, insert in front).
fn c4_trace(
    ops: &mut dyn FilletSolverOps,
    cx: &mut Cx<'_>,
    st: &mut CalcState,
    fi: &mut TraceIntersector<'_>,
) -> Result<()> {
    let _ = (ops, cx, st, fi);
    // P4-CHUNK(C4): SmFilletSolver.cpp:2093-2467
    Err(SolidError::Unported {
        upstream: "SmFilletSolver::CalcFilletGeom (p4 chunk C4 pending)",
        unit: "blend.fillet_solver_p4",
    })
}

// P4-CHUNK(C4) END

// =====================================================================================
// P4-CHUNK(C5) BEGIN: post-processing, SmFilletSolver.cpp:2470-2671
// (chunk-local imports go inside the functions: never edit the shared `use` block)
// =====================================================================================

/// :2470-2671: with one geometry and the split flag, `TopologySplit`, the rail and UV
/// curve copies, the `DropPoint` ladder, trims, the new start vertex and its rail
/// vertex uses, and the fillet surface split at the rail parameter; with several
/// geometries, `CalcCrossSection` of their cross-section edges and
/// `CalcBlendingGeom` of the blends; after a reversed trace, the rails' 3D and UV
/// curves, the centre line and the fillet surface reversed.
fn c5_post_process(
    ops: &mut dyn FilletSolverOps,
    cx: &mut Cx<'_>,
    st: &mut CalcState,
) -> Result<()> {
    let _ = (ops, cx, st);
    // P4-CHUNK(C5): SmFilletSolver.cpp:2470-2671
    Err(SolidError::Unported {
        upstream: "SmFilletSolver::CalcFilletGeom (p4 chunk C5 pending)",
        unit: "blend.fillet_solver_p4",
    })
}

// P4-CHUNK(C5) END

// =====================================================================================
// P4-CHUNK(C6) BEGIN: ComputePointValues, SmFilletSolver.cpp:3544-3768
// (chunk-local imports go inside the functions: never edit the shared `use` block)
// =====================================================================================

/// :3544-3768, the whole body of `ComputePointValues` (see
/// [`FilletSolver::compute_point_values`]). The step is
/// `k.latch_fillet_point_step(if tangent { 1.0e-3 } else { 1.0e-5 })`.
fn c6_compute_point_values(
    ops: &mut dyn FilletSolverOps,
    k: &Kernel,
    uv_values: [Point2; 2],
    uv_domain1: &Extent2d,
    uv_domain2: &Extent2d,
    curve_trace_direction: f64,
    point: &mut TsectPnt,
    previous: Option<&TsectPnt>,
) -> Result<()> {
    let _ = (
        ops,
        k,
        uv_values,
        uv_domain1,
        uv_domain2,
        curve_trace_direction,
        point,
        previous,
    );
    // P4-CHUNK(C6): SmFilletSolver.cpp:3544-3768
    Err(SolidError::Unported {
        upstream: "SmFilletSolver::ComputePointValues (p4 chunk C6 pending)",
        unit: "blend.fillet_solver_p4",
    })
}

// P4-CHUNK(C6) END

// =====================================================================================
// STUBS: calls into units that are not ported (one helper each; chunks call these)
// =====================================================================================

/// `CalcCylinderFilletGeom(sTsectPnts)` (default `pFG = NULL`).
fn calc_cylinder_fillet_geom(
    ops: &mut dyn FilletSolverOps,
    cx: &mut Cx<'_>,
    store: &mut TsectStore,
    points: &[TsRef],
) -> Result<()> {
    let _ = (ops, cx, store, points);
    // UNPORTED(blend.fillet_solver_p2): SmFilletSolver::CalcCylinderFilletGeom
    Err(SolidError::unported(
        "SmFilletSolver::CalcCylinderFilletGeom",
        "blend.fillet_solver_p2",
    ))
}

/// `CalcTorusFilletGeom(sTsectPnts)` (default `pFG = NULL`).
fn calc_torus_fillet_geom(
    ops: &mut dyn FilletSolverOps,
    cx: &mut Cx<'_>,
    store: &mut TsectStore,
    points: &[TsRef],
) -> Result<()> {
    let _ = (ops, cx, store, points);
    // UNPORTED(blend.fillet_solver_p2): SmFilletSolver::CalcTorusFilletGeom
    Err(SolidError::unported(
        "SmFilletSolver::CalcTorusFilletGeom",
        "blend.fillet_solver_p2",
    ))
}

/// `RecordFilletError(eErrorType, sTsectPnts, pMessage)` (upstream ignores its
/// status: call through `ignore_status`, which still propagates `Unported`).
fn record_fillet_error(
    ops: &mut dyn FilletSolverOps,
    cx: &mut Cx<'_>,
    store: &TsectStore,
    points: &[TsRef],
    error_type: FilletErrorType,
    message: &'static str,
) -> Result<()> {
    let _ = (ops, cx, store, points, error_type, message);
    // UNPORTED(blend.fillet_solver_p2): SmFilletSolver::RecordFilletError
    Err(SolidError::unported(
        "SmFilletSolver::RecordFilletError",
        "blend.fillet_solver_p2",
    ))
}

/// `sm_FindFilletTsectPnt(this, pBaseSrf1, pBaseSrf2, rPnt, rTangent, rTsectPnt)`:
/// the fillet point on the plane through `point` normal to `tangent`, written to
/// `out`.
fn find_fillet_tsect_pnt(
    ops: &mut dyn FilletSolverOps,
    k: &Kernel,
    base1: &Surface,
    base2: &Surface,
    point: Point3,
    tangent: Vec3,
    out: &mut TsectPnt,
) -> Result<()> {
    let _ = (ops, k, base1, base2, point, tangent, out);
    // UNPORTED(blend.fillet_solver_p3): sm_FindFilletTsectPnt
    Err(SolidError::unported(
        "sm_FindFilletTsectPnt",
        "blend.fillet_solver_p3",
    ))
}

/// `StepSolve(rTsectPnt, dStepSize, crUVDomain1, crUVDomain2, rbFoundAnswer,
/// rbClipped, rbBoundaryHit, rNextTSP)`: the flags and the next point (only its UV
/// positions set; it starts as `SmTsectPnt()`). Changes the offset surfaces' offsets.
fn step_solve(
    ops: &mut dyn FilletSolverOps,
    k: &Kernel,
    point: &TsectPnt,
    step_size: f64,
    uv_domain1: &Extent2d,
    uv_domain2: &Extent2d,
) -> Result<(StepSolveFlags, TsectPnt)> {
    let _ = (ops, k, point, step_size, uv_domain1, uv_domain2);
    // UNPORTED(blend.fillet_solver_p3): SmFilletSolver::StepSolve
    Err(SolidError::unported(
        "SmFilletSolver::StepSolve",
        "blend.fillet_solver_p3",
    ))
}

/// `pFilletGeom->TopologySplit(rpNewFilletGeom, eType)`: the new geometry.
fn topology_split(
    cx: &mut Cx<'_>,
    geom: FilletGeomId,
    geom_type: FilletGeomType,
) -> Result<FilletGeomId> {
    let _ = (cx, geom, geom_type);
    // UNPORTED(blend.fillet_geom_p3): SmFilletGeom::TopologySplit
    Err(SolidError::unported(
        "SmFilletGeom::TopologySplit",
        "blend.fillet_geom_p3",
    ))
}

/// `pFilletGeom->CalcBlendingGeom()`.
fn calc_blending_geom(cx: &mut Cx<'_>, geom: FilletGeomId) -> Result<()> {
    let _ = (cx, geom);
    // UNPORTED(blend.fillet_geom_p4): SmFilletGeom::CalcBlendingGeom
    Err(SolidError::unported(
        "SmFilletGeom::CalcBlendingGeom",
        "blend.fillet_geom_p4",
    ))
}

/// `pFilletEdge->CalcCrossSection(m_crContext)` (a pseudo-B-rep edge).
fn calc_cross_section(cx: &mut Cx<'_>, edge: EdgeId) -> Result<()> {
    let _ = (cx, edge);
    // UNPORTED(blend.fillet_geom_p5): SmFilletEdge::CalcCrossSection
    Err(SolidError::unported(
        "SmFilletEdge::CalcCrossSection",
        "blend.fillet_geom_p5",
    ))
}

/// `SmEdgeuse::GetOrCreateUVTrimCurve(rpUVTrimCurve)` (defaults: no distances,
/// `bUVSpaceOkay = TRUE`): a copy of the trim curve of the edge use or its mate
/// (`GetUVTrimCurve`), failing with the edge use's drop-curve failure status when it
/// has one; building a missing curve (`CreateUVTrimCurve`, which would attach it to
/// `brep`) is unported.
fn get_or_create_uv_trim_curve(brep: &mut Brep, edgeuse: EdgeuseId) -> Result<Curve> {
    let (curve, _) = eu::uv_trim_curve(brep, edgeuse)?;
    let Some(curve) = curve.cloned() else {
        // UNPORTED(topology.edgeuse_p3): SmEdgeuse::GetOrCreateUVTrimCurve
        return Err(SolidError::unported(
            "SmEdgeuse::GetOrCreateUVTrimCurve",
            "topology.edgeuse_p3",
        ));
    };
    let (fail, _) = eu::drop_curve_fail(brep, edgeuse)?;
    if let Some(fail) = fail {
        if fail.fail_status != SM_SUCCESS {
            return Err(SolidError::Upstream {
                operation: OP,
                code: fail.fail_status,
                entity: None,
            });
        }
    }
    Ok(curve)
}

/// The outputs of `SmEdgeuse::EvaluateBinormal(t, FALSE, rPnt, rBinVec, &tangent,
/// NULL, &uv)`.
#[derive(Debug, Clone, Copy)]
struct EdgeuseBinormal {
    /// `rBinormalPoint`: the 3D point on the edge.
    point: Point3,
    /// `rBinormalVec`: unit vector into the face.
    binormal: Vec3,
    /// `pOptTangent`: the (non-unit) edge-use tangent.
    tangent: Vec3,
    /// `pOptUVParameter`: the surface parameter of the point.
    uv: Point2,
}

/// `SmEdgeuse::EvaluateBinormal(t, FALSE, ..)` (`FALSE`: get or create the UV trim
/// curve, so the target B-rep may change).
fn evaluate_binormal(
    k: &Kernel,
    brep: &mut Brep,
    edgeuse: EdgeuseId,
    t: f64,
) -> Result<EdgeuseBinormal> {
    let _ = (k, brep, edgeuse, t);
    // UNPORTED(topology.edgeuse_p3): SmEdgeuse::EvaluateBinormal
    Err(SolidError::unported(
        "SmEdgeuse::EvaluateBinormal",
        "topology.edgeuse_p3",
    ))
}

/// The inputs of `sFI.DoPointIntersection(m_crContext, sTsectPnts, bExtendBefore,
/// bExtendAfter, &sDir, NULL, &dThisApproxTol3d, &dAngleTol, ..)`.
#[derive(Debug, Clone, Copy)]
struct PointIntersectionArgs {
    /// `bExtendBefore`.
    extend_before: bool,
    /// `bExtendAfter`.
    extend_after: bool,
    /// `&sDir` (the end direction is `NULL`).
    start_dir: Option<Vec3>,
    /// `&dThisApproxTol3d`.
    approx_tol_3d: f64,
    /// `&dAngleTol`.
    ang_tol_rad: f64,
}

/// `sFI.DoPointIntersection(..)` over `points` (a snapshot of `sTsectPnts`, written
/// in place as upstream writes through the pointers; the curve outputs are unused by
/// `CalcFilletGeom`; `FlushCurve` fills the solver's last fillet geometry).
fn do_point_intersection(
    fi: &mut TraceIntersector<'_>,
    ops: &mut dyn FilletSolverOps,
    cx: &mut Cx<'_>,
    points: &mut Vec<TsectPnt>,
    args: PointIntersectionArgs,
) -> Result<()> {
    let _ = (fi, ops, cx, points, args);
    // DISPATCH GAP (intersect.surface_intersector_p4): SmSurfaceIntersector::DoPointIntersection
    // runs the basic intersector's ComputeNextPoint / FlushCurve / ComputePointValues /
    // TestSpanAccuracy; the SmFilletIntersector overrides need a dispatch hook there.
    Err(SolidError::unported(
        "SmSurfaceIntersector::DoPointIntersection (SmFilletIntersector dispatch)",
        "intersect.surface_intersector_p4",
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resolver_vec_maps_and_registers() {
        let mut map = vec![FilletGeomId(7)];
        assert_eq!(
            map.geometry_id(FilletGeomRef { index: 0 }).unwrap(),
            FilletGeomId(7)
        );
        assert!(map.geometry_id(FilletGeomRef { index: 1 }).is_err());
        map.register_geometry(FilletGeomRef { index: 1 }, FilletGeomId(9))
            .unwrap();
        assert_eq!(
            map.geometry_id(FilletGeomRef { index: 1 }).unwrap(),
            FilletGeomId(9)
        );
        assert!(map
            .register_geometry(FilletGeomRef { index: 5 }, FilletGeomId(1))
            .is_err());
    }

    #[test]
    fn store_handles_write_back_to_their_storage() {
        let mut fillets = FilletBrep::new();
        let mut store = TsectStore::default();
        let a = store.alloc(TsectPnt::new());
        let b = store.alloc(TsectPnt::new());
        assert_eq!((a, b), (TsRef::Local(0), TsRef::Local(1)));
        let mut snap = store.snapshot(&fillets, &[b, a]).unwrap();
        snap[0].curve_parameter = 2.0;
        snap[1].curve_parameter = 1.0;
        store.write_back(&mut fillets, &[b, a], &snap).unwrap();
        assert_eq!(store.get(&fillets, a).unwrap().curve_parameter, 1.0);
        assert_eq!(store.get(&fillets, b).unwrap().curve_parameter, 2.0);
        assert!(store.get(&fillets, TsRef::Local(2)).is_err());
        assert!(store.write_back(&mut fillets, &[a], &snap).is_err());
    }

    #[test]
    fn placeholders_report_the_pending_chunk() {
        let k = Kernel::default();
        let mut solver = FilletSolver::new(1.0e-5, 0.01, 0.01);
        let mut p = TsectPnt::new();
        let d = Extent2d::from_bounds(0.0, 0.0, 1.0, 1.0);
        let e = FilletSolver::compute_point_values(
            &mut solver,
            &k,
            [Point2::new(0.0, 0.0); 2],
            &d,
            &d,
            1.0,
            &mut p,
            None,
        )
        .unwrap_err();
        assert!(e.is_unported());
    }
}
