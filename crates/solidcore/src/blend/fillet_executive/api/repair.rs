//! Opt-in repair of displaced split vertices in completed fillets (#1999).
//!
//! Upstream: `source/SM_API/src/SmApiBrep.cpp` (`SmApiFilletEdges`);
//! this validation and repair have no upstream equivalent.
//!
//! Deviations: only a degree-two topological vertex whose two curve endpoints agree
//! within scaled zero may move, and only within both existing vertex/edge intersection
//! tolerances. Curves, parameters, topology and tolerances are unchanged. The public
//! wrapper rolls back if the whole result still fails Level2. No recursion.

use crate::error::{EntityRef, Result, SolidError};
use crate::geom::tol::{scaled_zero2, xsect_tol_3d_from_zones};
use crate::kernel::Kernel;
use crate::math::Point3;
use crate::topology::brep::p3::with_kernel_callbacks;
use crate::topology::brep::validation::{check_brep, CheckLevel};
use crate::topology::{edge, vertex, Brep, EdgeId, VertexId};

fn valid(k: &mut Kernel, b: &Brep) -> Result<bool> {
    // AssertValid can update gap/tolerance metadata and caches. Keep those
    // validation side effects out of the accepted operation's actual result.
    let mut probe = b.clone_for_rollback()?;
    Ok(check_brep(k, &mut probe, CheckLevel::Level2)?.is_valid())
}

/// The end of a non-closed edge at `v`; a missing curve cannot support a repair.
fn endpoint(k: &Kernel, b: &Brep, e: EdgeId, v: VertexId) -> Result<Option<Point3>> {
    let (start, end) = edge::start_end_vertices(b, e)?;
    let iv = edge::interval(b, e)?;
    let t = match (start == Some(v), end == Some(v)) {
        (true, false) => iv.min(),
        (false, true) => iv.max(),
        _ => return Ok(None),
    };
    edge::curve(b, e)?
        .map(|c| c.evaluate_point(k, t))
        .transpose()
}

fn replacement(k: &Kernel, b: &mut Brep, v: VertexId) -> Result<Option<Point3>> {
    let edges = vertex::edges(b, v, None)?;
    let [first, second] = edges.as_slice() else {
        return Ok(None);
    };
    if first == second {
        return Ok(None);
    }
    let (Some(p), Some(q)) = (endpoint(k, b, *first, v)?, endpoint(k, b, *second, v)?) else {
        return Ok(None);
    };
    if ![p.x, p.y, p.z, q.x, q.y, q.z]
        .into_iter()
        .all(f64::is_finite)
        || p.distance(q) > scaled_zero2(&p, &q)
    {
        return Ok(None);
    }
    let gap = vertex::point(b, v)?.distance(p);
    let vertex_tol = vertex::tolerance(b, v)?;
    let limit = xsect_tol_3d_from_zones(k, vertex_tol, edge::tolerance(b, *first)?).min(
        xsect_tol_3d_from_zones(k, vertex_tol, edge::tolerance(b, *second)?),
    );
    if !gap.is_finite()
        || !limit.is_finite()
        || limit <= 0.0
        || gap <= vertex::corner::xsect_tol_3d_of_brep(k, b)
        || gap > limit
        || !vertex::is_topological_vertex(k, b, v)?
    {
        return Ok(None);
    }
    // Choose the first curve endpoint in native adjacency order, rather than
    // averaging points off both curves. The second agrees within scaled zero.
    Ok(Some(p))
}

pub(super) fn repair_result(k: &mut Kernel, b: &mut Brep) -> Result<()> {
    if valid(k, b)? {
        return Ok(());
    }
    let vertices: Vec<_> = b.vertices().ids().collect();
    let mut repairs = Vec::new();
    for v in vertices {
        if let Some(p) = replacement(k, b, v)? {
            repairs.push((v, p));
        }
    }
    if !repairs.is_empty() {
        with_kernel_callbacks(k, b, |_k, b, env| {
            for (v, p) in repairs {
                let this = (env.key_of)(EntityRef::Vertex(v));
                vertex::set_point(
                    b,
                    v,
                    p,
                    env.records,
                    env.attributes,
                    env.context.as_deref_mut(),
                    this,
                    &|vu| (env.key_of)(EntityRef::Vertexuse(vu)),
                )?;
            }
            Ok(())
        })?;
        b.refresh_cached_topology();
        if valid(k, b)? {
            return Ok(());
        }
    }
    Err(SolidError::sm_err(
        "blend::fillet_executive::api::invalid_result",
    ))
}
