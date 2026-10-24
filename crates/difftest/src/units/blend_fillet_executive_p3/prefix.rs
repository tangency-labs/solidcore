//! Actual shared parameter setup, corner creation and direct first-Nx1 topology.
//! The direct topology checkpoint is not a full `CreateFilletBrep` prelude: neither
//! side runs smooth-sequence processing. With the one actual solver in these
//! fixtures, its native `while (count > 1)` loop performs no work.
use super::{label, section, status, take_section, Checkpoint, RealSolid};
use crate::{HarnessError, Output};
use solidcore::blend::fillet_executive::session::{FilletObject, FilletSession, SolverObject};
use solidcore::kernel::Kernel;
use solidcore::topology::brep::p1;
use solidcore::topology::{edge, edgeuse::p1 as eu, Brep, EdgeId};
use solidcore::{EntityRef, Result as SolidResult, SolidError};
use solidcore_ffi::units::blend_fillet_executive_p3 as ffi;

fn invalid(reason: &'static str) -> SolidError {
    SolidError::invalid_input("harness::fillet_prefix", reason)
}

fn snapshot(session: &FilletSession<'_>, parameters: &mut Vec<f64>) -> SolidResult<Vec<f64>> {
    let exec = session.executive();
    let target = session.target();
    let target_vertices = p1::vertices(target, None)?;
    let target_edges = p1::edges(target, None)?;
    let target_faces = p1::faces(target, None)?;
    let mut target_uses = Vec::new();
    for edge in target_edges {
        target_uses.extend(edge::edgeuses(target, edge, None, None)?);
    }
    let pseudo = session.pseudo_brep();
    let vertices = p1::vertices(pseudo.brep(), None)?;
    let edges = pseudo.edge_list();
    let mut corners = Vec::new();
    let mut corner_ids = Vec::new();
    for object in &exec.corners {
        let FilletObject::Corner(id) = object else {
            return Err(invalid("actual corner required"));
        };
        corner_ids.push(*id);
        corners.push(session.corner(*id)?);
    }
    let corner_vertices: Vec<_> = corners
        .iter()
        .flat_map(|c| c.vertices.iter().copied())
        .collect();
    let mut geometry_ids = Vec::new();
    for object in &exec.solvers {
        let FilletObject::Solver(id) = object else {
            return Err(invalid("actual solver required"));
        };
        geometry_ids.extend(session.solver(*id)?.base().fillet_geoms.iter().copied());
    }
    let mut state = vec![
        f64::from(exec.options.global_merge),
        exec.solvers.len() as f64,
        corners.len() as f64,
        vertices.len() as f64,
        edges.len() as f64,
    ];
    for corner in corners {
        parameters.extend([corner.approx_tol_3d, corner.tangency_tol_radians]);
        state.extend([
            f64::from(corner.corner_type() as u32),
            label(&target_vertices, corner.vertex),
            corner.vertices.len() as f64,
            corner.edges.len() as f64,
        ]);
        for &id in &corner.vertices {
            let vertex = pseudo.fillet_vertex(id)?;
            let object = vertex.point_class_object();
            let object_label = match object {
                None => -1.0,
                Some(EntityRef::Edgeuse(eu)) => label(&target_uses, Some(eu)),
                Some(_) => -2.0,
            };
            state.extend([
                label(&corner_vertices, Some(id)),
                label(&vertices, Some(id)),
                f64::from(vertex.fillet_vertex_type() as u32),
                label(&corner_ids, vertex.fillet_corner()),
                f64::from(vertex.point_class() as u32),
                object.map_or(0.0, |object| {
                    session.target_identity(object).type_code as f64
                }),
                object_label,
            ]);
        }
        for &edge in &corner.edges {
            state.push(label(edges, Some(edge)));
        }
    }
    for object in &exec.solvers {
        let FilletObject::Solver(id) = object else {
            return Err(invalid("actual solver required"));
        };
        let solver = session.solver(*id)?;
        let constant = if let SolverObject::ConstantRadius(solver) = solver {
            Some(solver)
        } else {
            None
        };
        let solver_type = match solver {
            SolverObject::ConstantRadius(s) => s.constant_radius_solver_type(),
            SolverObject::ConstantDistance(s) => s.solver_type(),
            SolverObject::VariableRadius(s) => s.solver_type(),
            SolverObject::SurfaceSurface(s) => s.solver_type(),
        };
        state.extend([
            f64::from(solver_type as u32),
            if constant.is_some() { 2.0 } else { 0.0 },
        ]);
        if let Some(constant) = constant {
            parameters.extend([constant.fillet_radius(0)?, constant.fillet_radius(1)?]);
        }
        let base = solver.base();
        parameters.extend([
            base.this_approx_tol_3d(),
            base.this_ang_tol_rad(),
            base.tangency_tolerance(),
        ]);
        state.push(base.fillet_geoms.len() as f64);
        for &geometry in &base.fillet_geoms {
            for rail in 0..2 {
                state.push(label(edges, pseudo.geom(geometry)?.rail(rail)));
            }
        }
    }
    for &id in edges {
        let edge = pseudo.fillet_edge(id)?;
        state.extend([
            f64::from(edge.fillet_edge_type() as u32),
            label(&target_faces, edge.original_face()),
            label(&corner_ids, edge.fillet_corner()),
        ]);
        let uses = edge::edgeuses(pseudo.brep(), id, None, None)?;
        state.push(uses.len() as f64);
        for use_ in uses {
            let edgeuse = pseudo.brep().edgeuse(use_)?;
            let vertex = edgeuse
                .vertexuse
                .map(|vu| pseudo.brep().vertexuse(vu).map(|vu| vu.vertex))
                .transpose()?
                .flatten();
            state.extend([
                f64::from(eu::orientation(pseudo.brep(), use_)? as i32),
                label(&corner_vertices, vertex),
                label(&vertices, vertex),
                label(&geometry_ids, pseudo.fillet_edgeuse(use_)?.fillet_geom()),
            ]);
        }
    }
    Ok(state)
}

pub(super) fn run(
    k: &mut Kernel,
    body: &mut Brep,
    selected: &[EdgeId],
    input: &RealSolid,
    before: &[f64],
) -> SolidResult<Vec<f64>> {
    let mut session = FilletSession::new(k, body)?;
    let mut statuses = [-1.0; 3];
    let prepare = session.prepare_constant_radius_edges(
        k,
        selected,
        input.radius,
        input.xsect,
        input.continuity,
        input.thumbweight,
    );
    statuses[0] = match prepare {
        Ok(()) => 0.0,
        Err(error) => status(error)?,
    };
    let mut prior = Vec::new();
    if statuses[0] == 0.0 && input.checkpoint != Checkpoint::Prepared {
        statuses[1] = match session.create_corners(k, None) {
            Ok(()) => 0.0,
            Err(error) => status(error)?,
        };
        if statuses[1] == 0.0 && input.checkpoint == Checkpoint::FirstTopology {
            prior = snapshot(&session, &mut Vec::new())?;
            if let Some(FilletObject::Corner(id)) = session.executive().corners.first().copied() {
                statuses[2] = match session.make_nx1_corner_topology(k, id) {
                    Ok(()) => 0.0,
                    Err(error) => status(error)?,
                };
            }
        }
    }
    let mut parameters = Vec::new();
    let state = snapshot(&session, &mut parameters)?;
    let mut out = Vec::new();
    section(&mut out, before);
    out.extend(statuses);
    section(&mut out, &prior);
    section(&mut out, &state);
    section(&mut out, &parameters);
    Ok(out)
}

pub(super) fn output(mut values: &[f64]) -> Result<Output, HarnessError> {
    let mut out = Output::new();
    out.values.insert(
        "selected_edge_before".into(),
        take_section(&mut values)?.to_vec(),
    );
    let statuses = values
        .get(..3)
        .ok_or_else(|| HarnessError::new("missing prefix statuses"))?;
    if statuses.iter().any(|v| !v.is_finite() || v.fract() != 0.0) {
        return Err(HarnessError::new("invalid prefix statuses"));
    }
    out.ints.insert(
        "prefix_statuses".into(),
        statuses.iter().map(|v| *v as i64).collect(),
    );
    values = values
        .get(3..)
        .ok_or_else(|| HarnessError::new("missing prefix state"))?;
    for key in ["before_direct_topology", "prefix_state"] {
        let state = take_section(&mut values)?;
        if state.iter().any(|v| !v.is_finite() || v.fract() != 0.0) {
            return Err(HarnessError::new("invalid prefix state"));
        }
        out.ints
            .insert(key.into(), state.iter().map(|v| *v as i64).collect());
    }
    out.values.insert(
        "actual_solver_corner_parameters".into(),
        take_section(&mut values)?.to_vec(),
    );
    if !values.is_empty() {
        return Err(HarnessError::new("trailing prefix data"));
    }
    Ok(out)
}
pub(super) fn cpp(input: &RealSolid) -> Result<Output, HarnessError> {
    let values = ffi::real_prefix(&input.wire(), input.checkpoint.code())
        .map_err(|e| HarnessError::new(e.to_string()))?;
    output(&values)
}
