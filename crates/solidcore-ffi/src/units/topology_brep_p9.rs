//! Native `topology.brep_p9` entrypoints over validated real B-rep fixtures.
//!
//! The export call builds exactly the same initial graph as the operation call.
//! For insertion it adds two isolated interior loop vertices as fixture setup.
//! Outputs are status, record length, traversal-index record, raw graph export.
use super::topology_vertexuse::{self as fixture, Op};
use crate::{call_with_buffer, FfiError};
use std::ffi::CString;

mod sys {
    extern "C" {
        pub fn sc_topology_brep_p9_eval(
            op: i32,
            export_only: i32,
            selector: i32,
            variant: i32,
            editing: i32,
            same: i32,
            path: *const std::ffi::c_char,
            input: *const f64,
            n: i64,
            parameter: f64,
            point: *const f64,
            start_tol: f64,
            end_tol: f64,
            out: *mut f64,
            cap: i64,
            len: *mut i64,
        ) -> i32;
    }
}

/// One native operation and its fixture construction program.
#[derive(Debug, Clone, Copy)]
pub struct Request<'a> {
    /// 0 infinite region, 1 orient sheet, 2 shell vertex, 3 vertex loop,
    /// 4 split edge, 5 insert edge in face, 6/7 point-class merge into edge/vertex.
    pub op: u8,
    /// Entity selector, resolved modulo the actual count.
    pub selector: u8,
    /// Insertion refusal/default-region selector (0..=5).
    pub variant: u8,
    /// Caller editing state.
    pub editing: bool,
    /// Same/opposite orientation.
    pub same: bool,
    /// Upstream file path, empty for procedural fixtures.
    pub path: &'a str,
    /// Existing `topology_vertexuse` fixture program.
    pub block: &'a [f64],
    /// Split parameter in the curve's actual domain.
    pub parameter: f64,
    /// Point for vertex creation.
    pub point: [f64; 3],
    /// Raw endpoint tolerances, converted by native `SmZoneTol3d` before assignment.
    pub tolerances: [f64; 2],
}

/// Export the initial graph (with its prepared VU indices) or run the operation.
/// Refuses only invalid layouts, missing required fixture entities and nonfinite
/// scalars before native dereferences. Safe NULL refusal branches remain callable.
pub fn eval(request: Request<'_>, export_only: bool) -> Result<Vec<f64>, FfiError> {
    let bad = || FfiError::InvalidArgument("brep_p9 request or fixture".into());
    if request.op > 7
        || request.variant > 5
        || !request.parameter.is_finite()
        || request.point.iter().any(|v| !v.is_finite())
        || request
            .tolerances
            .iter()
            .any(|v| !v.is_finite())
    {
        return Err(bad());
    }
    let raw = fixture::eval(Op::Export, request.path, request.block)?;
    let counts = fixture::counts(&raw).ok_or_else(bad)?;
    let required = match request.op {
        2 => counts[0] > 0,
        3 => counts[3] > 0,
        4 | 6 | 7 => counts[7] > 0,
        5 if request.variant == 0 && request.editing => counts[3] > 0,
        _ => true,
    };
    if !required {
        return Err(bad());
    }
    let path = CString::new(request.path).map_err(|_| bad())?;
    call_with_buffer(4096, |out, len| {
        // SAFETY: fixture::eval validates the entire construction program and
        // path; scalar ranges and required entities are checked above. Native
        // additionally checks the selected edge's curve and endpoints. Buffers
        // and C strings stay live; all native allocations stay inside the shim.
        unsafe {
            sys::sc_topology_brep_p9_eval(
                i32::from(request.op),
                i32::from(export_only),
                i32::from(request.selector),
                i32::from(request.variant),
                i32::from(request.editing),
                i32::from(request.same),
                path.as_ptr(),
                request.block.as_ptr(),
                request.block.len() as i64,
                request.parameter,
                request.point.as_ptr(),
                request.tolerances[0],
                request.tolerances[1],
                out.as_mut_ptr(),
                out.len() as i64,
                len,
            )
        }
    })
}
