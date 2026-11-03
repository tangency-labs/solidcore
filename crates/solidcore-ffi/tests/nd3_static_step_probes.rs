//! ND-3: on a process that latched nothing, the two `sdStaticStep` probes read unlatched,
//! so `check_parent_static_step` passes a clean parent. Against the SMLib fork the difftest
//! is pinned to, which picks the step for each point (no `sdStaticStep`, tangency-labs/smlib
//! 9e17f64), both read `1e-5`: the crossing point's own step, whatever ran before it. (v1.0.1
//! read `1e-3` for the second, latched by the tangent point.) Its own test binary: no other
//! test runs C++ in this process first.

#![cfg(feature = "cpp")]

use solidcore_ffi::isolate::fresh_process;
use solidcore_ffi::units::blend_fillet_solver_p4::{
    check_parent_static_step, static_step, static_step_after_tangent,
};

#[test]
fn a_clean_process_reads_both_probes_unlatched() {
    let plain = fresh_process(|| Ok(vec![static_step()?])).expect("plain probe");
    let tangent = fresh_process(|| Ok(vec![static_step_after_tangent()?])).expect("tangent probe");
    eprintln!("ND-3 probes on a clean process: plain {plain:?}, after a tangent point {tangent:?}");
    assert_eq!(plain, [1.0e-5], "the crossing point's own step");
    assert_eq!(tangent, [1.0e-5], "its own step after a tangent point too");
    check_parent_static_step("in the self-check").expect("a clean parent reads unlatched");
}
