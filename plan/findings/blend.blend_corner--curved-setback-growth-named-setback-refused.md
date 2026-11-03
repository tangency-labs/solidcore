# blend.blend_corner: the setback a curved Setback's growth names as building is refused when asked for (2026-11-03)

- **Unit:** `blend.blend_corner` (the Setback beside curved faces, #2343, PR #2366; `Kernel::fixes.fillet_blend_corner`).
- **Inputs:** the D prism (R 10, flat a 4, H 10), its top corner (4, sqrt 84, 10), the chord, the upright and the arc rounded 2, 2.5 and 3 (`crates/solidcore/tests/issue_2343_curved_setback.rs`, `d_prism_case`).
  - **Asked 3.05:** the corner is refused: "a setback of 3.05 is too short for the setback corner's fill to stay no more concave than its faces: a setback of 3.24062 builds it". `fit_setback` is 3.2406249999999996 (3.05 x 1.0625).
    - The growth's trials at 2, 1.5, 1.25, 1.125 and 1.0625 certify; those at 1.03125, 1.046875 and 1.0546875 do not.
    - At 1.0625 the fill certifies only at degree 9: degrees 4 to 7 are infeasible, and degree 8 solves but is not certified.
  - **Asked 3.2406249999999996:** refused in turn, naming 3.26594.
- **Expected:** a setback the corner names as building builds when asked for (`fit_setback`'s contract).
- **Suspected cause:** the growth's trials are not independent of their order. Each trial places the setback vertices again (`CalcSetBackVert`), starting from the previous trial's state. That state is the vertex use's `TsectPnt` and whatever the solver keeps between solves (`point_on_plane_solve` writes its rail offsets back). A fresh build at the same setbacks starts from the trace's state. At 1.0625 the certificate is marginal (degree 9 only), so a placement that differs by rounding certifies in the trial and not in the fresh build.
  - The fork's work in progress notes the same class open at 37c9e1d: "the cone still differs grown vs scaled with identical round rows".
- **Not a hang:** the growth is bounded: at most nine fills, seven degrees each, at most 200 row rounds per quadratic programme, each with at most `50 (m + n) + 1000` dual steps. Asked 3.05 takes 136 s in release, and most of an hour in a debug build. Each solve at degree 9 near the threshold takes 8 to 27 s in release, mostly the dual active-set steps over a working set of up to 12000 rows. A debug-build probe of it was first taken for a hang.
- **State:** open. The test `d_prism_asked_too_short_names_the_setback_that_builds` checks the refusal and its named setback. `the_named_setback_builds` is ignored as known failing, with this file as its reason. Both are ignored as slow.
