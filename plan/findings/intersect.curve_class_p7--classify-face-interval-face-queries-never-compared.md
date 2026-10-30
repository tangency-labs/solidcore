# intersect.curve_class_p7: ClassifyFaceInterval's face queries are never compared against C++ (2026-10-30)

Unit: `intersect.curve_class_p7` (runner `crates/difftest/src/units/intersect_curve_class_p7.rs`,
shim `harness/cpp/shim/intersect_curve_class_p7.cpp`).

## Gap

The shim's `OP_CLASSIFY_FACE` returns `SC_ERR_ARG` for every `ClassifyFaceInterval` that
gets past its early returns (middle unknown, or a face one given an old face) and has a
vertex or edge end, or asks for the point classification:

```cpp
if (!early && (is_ve(ivl.m_vStart) || is_ve(ivl.m_vEnd) || do_point))
    return SC_ERR_ARG;
```

The guard exists because the bare stand-in faces have no loops for the face queries to
walk, but it also applies when the face is the real one-face B-rep (`has_face`, whose
topology the pool marks `real`). The runner's `generate` drops every script the shim
rejects (`shim_accepts`), so no generated or edge case in this unit runs:

- `SmFace::FindVertexSector` (vertex end; Rust `Mutable` / `ReadOnly::find_vertex_sector`
  in `crates/solidcore/src/intersect/curve_class/p7.rs`),
- `SmFace::EdgeCurveOnClassify` (edge end),
- `SmFace::PointClassify` / `Point3DClassify` (point classification),

so none of them, on the mutable or the read-only path, has been compared against C++ by
this unit. The unit's green status covers only the early-return paths of
`ClassifyFaceInterval`. (Library callers reach the mutable path through booleans and
fillets end to end; whether those units' fixtures exercise these queries was not measured.)

## How it was found

Wiring the read-only `SmFace::FindVertexSector` (dropped, stub kept), read-only twin
cases with a vertex start on the real B-rep were rejected by this guard, leaving 0 cases
(1516 cases before and after, `--filter read_only` 0, 0 breakpoint hits in the Rust
sector search).

## Available follow-up

Lift the guard when the face is the real B-rep's face and every vertex or edge end is one
of that B-rep's own entities (`pool.is_real`), with the point classification too once its
loop walk on the real face is checked. This changes which generated scripts the shim
accepts (`generate` retries on rejection), so the generated case set changes; any
mismatches it exposes stay red with findings. With it in place, the read-only
`FindVertexSector` stub (`UNPORTED(topology.face_p4)`, `ReadOnly::find_vertex_sector`) can
be wired against a reaching case.
