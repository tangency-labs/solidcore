Updates: WR-599

#### WR-599 (update, 2026-11-03). Generic topological cleanup completes the merged face's seam

- **Location:** `topology::brep::p5::remove_topological_with_notifications`, after
  `DeleteEdge` and the merged UV domain update. `boolean::merge::p1` already completes
  the seam at the corresponding point of boolean postprocessing.
- **What happens:** `cover_wraps_across_closed_seam` allows both cleanup paths to merge
  faces across a periodic surface's seam. Only boolean postprocessing reconstructed
  the missing seam and split boundary edges crossing it. Generic cleanup could return
  a manifold body with invalid UV trims. SolidCraft's pre-offset healing reaches this
  path before the offset itself (SolidCraft #11582).
- **Reproduction:** construct the existing WR-599 peg with wrapped coverage disabled,
  then enable it and call `remove_topological_edges_and_vertices`. The input passes
  level-2 validation; cleanup makes the joined cylinder wall fail with a UV endpoint
  gap of one period, `50.265482` at radius 8. No screw sweep is involved.
- **Construction:** both merge paths call the same face-layer helper after setting
  the merged domain. It classifies the exact trim loops and calls `SplitAtSeam` when
  a seam is crossed or missing, under the existing option. The caller's notification
  environment is passed through, so boundary splits retain their provenance. There
  is no change to cover acceptance or the option-disabled path.
- **Regression:** the existing
  `crates/solidcore/tests/fork_16bf1cb_full_turn_cover_across_seam.rs` checks generic
  cleanup's validity, seven faces, the new seam endpoint, closed-form volume,
  repeated-cleanup stability and the boundary split's caller-owned notifications.
- **Scope:** a Rust caller-completion defect demonstrated locally. No new C++ parity
  claim; the existing boolean seam construction is shared without numerical changes.
