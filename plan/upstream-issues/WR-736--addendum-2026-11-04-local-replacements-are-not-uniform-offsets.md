Updates: WR-736

The vertex-image construction applies to uniform normal offsets. A local operation
supplies independent replacement surfaces; its executive's `offset_distance` is
half the bounding-sphere radius, used to size support extensions. It is not the
replacement displacement. Applying `V + d n(V)` to local lamina edges therefore
places corners at unrelated points and can fail gap closing with status 1001.

`trim_extended_edge` now dispatches on the existing operation kind: uniform offsets
retain WR-736's exact vertex-image construction, while local operations use the
replacement surfaces' trim constraints. No new switch, retry, tolerance change or
validation exception is introduced.

Evidence: both existing `issue_2048_check_offset_surface_face` eight-turn sheet
regressions fail in `step_gaps` with all fixes enabled before this correction.
A standalone probe linked to that same library passes both tests (including both
signed support-continuation cases) when only WR-736's existing flag is disabled.
The committed regressions explicitly keep that flag enabled and check exact
surface evaluations, topology preservation and strict result validity. Existing
`issue_2611_thicken_arc_trim`, `issue_2174_cap_rim_two_edges` and
`issue_1502_z5_cap_thicken` cover the uniform-offset construction which remains
unchanged. This is a Rust integration defect; no C++ reproduction is claimed.
