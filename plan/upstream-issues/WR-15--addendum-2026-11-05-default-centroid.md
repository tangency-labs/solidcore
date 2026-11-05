Updates: WR-15

#### WR-15 (update, 2026-11-05). All centroid blocks by default

`PointSet3d::centroid` now uses the former all-blocks path and caches its result.
For 100 copies of `(2, 4, 6)` queried before axes, the former default returned
`(1, 2, 3)`; the mean is `(2, 4, 6)`. Empty sets still return zero. The existing
weighted block-sum arithmetic and its IEEE overflow behavior are unchanged.
Axes computed without a supplied normal can replace the cached centroid with the
best-fit plane centroid; arbitrary floating-point inputs are not promised bitwise
call-order identity.

The flag `point_set_centroid_all_blocks` and methods `centroid_opt_in` and
`centroid_with_kernel` are removed; callers should use `centroid()`. This is a
source-breaking API change. The historical C++ finding remains confirmed; any
divergence on current native fixtures has not been verified in this change.
