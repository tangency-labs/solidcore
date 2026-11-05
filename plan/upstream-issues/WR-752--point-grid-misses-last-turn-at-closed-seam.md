#### WR-752. `SmPointGrid::SetUpMonotonic` omits the last turn of a closed row or column

- **Severity:** Wrong result
- **Area:** Surfaces and volumes
- **Location:** `source/SMLib/src/SmPointSet.cpp:2179`, `SmPointGrid::SetUpMonotonic`; row walk `:2210-2281`, column walk `:2285-2336` at v1.0.1 (`3fd5b9548d19aa8c12c29bc9637a96d42a06f8a5`).
- **What happens:** Each closed walk starts with the difference from its last point to its first, then visits the ordinary differences in order. It never compares the final ordinary difference against that closing difference. A genuine reversal there is accepted. Rotating the starting point of the same closed line changes the answer.
- **Trigger:** A closed line `(0,0,0), (0,1,0), (1,0,0), (2,0,0)`. Its ordinary differences are `(0,1,0), (1,-1,0), (1,0,0)` and its closing difference is `(-2,0,0)`. The only reversal is the omitted final comparison: 180 degrees, above the unchanged 170-degree criterion. The other three turns are 90, 135 and 45 degrees. Duplicate this line at z = 10 to make an eight-point grid with row size 4 and `bClosedRow = TRUE`; interleave each point with its duplicate at z = 10 to make row size 2 and `bClosedCol = TRUE`. The other grid direction is open and has only one step.
- **How found:** Issue [#2683](https://github.com/tangency-labs/solidcore/issues/2683), confirmed by a standalone probe of the pinned C++ library during [PR #2784](https://github.com/tangency-labs/solidcore/pull/2784).
- **Confirmed:** yes. `SmPointGrid::HasNonMonotonicPoints()` returns the following on the exact points above. `start` rotates the line before constructing either grid; all four starts describe the same closed line.

  | Grid direction | start 0 | start 1 | start 2 | start 3 | Correct result, every start |
  |---|---|---|---|---|---|
  | closed row | false | true | true | true | true |
  | closed column | false | true | true | true | true |

  Built against the precompiled reference `v1-3fd5b9548d19-726643fb17a52b7d`, with that revision's SMLib/SM_API headers, `-O2 -ffp-contract=off -fno-fast-math`, `NDEBUG`, `SM_GRAPHICS_CALLBACKS=1`, and the C++11 libstdc++ ABI. No C++ source changes.

The Rust regressions are the `issue_2683_*` tests in `surface/surface/assert_valid.rs`. This is separate from WR-737's criterion for whether a step carries a direction: completing a cycle must use the existing zero criterion. In the rounding-aware walk, repeated or rounding-only steps at the seam must also be passed over when comparing its first and last real directions.
