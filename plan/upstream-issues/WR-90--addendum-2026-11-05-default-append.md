Updates: WR-90

#### WR-90 (update, 2026-11-05). Correct maximum on every append

- `GapArray::append` now maintains the earliest longest gap unconditionally, using the cached maxima of both arrays. The WR-90 flag and the public `append_opt_in` and `append_with_kernel` methods are removed; this is a source-breaking API cleanup.
- The independent Rust regression scans stored nonnegative finite lengths, including empty arrays, ties, extreme scales and a chained append. For `[0.5]` appended with `[1, 3, 2]`, the correct maximum is length 3 at index 2 rather than the previous length 0.5 at index 0.
- The existing WR-90 record's C++ confirmation is historical, not a new native execution. The `heal.gap` diagnostic emits the appended array's maximum index and length; cases with a longer appended maximum may now disagree with C++. No diagnostic cases or expectations were removed or marked, and this change has not itself run the C++ diagnostic.
