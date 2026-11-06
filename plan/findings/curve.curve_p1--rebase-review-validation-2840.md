# PR #2840 rebase and review-test validation

Rebased onto main `513e6fcfcd049d41e9487599880479de4add63f2` at the user's
explicit request. The pre-rebase head was `ff3089ca9`; the rebased code head is
`65f2f34fd`, followed by review tests and helper documentation in `0b9ffdd28`.

The sole conflict was `plan/units/curve.curve_p5.toml`: main added the curve-length
accuracy findings and changed its status to partial. Resolution preserves those
findings while retaining this PR's blocked status for its unresolved projected
query disagreements. Range-diff confirms the other nineteen replayed commits
retain their patches. The three modified production solver files were identical
before and after rebase; subsequent edits to them add helper documentation only.

## Checks

- Workspace formatting, workspace/all-target clippy with warnings denied, and
  the wasm check pass on the rebased code.
- Library tests: **3,569 passed, four ignored**.
- Focused regression suites: **20 passed**, with one diagnostic export ignored.
- Both accuracy-refusal tests now require caller-facing `NotWithinTolerance`.
- The additional public-query regression covers ordinary/reversed bounded lines
  with distinct natural and analytic domains, and ordinary/reversed infinite
  lines. Evaluation and projected extrema agree with independent closed forms.
- The final review-test file was checked against the executed file by SHA-256;
  final formatting and targeted clippy pass after the review changes.

The full integration suite remains incomplete. The earlier native comparison
belongs to its recorded older base; this rebase does not establish a fresh native
comparison or resolve the documented zero-tolerance, domain, surface, curve-pair
and trim questions. This remains a draft with the affected units blocked.
