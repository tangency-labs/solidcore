# heal.gap: default append changes maximum state (2026-11-05)

PR #2803, follow-up to #2760 and urgent main sweep #2794. Seven ordinary diagnostic mismatches are recorded; none is parked or deferred. The expected status is partial, not green. A current-head native run remains pending.

## Evidence and limits

The official push/main CI run [41965907239, attempt 1](https://github.com/tangency-labs/solidcore/actions/runs/41965907239), head `4b10d7d6e22382e37e8239b0f395c7326807e6c6`, failed: `heal.gap` had 749 cases, 742 passed, 7 failed, 0 deferred. This is the old report, not a current-head run or a passed baseline. The sweep command was `cargo run --profile difftest -p difftest --features cpp -- --green --jobs 0`; 16 seeds, at least 400 generated cases, relative tolerance `1e-9`, no absolute floors. Artifact 13388006545 ZIP SHA-256 `8956146fd500301af25c57cdd292d65bd6ca1a41c67a8763a38318da31629966`; its `heal.gap.json` SHA-256 `0e7bd4c0f666e19951f1fbacdc12678fb6919bcb77854d8879ff7349c34acd54`. The report records mismatches, not complete outputs. The archive has no provenance JSON; the source pin is `3fd5b9548d19aa8c12c29bc9637a96d42a06f8a5` (`scripts/upstream-pin.env`) and CI uses the C++ reference action. Neither a build manifest nor binary reproducibility is certified here; no local C++ execution is claimed. The failed artifact is diagnostic only, never a baseline.

#2760 merged as `3b0019128086cad4781a64f7f86985ede2503fa6`; its first parent is `159cd0b849a6faf43f83b11730e61958e5685657`. The previous default Rust path can be reproduced from `git show 159cd0b849a6faf43f83b11730e61958e5685657:crates/solidcore/src/heal/gap.rs` (Git blob `af94b20c7055abff5d5d12ba1a1c31ff73ea3afe`), specifically `GapArray::append` default behavior. #2760 added regression tests before changing append to sequential `add`; see `plan/upstream-issues/WR-90--addendum-2026-11-05-default-append.md`. This is source-path reproduction, not a new run of that parent.

## Independent input oracle

The runner builds gaps from absolute input values in order and appends A to destination B (`gb.append(&ga)`). Before append the destination maximum is the first argmax of `abs(B)`, or index -1 and length 0 if B is empty. The old default append leaves that cache unchanged; the new sequential add gives the first argmax of `abs(B)` followed by `abs(A)`. Every input here is finite, and each combined array has a unique maximum; B's first maximum is used when B ties. Therefore the following 14 reported mismatch components are independently obtained from the inputs, not by taking the Rust output as the oracle. `values.state[4]` is a serialized float index, not an exact integer field; `values.state[5]` is length. Strict `>` retains the earlier maximum on finite ties. Nonfinite values are outside these seven cases: no general NaN policy is inferred from them; sequential insertion follows the existing `add` comparison behavior.

| case | source | C++ state[4], state[5] | Rust state[4], state[5] |
|---|---|---|---|
| `gap_array_append` | `fixture:fixtures/heal.gap/upstream.json#37 (upstream source/SMLib/inc/SmGap.h:875)` | 0e0, 1e-1 | 3e0, 2e0 |
| `gen0017` | `generated:seed=6223471294032734338` | 4e0, 7.215699039766319e-1 | 5e0, 9.029522644380432e-1 |
| `gen0179` | `generated:seed=6223471294032734500` | 1e0, 8.196507478065342e-1 | 6e0, 8.932539832416029e-1 |
| `gen0289` | `generated:seed=6223471294032734610` | 0e0, 4.130508541544686e-1 | 3e0, 9.20731530769151e-1 |
| `gen0369` | `generated:seed=6223471294032734690` | 0e0, 7.557292741382693e-1 | 4e0, 8.599124624438688e-1 |
| `gen0389` | `generated:seed=6223471294032734710` | 2e0, 8.792184173358797e-1 | 8e0, 9.529596725468907e-1 |
| `gen0393` | `generated:seed=6223471294032734714` | -1e0, 0e0 | 0e0, 8.995292105893602e-1 |

The exact inputs are keyed by case name below. Canonical input SHA-256 is `1a09fdfc3d9b7512a2cff17ef0d7bb7992dc6e99663b4c92d5c541f74217767c`: UTF-8 bytes of the single JSON line, with keys sorted recursively, compact separators `,` and `:`, ASCII escaping, and no trailing newline.

```json
{"gap_array_append":{"a":[0.5,-2,1],"b":[0.1,0.1],"op":"gap_array"},"gen0017":{"a":[0.9029522644380432,0.38458932966695003,0.2566099294854427,0.06175405234114262],"b":[0.31739121473222864,0.3759323387661109,0.6181548613594894,0.20927078507788166,-0.7215699039766319],"op":"gap_array"},"gen0179":{"a":[0.6002375791562147,0.8932539832416029],"b":[0.25726814861949165,0.8196507478065342,0.09811641206343635,0.5828907401188577,0.22796898554308842],"op":"gap_array"},"gen0289":{"a":[0.5686734207377255,-0.920731530769151,-0.29513502147103376,-0.7244026310902894],"b":[-0.4130508541544686,0.19864895566768404],"op":"gap_array"},"gen0369":{"a":[-0.8599124624438688,0.5137855628350649],"b":[-0.7557292741382693,0.6905428569134193,-0.3560720485350024,0.05134780441645104],"op":"gap_array"},"gen0389":{"a":[-0.8735662808114504,0.7573246285521844,0.5398312625850443,0.9529596725468907,0.7950400257368389],"b":[-0.48381523102265267,-0.12364977104895614,-0.8792184173358797,-0.5937021579506383,0.14200498750575852],"op":"gap_array"},"gen0393":{"a":[0.8995292105893602,0.4262820965714187],"b":[],"op":"gap_array"}}
```

The observed stale C++ values are consistent with WR-90; this record does not add a `cpp_defect` marker or claim new confirmation of the reference build. Ordinary known entries still run and must continue to mismatch or fail as stale. No case, tolerance, comparator, seed or fixture is changed.
