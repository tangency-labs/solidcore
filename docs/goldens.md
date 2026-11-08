# Golden data: where it lives and how a golden changes

The goldens (`expected/<unit>.jsonl`, `expected/known_failing.json`) and the per-unit fixtures
(`fixtures/<unit>/`) live in [tangency-labs/solidcore-goldens](https://github.com/tangency-labs/solidcore-goldens),
not in this repository. `goldens.pin` at the repository root names the one sha of it that this
commit is checked against. The SolidCraft case registry, `fixtures/solidcraft/`, stays here.

## Getting it

```sh
cargo xtask goldens-fetch           # solidcore-goldens at goldens.pin, into .goldens/ (gitignored)
cargo xtask goldens-fetch --check   # only check that the pinned sha is on solidcore-goldens
cargo xtask goldens-pin <sha>       # move the pin (checked), then fetch
```

The fetch is shallow and idempotent: a checkout already at the pin is left alone, one at another
sha is moved to it, and one with changes (a bless in progress) is refused, naming them.
`SOLIDCORE_GOLDENS_DIR` points the readers at another checkout instead of `.goldens/`
(`SOLIDCORE_GOLDENS_URL` points the fetch at a mirror).

Every reader goes through `crates/goldens-data`: the goldens binary, the units' fixtures and the
tests that read them. Missing data, or a `.goldens/` at another sha than the pin, is an error
naming `cargo xtask goldens-fetch`. Nothing skips: a test that needs data and has none fails.

## Changing a golden (the bless flow)

A golden proves a result unchanged, not correct. A change that moves one on purpose:

1. Bless into the checkout: `cargo run --profile goldens -p goldens -- --unit <id> --bless`
   writes `.goldens/expected/<id>.jsonl` (and a fixture generator writes `.goldens/fixtures/`).
2. Propose them: `cargo xtask goldens-bless-pr <branch> "<title>" <justification.md>` makes a
   branch of solidcore-goldens from the pin, commits every change in `.goldens/` (the title its
   subject, the file its message), pushes it and opens the pull request with `gh`, the file its
   description. The justification is what AGENTS.md asks for: every case that changed with its
   old and new values, why the new ones are right, and pictures before and after for geometry,
   and it links the solidcore pull request that needs it. By hand it is the same: `git -C
   .goldens switch -c <branch>`, `git -C .goldens add -A && git -C .goldens commit`, push the
   branch to solidcore-goldens, open the pull request.
3. Once it merges, `cargo xtask goldens-pin <merge sha>` in the solidcore pull request: it
   checks the sha is on solidcore-goldens, writes `goldens.pin` and fetches it. Commit the pin
   with the code change, and repeat the justification there (or link it).

CI checks the pin on every pull request (the `goldens-pin` job): a pin that is not on
solidcore-goldens fails it, so a solidcore pull request cannot land ahead of its data. The
goldens, the workspace tests and the site read the data at the pin.
