# Progress site

`site/build.py` generates the GitHub Pages site from data already in the repository;
nothing on it is maintained by hand except the short intro paragraph in the generator.

```sh
python3 site/build.py && python3 -m http.server -d bin/site 8000   # build (~3 s), then open http://localhost:8000
```

Python 3.11+ standard library only (`tomllib`). Output goes to `bin/site/` (gitignored);
`--out DIR` writes elsewhere. Output is deterministic: dates come from the commit
date of `HEAD` and from measured-report provenance.

| Page | Data |
|------|------|
| Overview (`index.html`) | `plan/units/*.toml` via `harness/plan/progress.py` (`summary`, `is_ready`), claims from the "In progress" section of `plan/PROGRESS.md`, `reports/ssi/ssi-corpus.json`, and the plan along `main`'s first-parent git history for the timeline (`plan/units/`, or `plan/modules.toml` before the per-unit split). The history is read once; overall, per-module and per-unit series come from the same pass |
| Module pages (`areas/<module>.html`) | `plan/units/` and its history: the same chart as the overview, over that module's units only |
| Units (`units/index.html`) | `plan/units/*.toml`, `harness/callgraph/units.json`, `reports/difftest/*.json` (`cases_*`, `deferred_by_unit`; copied to `reports/<unit>.json` on the site). Reports are not committed: CI fetches main's latest set with `scripts/fetch-reports.sh` before building; run it locally too |
| Unit pages (`units/<id>.html`) | the above, plus the `//! | Rust | Upstream |` table in the module's file header and the unit's status history from the timeline |
| Dependency graph (`graph.html`, `assets/graph.json`) | `plan/units/*.toml` (`depends_on`, module, status) and claims, laid out as one force-directed cluster per module; lists of ready and most-blocking units are rendered server-side, the graph by `site/assets/graph.js` |
| Upstream issues (`issues.html`) | `plan/upstream-issues/*.md`, one `#### <SERIES>-<n>.` entry per file, in name order; the intro is that directory's `README.md` |
| SSI corpus (`ssi.html`) | `reports/ssi/ssi-corpus.json` and `ssi-provenance.json`, fetched from the `SSI corpus` workflow with `python3 scripts/ssi_report.py fetch`. The four engines have separate scores; only C++/Rust are compared for parity. The page shows the measured SHA, SMLib pin, date and run link; absent or invalid artifacts show “No SSI report available” |
| Solids (`solids.html`; tile on the overview) | `reports/solids/solids.json`: the `primitives` example's `--json` output (schema in `crates/solidcore/examples/primitives/main.rs`; copied to `reports/solids.json` on the site). Not committed: CI's test job uploads it from `main` as the `solids-report` artifact and the Pages build fetches it with `scripts/fetch-solids.sh`. Without it the page says "no data yet". The "Corpus (real models)" section is a placeholder until `corpus.models` is filled |

Text from the repository is HTML-escaped before a small markdown subset (inline code,
bold, italic, http(s) links, nested `- ` lists) is applied, so raw HTML written in a PR
never reaches the page. Charts are inline SVG; `site/assets/app.js` only sorts and
filters tables, drives the chart tooltips and the overview's module picker (without JS the overview
shows the all-units chart). The only third-party code is d3-force (ISC), vendored under
`site/assets/vendor/` and used only by the graph page (see `ATTRIBUTION.md`). No external assets.

Deployment: `.github/workflows/pages.yml` rebuilds on every push to `main`, after
completed `CI` or `SSI corpus` runs on `main` (fresh measured reports), hourly
(claims from open `port:` PRs) and on manual dispatch. It first generates the unit
checklist (`harness/plan/progress.py --prs`), which
the site publishes as `PROGRESS.md`; nothing is committed back to `main`. GitHub Pages
must be enabled with source "GitHub Actions".

The SSI benchmark runs daily at 03:23 UTC and on manual dispatch of
`.github/workflows/ssi-corpus.yml` on `main`. It uploads the complete 119-case,
four-engine report and its provenance as `ssi-corpus-report` (90 days); raw engine
outputs and logs are retained for 14 days. The report JSON stays deterministic;
measurement time and source/upstream/corpus revisions live in the provenance file.
Pages downloads only artifacts from successful runs of that workflow in this
repository on `main`, validates their summaries and binds provenance to the Actions
run SHA. It never falls back to the old committed SSI snapshot. Both JSON files
are downloadable from the generated site under `reports/`.

To fetch the latest measured SSI results before a local site build:

```sh
python3 scripts/ssi_report.py fetch
python3 site/build.py
```

The site revision and measurement revision may differ between nightly runs. An
upstream pin change is also shown explicitly until a new benchmark is available.
After first merging the workflow, dispatch `SSI corpus` once to populate the site.
