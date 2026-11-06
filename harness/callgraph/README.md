# harness/callgraph

Function-level call graph of upstream SMLib + SM_API and the port plan derived from it.
Owner: agent B. Only licence-cleared trees are parsed (`source/SMLib`, `source/SM_API`,
`tests/SM_API_test`, `tests/prog_test`).

```
harness/callgraph/run.sh          # full regeneration (about 3-4 min on 6 cores)
python3 harness/callgraph/why.py 'SmFilletIntersector::FlushCurve'   # shortest call path from an entry point
python3 harness/callgraph/why.py SmPoly.cpp --group boolean           # restrict to one entry-point group
```

| step | script | output (bin/callgraph/ unless noted) |
|---|---|---|
| 1 | `extract.py` | `tu/*.json`: per-TU definitions, call edges, classes, globals (libclang with `compile_commands.py` analysis metadata) |
| 2 | `build_graph.py` | `graph.json`, `closure.json` (CHA + RTA closures from `entrypoints.py` roots), `files.json`, `resolve.json`, `overriders.json` |
| 3 | `plan_units.py` | `plan/units/<id>.toml`, `plan/order.md`, `harness/callgraph/units.json`, `plan_summary.json` |
| 4 | `write_closure.py` | `plan/closure.md`, `harness/callgraph/unused_functions.txt` |

Other files: `entrypoints.py` (in-scope SM_API entry points and kernel roots), `pp.py` (`#if`
evaluator for text scans), `scc_probe.py` (function-level SCC statistics), `why.py`.

`units.json` maps every unit id in `plan/units/` to the exact upstream functions (qualified
name, file, line, line count, clang USR) it owns. A unit that names a file in `upstream` owns only
the listed functions of that file; other functions of the same file belong to other units or are
out of scope (`unused_functions.txt`).

Tooling: libclang 21 from Homebrew LLVM (`LIBCLANG`, `CLANG_RESOURCE_DIR` override the paths) with
the `libclang` pip bindings in `bin/callgraph-venv`; libstdc++ 13 headers (`CXX_SYS_INCLUDES`).

Extraction needs no CMake configuration, `bin/cpp/`, or C++ build/link/run. The Python-only
`compile_commands.py` enumerates `src/*.cpp` directly in the four cleared trees above,
preserving the former release definitions (`NDEBUG`, `NOMINMAX`, `SM_GRAPHICS_CALLBACKS=1`,
`_GLIBCXX_USE_CXX11_ABI=1`), C++17, and each target's include order. SM_API sources add
`SM_API_EXPORTS`; `SM_API_test.cpp` adds `SM_API_TEST_STANDALONE`. That test root is analysed
once, not both with and without `main` as in the old database. The two tests requiring
uncleared Examples sources remain excluded; deleted generated drivers/stubs are not recreated.
Missing or empty source groups refuse before parsing. An `upstream/` symlink is supported;
symlinked source-group directories or individual sources are refused.

Direct extraction constructs the metadata in memory. To inspect it without libclang or a
compiler, run `python3 harness/callgraph/compile_commands.py` (optional `--root PATH`). It prints
Clang's [JSON compilation database](https://clang.llvm.org/docs/JSONCompilationDatabase.html)
`arguments` arrays, keeping spaces unescaped and intact, and never executes those commands.
The existing libclang library/resource/standard-header configuration is still required to
parse an AST. Metadata generation alone is not evidence of an error-free parse or unchanged graph.

The full wrapper also rewrites plan files; use direct `extract.py --filter SUBSTR --jobs N`
for bounded analysis, and do not regenerate the shared plan as part of an unrelated fix.
