#!/usr/bin/env python3
"""Parked difftest cases against upstream's latest C++ (the parked-upstream workflow).

A known-mismatch entry with `cpp_defect = "<id>"` is parked: the gate does not run it
(crates/difftest/README.md, "Known mismatches and parked C++ defects"). When upstream
usd-brep changes (a new tag or a new HEAD), the workflow builds the C++ reference from
that commit and runs `difftest --run-parked` on the units this script lists; a parked case
that now agrees with Rust is reported as an ordinary `pass`, and this script writes the
title and body of the one tracking issue,
"usd-brep <ref> released: N parked cases now agree, adopt?" (or "... none changed").

A `pass` alone is not adoption evidence for a nondeterminism (`ND-*`) defect: some C++
results depend on heap addresses (ND-2, ND-258), so a case can pass in one run and fail
in the next (issue #1405 had five such parked sweep cases advertised as agreement
although nothing had changed). For units with ND-* cases the workflow reruns C++ with
the pads below and `difftest` records per-case `heap_qualification`. An ND-* pass is
eligible only with complete evidence for the same unit and case, every pad observed
on a really padded forked child, and no divergence. Every report must identify the
requested unit. Non-ND units keep their original single-run path; valid non-ND reports
keep the plain `pass` rule.

  parked_upstream.py units [--root DIR]
      the unit ids with at least one parked case, one per line
  parked_upstream.py pads [--unit UNIT] [--root DIR]
      the required heap pads (comma-separated); empty for a unit with no ND-* cases
  parked_upstream.py issue --report-dir DIR --upstream-ref REF --upstream-sha SHA
                          --title-out FILE --body-out FILE [--root DIR]
      writes the issue title and body; exit 2 when a unit's report is missing (the run
      did not produce it)
"""

import argparse
import json
import pathlib
import sys
import tomllib

KNOWN_DIR = "crates/difftest/known_mismatches"

# Every title starts with this and contains TITLE_MARK, so the workflow finds the one
# open tracking issue again whatever the ref and count.
TITLE_PREFIX = "usd-brep "
TITLE_MARK = " released: "

# The maintainer-recommended pad set on #1405 includes the native-self flips observed
# at the pinned commit. Qualification requires all of these observations, but is not
# a proof that every possible heap layout is address-independent.
REQUIRED_PADS = (1, 2, 3, 5, 8, 13, 21, 34, 55, 89)

# Defects whose result depends on memory addresses, so a single pass proves nothing.
HEAP_DEPENDENT_PREFIX = "ND-"


def parked_cases(root):
    """{unit: {case: cpp_defect id}} for every entry that has a `cpp_defect`."""
    out = {}
    for path in sorted((root / KNOWN_DIR).glob("*.toml")):
        with path.open("rb") as f:
            data = tomllib.load(f)
        for entry in data.get("mismatch", []):
            defect = entry.get("cpp_defect")
            if defect:
                out.setdefault(entry["unit"], {})[entry["case"]] = defect
    return out


def heap_qualified(report, case):
    """True when `report` carries complete `--heap-pads` evidence for `case`.

    Complete means: the report's pad set is exactly `REQUIRED_PADS`; one observation
    names this case; every pad was rerun (`observed`) against a forked child that really
    padded (`padded`); and no rerun diverged or failed. A report without the block (an
    older run, or one made without pads) fails closed.
    """
    q = report.get("heap_qualification")
    if (
        not isinstance(q, dict)
        or q.get("pads") != list(REQUIRED_PADS)
        or any(type(pad) is not int for pad in q["pads"])
    ):
        return False
    observations = q.get("cases")
    if not isinstance(observations, list):
        return False
    entries = [
        o for o in observations if isinstance(o, dict) and o.get("case") == case
    ]
    if len(entries) != 1:
        return False
    o = entries[0]
    return (
        type(o.get("observed")) is int
        and type(o.get("padded")) is int
        and o["observed"] == len(REQUIRED_PADS)
        and o["padded"] == len(REQUIRED_PADS)
        and o.get("divergence") is None
    )


def agreeing(report_dir, parked):
    """[(unit, case, id)] adoptable parked cases and the units without a report.

    A case is adoptable when its `--run-parked` outcome is `pass` and, for an `ND-*`
    defect, its report carries complete heap qualification. Passing cases without
    matching report identity or the required heap qualification are returned separately
    and must not be advertised as agreement.
    """
    found, missing, unqualified = [], [], []
    for unit, cases in sorted(parked.items()):
        path = report_dir / f"{unit}.json"
        if not path.is_file():
            missing.append(unit)
            continue
        report = json.loads(path.read_text())
        for c in report.get("cases", []):
            if c.get("outcome") != "pass" or c.get("name") not in cases:
                continue
            defect = cases[c["name"]]
            if report.get("unit") != unit or (
                defect.startswith(HEAP_DEPENDENT_PREFIX)
                and not heap_qualified(report, c["name"])
            ):
                unqualified.append((unit, c["name"], defect))
                continue
            found.append((unit, c["name"], defect))
    return found, missing, unqualified


def title(found, ref):
    if not found:
        return f"{TITLE_PREFIX}{ref}{TITLE_MARK}none changed"
    n = len(found)
    cases = "case now agrees" if n == 1 else "cases now agree"
    return f"{TITLE_PREFIX}{ref}{TITLE_MARK}{n} parked {cases}, adopt?"


def body(found, parked, unqualified, ref, sha):
    total = sum(len(c) for c in parked.values())
    lines = [
        f"Upstream usd-brep changed (`{ref}`, commit `{sha}`). The parked-upstream "
        f"workflow built the C++ reference from it and ran the {total} parked difftest "
        "case(s) (known mismatches with a `cpp_defect`).",
        "",
    ]
    if found:
        lines += [
            "These now agree with Rust:",
            "",
            "| unit | case | upstream issue |",
            "|---|---|---|",
        ]
        lines += [f"| `{u}` | `{c}` | `{i}` |" for u, c, i in found]
        lines += [
            "",
            "To adopt a fix once our pin moves to that upstream commit: delete the case's "
            "`cpp_defect` line (or its whole entry) in "
            "`crates/difftest/known_mismatches/<unit>.toml`, and update the upstream issue "
            "record in `plan/upstream-issues/`.",
        ]
    elif unqualified:
        lines.append(
            "Nothing is adoptable yet: the case(s) below passed but lack matching "
            "report identity or required heap qualification."
        )
    else:
        lines.append("None of them agrees with Rust yet: nothing to adopt.")
    if unqualified:
        lines += [
            "",
            "Passed once, not adopted (report identity mismatch or `ND-*` heap "
            "qualification incomplete; this is not agreement evidence for the same "
            "unit and case): "
            + ", ".join(f"`{u}`/`{c}` ({i})" for u, c, i in unqualified)
            + ".",
        ]
    lines += ["", "This issue is rewritten each time upstream changes."]
    return "\n".join(lines) + "\n"


def main(argv=None):
    # (no docstring under `python3 -OO`)
    summary = (__doc__ or "").strip().splitlines()
    p = argparse.ArgumentParser(description=summary[0] if summary else None)
    p.add_argument("command", choices=["units", "pads", "issue"])
    p.add_argument("--root", type=pathlib.Path, default=pathlib.Path("."))
    p.add_argument("--unit", help="only request pads when this unit has an ND-* case")
    p.add_argument("--report-dir", type=pathlib.Path)
    p.add_argument("--upstream-ref", default="HEAD")
    p.add_argument("--upstream-sha", default="unknown")
    p.add_argument("--title-out", type=pathlib.Path)
    p.add_argument("--body-out", type=pathlib.Path)
    a = p.parse_args(argv)
    parked = parked_cases(a.root)
    if a.command == "pads":
        if a.unit is None or any(
            defect.startswith(HEAP_DEPENDENT_PREFIX)
            for defect in parked.get(a.unit, {}).values()
        ):
            print(",".join(str(n) for n in REQUIRED_PADS))
        return 0
    if a.command == "units":
        for unit in sorted(parked):
            print(unit)
        return 0
    if a.report_dir is None or a.title_out is None or a.body_out is None:
        p.error("issue needs --report-dir, --title-out and --body-out")
    found, missing, unqualified = agreeing(a.report_dir, parked)
    if missing:
        print("parked_upstream: no report for " + ", ".join(missing), file=sys.stderr)
        return 2
    a.title_out.write_text(title(found, a.upstream_ref) + "\n")
    a.body_out.write_text(body(found, parked, unqualified, a.upstream_ref, a.upstream_sha))
    return 0


if __name__ == "__main__":
    sys.exit(main())
