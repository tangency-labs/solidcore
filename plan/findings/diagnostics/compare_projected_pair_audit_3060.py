#!/usr/bin/env python3
"""Compare exact-input JSONL audits without accepting changed outputs as goldens."""
import argparse
import json
from pathlib import Path


def load(path):
    rows = {}
    for line in Path(path).read_text(encoding="utf-8-sig").splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        key = (row["unit"], row["case"])
        if key in rows:
            raise ValueError(f"duplicate case {key}")
        rows[key] = row
    return rows


def status(row):
    if row.get("harness_error"):
        return {"kind": "harness_error", "reason": row["harness_error"]}
    out = row.get("output", {})
    for kind in ("error", "rejected", "deferred"):
        if out.get(kind) is not None:
            return {"kind": kind, "reason": out[kind]}
    code = out.get("ints", {}).get("status", out.get("values", {}).get("status"))
    return {"kind": "success" if code == [0] else "status", "code": code}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("before")
    parser.add_argument("after")
    args = parser.parse_args()
    before, after = load(args.before), load(args.after)
    if before.keys() != after.keys():
        raise ValueError("before/after case sets differ")
    counts = {"cases": len(before), "changed": 0, "new_refusals": 0, "decode_errors": 0}
    for key in sorted(before):
        old, new = before[key], after[key]
        if old["input"] != new["input"]:
            raise ValueError(f"inputs changed for {key}")
        old_status, new_status = status(old), status(new)
        changed = old.get("output") != new.get("output") or old.get("harness_error") != new.get("harness_error")
        refused = old_status["kind"] == "success" and new_status["kind"] != "success"
        decode_error = old.get("decode_error") or new.get("decode_error")
        counts["changed"] += changed
        counts["new_refusals"] += refused
        counts["decode_errors"] += bool(decode_error)
        print(json.dumps({"unit":key[0], "case":key[1], "changed":changed,
            "new_refusal":refused, "before_status":old_status, "after_status":new_status,
            "before_solutions":old.get("solutions"), "after_solutions":new.get("solutions"),
            "decode_error":decode_error}, sort_keys=True))
    print(json.dumps({"summary":counts}, sort_keys=True))


if __name__ == "__main__":
    main()
