"""Validate and hash the immutable #3110 baseline artifact directory."""
import argparse
import hashlib
import json
from pathlib import Path

PINNED = {
    "cf110-inputs.jsonl": "eb9347c55dd6d8b1c503f09db342adc0608829057dda36639c3da402bcb29398",
    "cf110-before/outputs.jsonl": "dcc2cd45d4071d9b0fd3efc200717c96c4ecf9cb8cdb96f812865932246a6a42",
}


def require(condition, message):
    """Reject invalid audit evidence even when Python assertions are disabled."""
    if not condition:
        raise ValueError(message)


def rows(path):
    """Read nonempty JSONL records in their original order."""
    return [json.loads(line) for line in path.read_text(encoding="utf-8-sig").splitlines() if line.strip()]


def focal(value):
    """Select the native curve/edge-face Single projected maximum fixtures."""
    return (value.get("op") == "static" and value.get("operation") == 15
            and value.get("request") == 0 and value.get("entry") in (6, 7))


def validate(root, write=True):
    """Check pinned baseline bytes, cohort identity and decoder completeness."""
    root = Path(root)
    hashes = {name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in PINNED}
    require(hashes == PINNED, "immutable baseline hash mismatch")
    inputs = rows(root / "cf110-inputs.jsonl")
    outputs = rows(root / "cf110-before/outputs.jsonl")
    require(len(inputs) == len(outputs) == 2575, "baseline case count changed")
    require([(r["unit"], r["case"], r["input"]) for r in inputs]
            == [(r["unit"], r["case"], r["input"]) for r in outputs], "input/output cases differ")
    require(len({(r["unit"], r["case"]) for r in inputs}) == 2575, "duplicate baseline case")
    require(all(r["unit"] == "intersect.topology_solver_p2" for r in inputs), "unexpected unit")
    require(sum(r["case"].startswith("gen") for r in inputs) == 600, "generated cohort changed")
    require(all(r["focal_native_pair"] == focal(r["input"]) for r in inputs + outputs), "focal flag mismatch")
    native = [r for r in outputs if focal(r["input"])]
    require(len(native) == 14, "native cohort changed")
    require(not any("decode_error" in r for r in outputs), "baseline decoder failed")
    require(all(r.get("output", {}).get("values", {}).get("status") == [0.0]
                and len(r.get("solutions") or []) == 1 for r in native), "native baseline is not fourteen successes")
    manifest = {
        "baseline_revision": "cf11091a9", "rustc": "1.98.0", "cases": 2575,
        "native_pair_cases": 14, "generated_tree_controls": 600, "fixture_cases": 1975,
        "inputs_match_outputs": True, "decode_errors": 0,
        "harness_errors": sum("harness_error" in r for r in outputs), "sha256": hashes,
    }
    for name in ("cf110-before/intersect.topology_solver_p2.debug", "affine_face_output_audit_3110.rs", "audit-cf110"):
        if (root / name).exists():
            manifest["sha256"][name] = hashlib.sha256((root / name).read_bytes()).hexdigest()
    if write:
        (root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        for name, data in (("cf110-native14-inputs.jsonl", [r for r in inputs if focal(r["input"])]),
                           ("cf110-before/native14-outputs.jsonl", native)):
            (root / name).write_text("".join(json.dumps(r, separators=(",", ":")) + "\n" for r in data), encoding="utf-8")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_directory", type=Path)
    print(json.dumps(validate(parser.parse_args().artifact_directory), indent=2))
