"""Replay #3110 using existing rlibs; never starts Cargo or changes goldens.

All compilation scratch, executables, logs and outputs stay under artifact-dir,
which must be inside the caller's explicitly supplied storage-root.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess

from freeze_affine_face_audit_3110 import require, rows, validate
from compare_projected_pair_audit_3060 import status


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "storage-root", "artifact-dir", "deps", "goldens-rlib", "solidcore-rlib", "serde-json-rlib"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--case-root", type=Path, help="fixture root, defaults to repo; use immutable input-root for baseline")
    parser.add_argument("--revision", required=True, help="revision of the already-built kernel")
    parser.add_argument("--proof-revision", required=True, help="reviewed native residual source revision")
    parser.add_argument("--label", required=True, help="new artifact subdirectory, never cf110-before")
    parser.add_argument("--rustc", default="rustup run 1.98.0 rustc")
    parser.add_argument("--reference-outputs", type=Path, help="optional immediate-parent frozen replay for attribution")
    parser.add_argument("--expect-certified", action="store_true")
    parser.add_argument("--expect-baseline", action="store_true")
    args = parser.parse_args()
    root = args.artifact_dir.resolve()
    require(root.is_relative_to(args.storage_root.resolve()), "artifacts must be on designated storage")
    require(args.label not in ("", ".", "..", "cf110-before") and Path(args.label).name == args.label,
            "invalid or reserved label")
    validate(root)
    source = Path(__file__).resolve().parent
    out = root / args.label
    out.mkdir(exist_ok=False)
    scratch = out / "scratch"
    scratch.mkdir()
    environment = dict(os.environ, TMPDIR=str(scratch), TMP=str(scratch), TEMP=str(scratch))
    adapter = subprocess.check_output(["git", "-C", str(args.repo), "show", args.revision
        + ":crates/goldens/src/units/intersect_topology_solver_p2.rs"])
    (out / "static_adapter_source.rs").write_bytes(adapter)
    adapter_hash = hashlib.sha256(adapter).hexdigest()
    require(adapter_hash == "ab07e2085d02e1025ea981d9456b54fd819e0feb42b32b0471752f7eeb83dcc8",
            "static adapter constructor contract changed; independently re-audit before certifying")
    contract = {"adapter_source_sha256": adapter_hash,
        "curve_controls": [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        "curve_degree": 1, "curve_knots": [0.0, 0.0, 1.0, 1.0],
        "curve_interval": [0.0, 1.0], "normal": [0.0, 0.0, 1.0],
        "best_answer_so_far_sq": 0.0,
        "origin": "pinned static_rust adapter; constants are not fields in fixture JSON"}
    (out / "adapter-contract.json").write_text(json.dumps(contract, indent=2) + "\n", encoding="utf-8")
    original = subprocess.check_output(["git", "-C", str(args.repo), "show", args.proof_revision
        + ":crates/solidcore/src/intersect/topology_solver/projected_residual.rs"])
    adapted = original.decode().replace(
        "use crate::blend::fillet_section::certificate::knot_cuts;", "use crate::knot_cuts;"
    ).replace("use crate::", "use solidcore::").replace("use solidcore::knot_cuts;", "use crate::knot_cuts;")
    (out / "projected_residual_audit_copy.rs").write_text(adapted, encoding="utf-8")
    require((out / "projected_residual_audit_copy.rs").read_text(encoding="utf-8") == adapted,
            "residual import-only copy changed")
    proof_commit = subprocess.check_output(["git", "-C", str(args.repo), "rev-parse",
        "--verify", args.proof_revision + "^{commit}"], text=True).strip()
    provenance = {"operator_supplied_kernel_revision": args.revision,
        "rlib_origin_verified_by_driver": False, "proof_source_commit": proof_commit,
        "adapter_source_sha256": adapter_hash, "adapter_contract_sha256": digest(out / "adapter-contract.json"),
        "residual_source_sha256": hashlib.sha256(original).hexdigest(),
        "residual_import_adapted_sha256": digest(out / "projected_residual_audit_copy.rs"),
        "goldens_rlib_sha256": digest(args.goldens_rlib), "solidcore_rlib_sha256": digest(args.solidcore_rlib),
        "serde_json_rlib_sha256": digest(args.serde_json_rlib)}
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    compiler = ["nice", "-n", "10"] + shlex.split(args.rustc)
    for stem, externs in (("affine_face_output_audit_3110", {"goldens": args.goldens_rlib}),
                          ("affine_face_native_audit_3110", {"solidcore": args.solidcore_rlib})):
        copied = out / (stem + ".rs")
        copied.write_bytes((source / copied.name).read_bytes())
        provenance[stem + "_source_sha256"] = digest(copied)
        (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
        command = compiler + ["--edition", "2021", str(copied), "-L", "dependency=" + str(args.deps.resolve()),
                              "--extern", "serde_json=" + str(args.serde_json_rlib.resolve()), "-o", str(out / stem)]
        for name, path in externs.items():
            command += ["--extern", name + "=" + str(path.resolve())]
        with (out / (stem + ".compile.log")).open("w", encoding="utf-8") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=environment, check=True)
        provenance[stem + "_executable_sha256"] = digest(out / stem)
        (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    executable = out / "affine_face_output_audit_3110"
    with (out / "current-inputs.jsonl").open("w", encoding="utf-8") as output:
        subprocess.run(["nice", "-n", "10", str(executable), str(args.case_root or args.repo), "freeze"], stdout=output, env=environment, check=True)
    frozen = rows(root / "cf110-inputs.jsonl")
    current = rows(out / "current-inputs.jsonl")
    require(frozen == current and (root / "cf110-inputs.jsonl").read_bytes()
            == (out / "current-inputs.jsonl").read_bytes(),
            "final fixture/generator inputs differ from frozen bytes; inspect before replay")
    with (out / "outputs.jsonl").open("w", encoding="utf-8") as output:
        subprocess.run(["nice", "-n", "10", str(executable), str(args.repo), "replay", str(root / "cf110-inputs.jsonl")],
                       stdout=output, env=environment, check=True)
    before = rows(root / "cf110-before/outputs.jsonl")
    after = rows(out / "outputs.jsonl")
    require([(r["unit"], r["case"], r["input"]) for r in before]
            == [(r["unit"], r["case"], r["input"]) for r in after], "replay cases changed")
    if args.expect_baseline:
        require(before == after and digest(out / "outputs.jsonl")
                == digest(root / "cf110-before/outputs.jsonl"), "baseline replay is not bit-identical")
    changes = []
    for old, new in zip(before, after):
        changes.append({"case": old["case"], "changed": old.get("output") != new.get("output")
            or old.get("harness_error") != new.get("harness_error"),
            "new_refusal": status(old)["kind"] == "success" and status(new)["kind"] != "success",
            "before_status": status(old), "after_status": status(new),
            "decode_error": new.get("decode_error"), "before_solutions": old.get("solutions"),
            "after_solutions": new.get("solutions")})
    (out / "comparison.jsonl").write_text("".join(json.dumps(r) + "\n" for r in changes), encoding="utf-8")
    with (out / "native.jsonl").open("w", encoding="utf-8") as output:
        subprocess.run(["nice", "-n", "10", str(out / "affine_face_native_audit_3110"), str(out / "outputs.jsonl"), str(out / "adapter-contract.json")],
                       stdout=output, env=environment, check=True)
    native = rows(out / "native.jsonl")
    require(len(native) == 14 and {r["case"] for r in native}
            == {r["case"] for r in after if r["focal_native_pair"]}, "native audit case set changed")
    obligations = ("admissible", "report_in_own_enclosure", "native_metric_in_own_enclosure",
                   "left_domain_proved", "right_domain_proved", "reduced_threshold_refused")
    failed = [r["case"] for r in native if "native_audit_error" in r
              or not all(r.get("native", {}).get(key) is True for key in obligations)]
    summary = {"cases": len(after), "changed": sum(r["changed"] for r in changes),
        "new_refusals": sum(r["new_refusal"] for r in changes),
        "decode_errors": sum(r["decode_error"] is not None for r in changes),
        "native_cases": len(native), "native_unproved_or_invalid": failed,
        "after_outputs_sha256": digest(out / "outputs.jsonl"), "provenance": provenance}
    if args.reference_outputs:
        reference = rows(args.reference_outputs)
        require([(r["unit"], r["case"], r["input"]) for r in reference]
                == [(r["unit"], r["case"], r["input"]) for r in after], "immediate-parent cases differ")
        parent_changes = []
        for old, new in zip(reference, after):
            parent_changes.append({"case": old["case"],
                "changed": old.get("output") != new.get("output")
                    or old.get("harness_error") != new.get("harness_error"),
                "new_refusal": status(old)["kind"] == "success" and status(new)["kind"] != "success",
                "before_status": status(old), "after_status": status(new),
                "before_solutions": old.get("solutions"), "after_solutions": new.get("solutions")})
        (out / "immediate-parent-comparison.jsonl").write_text(
            "".join(json.dumps(r) + "\n" for r in parent_changes), encoding="utf-8")
        summary["immediate_parent"] = {"outputs_sha256": digest(args.reference_outputs),
            "changed": sum(r["changed"] for r in parent_changes),
            "new_refusals": sum(r["new_refusal"] for r in parent_changes)}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    require(summary["decode_errors"] == 0, "output decoding failed")
    if args.expect_certified:
        require(not failed and summary["new_refusals"] == 0
                and summary.get("immediate_parent", {}).get("new_refusals", 0) == 0,
                "new refusal or native proof obligation needs investigation; no goldens may be blessed")


if __name__ == "__main__":
    main()
