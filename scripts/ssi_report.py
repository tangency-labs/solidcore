#!/usr/bin/env python3
"""Validate, record and fetch the complete nightly SSI report (standard library only)."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
ENGINES = ("cpp", "rust", "ssi", "solidcore-ssi")
LABELS = {"cpp": "C++", "rust": "Rust", "ssi": "ssi", "solidcore-ssi": "solidcore-ssi"}
VERDICTS = {"PASS", "INACCURATE", "WRONG TOPOLOGY", "INCOMPLETE", "ERROR", "CRASH", "N/A", "DEFERRED"}
WORKFLOW = ".github/workflows/ssi-corpus.yml"
ARTIFACT = "ssi-corpus-report"
REPORT = "ssi-corpus.json"
PROVENANCE = "ssi-provenance.json"
CASE_COUNT = 119


class ReportError(ValueError):
    """An unavailable, incomplete or untrusted report."""


def require(condition, message):
    if not condition:
        raise ReportError(message)


def read_json(path):
    try:
        value = json.loads(Path(path).read_text())
        require(isinstance(value, dict), f"{path}: expected an object")
        return value
    except (OSError, ValueError) as error:
        raise ReportError(f"{path}: {error}") from error


def corpus_pin(root=ROOT):
    text = (root / "scripts/fetch-ssi-corpus.sh").read_text()
    match = re.search(r'^COMMIT="([0-9a-f]{40})"$', text, re.MULTILINE)
    require(match is not None, "missing corpus pin in fetch-ssi-corpus.sh")
    return match.group(1)


def validate_report(report, pin):
    require(report.get("corpus") == {"repo": "https://github.com/jyviko/ssi-corpus",
                                    "commit": pin, "glob": None}, "wrong corpus pin or filtered report")
    require(report.get("approx_tol") == "1e-6" and report.get("acc_tol") == 10.0 * 1e-6
            and report.get("samples_per_curve") == 4001, "wrong SSI sampling or tolerances")
    cases, summary = report.get("cases"), report.get("summary")
    require(isinstance(cases, list) and len(cases) == CASE_COUNT, "expected all 119 SSI cases")
    require(isinstance(summary, dict) and set(summary) == {*ENGINES, "parity"},
            "expected all four engines and C++/Rust parity")
    seen, counts = set(), {engine: Counter() for engine in ENGINES}
    parity = dict(compared=0, identical=0, mismatched=0, same_verdict=0)
    for case in cases:
        require(isinstance(case, dict), "invalid SSI case")
        case_id = case.get("id")
        require(isinstance(case_id, str) and case_id and case_id not in seen, "missing or duplicate case id")
        seen.add(case_id)
        for engine in ENGINES:
            grade = case.get(engine)
            require(isinstance(grade, dict) and isinstance(grade.get("verdict"), str)
                    and grade["verdict"] in VERDICTS,
                    f"{case_id}: missing or invalid {engine} grade")
            counts[engine][grade["verdict"]] += 1
        cpp, rust = (case[engine]["verdict"] for engine in ("cpp", "rust"))
        parity["same_verdict"] += int(rust != "DEFERRED" and cpp == rust)
        comparable = cpp not in {"CRASH", "DEFERRED"} and rust not in {"CRASH", "DEFERRED"}
        result = case.get("parity")
        require((result is not None) == comparable, f"{case_id}: missing or unexpected parity")
        if comparable:
            require(isinstance(result, dict) and type(result.get("pass")) is bool
                    and type(result.get("mismatches")) is int and result["mismatches"] >= 0,
                    f"{case_id}: invalid parity result")
            require(result["pass"] == (result["mismatches"] == 0), f"{case_id}: inconsistent parity result")
            parity["compared"] += 1
            parity["identical"] += int(result["pass"])
            parity["mismatched"] += int(not result["pass"])
    for engine in ENGINES:
        actual = summary[engine]
        require(isinstance(actual, dict) and all(k in VERDICTS and type(v) is int and v >= 0
                                               for k, v in actual.items()), f"invalid {engine} summary")
        require({k: v for k, v in actual.items() if v} == dict(counts[engine]),
                f"{engine} summary disagrees with cases")
    require(summary["parity"] == parity, "parity summary disagrees with cases")
    return report


def trusted_run(run, repository):
    return (isinstance(run, dict) and isinstance(run.get("repository"), dict)
            and isinstance(run.get("head_repository"), dict)
            and run.get("path") == WORKFLOW and run.get("head_branch") == "main"
            and isinstance(run.get("event"), str) and run["event"] in {"schedule", "workflow_dispatch"}
            and run.get("status") == "completed" and run.get("conclusion") == "success"
            and run.get("repository", {}).get("full_name") == repository
            and run.get("head_repository", {}).get("full_name") == repository)


def validate_bundle(directory, repository, root=ROOT, run=None):
    report = validate_report(read_json(directory / REPORT), corpus_pin(root))
    provenance = read_json(directory / PROVENANCE)
    require(provenance.get("schema_version") == 1 and provenance.get("repository") == repository,
            "wrong SSI provenance schema or repository")
    for key in ("source_sha", "upstream_commit", "corpus_commit"):
        require(isinstance(provenance.get(key), str) and re.fullmatch(r"[0-9a-f]{40}", provenance[key]),
                f"invalid provenance {key}")
    require(provenance["corpus_commit"] == report["corpus"]["commit"], "provenance corpus pin disagrees")
    require(provenance.get("ref") == "refs/heads/main"
            and isinstance(provenance.get("event"), str)
            and provenance["event"] in {"schedule", "workflow_dispatch"}, "report is not from main")
    run_id = provenance.get("run_id")
    require(type(run_id) is int and run_id > 0, "invalid provenance run id")
    require(provenance.get("run_url") == f"https://github.com/{repository}/actions/runs/{run_id}",
            "invalid provenance run URL")
    require(isinstance(provenance.get("upstream_version"), str) and provenance["upstream_version"],
            "missing upstream version")
    require(isinstance(provenance.get("generated_at"), str), "missing SSI report date")
    try:
        date = datetime.fromisoformat(provenance["generated_at"].replace("Z", "+00:00"))
        require(date.utcoffset() is not None, "SSI report date needs a timezone")
    except (KeyError, TypeError, ValueError) as error:
        raise ReportError("invalid SSI report date") from error
    if run is not None:
        require(trusted_run(run, repository) and run_id == run.get("id")
                and provenance["source_sha"] == run.get("head_sha")
                and provenance["event"] == run.get("event"), "provenance disagrees with Actions run")
    return report, provenance


def git(root, *args):
    return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()


def record(directory, root=ROOT):
    report = validate_report(read_json(directory / REPORT), corpus_pin(root))
    pin = dict(line.split("=", 1) for line in (root / "scripts/upstream-pin.env").read_text().splitlines()
               if line.startswith("UPSTREAM_"))
    source = git(root, "rev-parse", "HEAD")
    require(source == os.environ["GITHUB_SHA"], "checkout differs from workflow revision")
    require(not git(root, "diff", "HEAD", "--name-only"), "benchmark source has uncommitted changes")
    require(git(root / "upstream", "rev-parse", "HEAD") == pin["UPSTREAM_COMMIT"], "wrong upstream checkout")
    require(git(root / "corpus/ssi", "rev-parse", "HEAD") == report["corpus"]["commit"], "wrong corpus checkout")
    repository, run_id = os.environ["GITHUB_REPOSITORY"], int(os.environ["GITHUB_RUN_ID"])
    provenance = dict(schema_version=1, repository=repository, source_sha=source,
                      upstream_commit=pin["UPSTREAM_COMMIT"], upstream_version=pin["UPSTREAM_VERSION"],
                      corpus_commit=report["corpus"]["commit"], ref=os.environ["GITHUB_REF"],
                      event=os.environ["GITHUB_EVENT_NAME"], run_id=run_id,
                      run_url=f"https://github.com/{repository}/actions/runs/{run_id}",
                      generated_at=datetime.now(timezone.utc).isoformat())
    (directory / PROVENANCE).write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    validate_bundle(directory, repository, root)


def gh_json(path):
    return json.loads(subprocess.check_output(["gh", "api", path], text=True))


def fetch(directory, repository, root=ROOT):
    runs = gh_json(f"repos/{repository}/actions/workflows/ssi-corpus.yml/runs?branch=main&per_page=30")
    for run in runs.get("workflow_runs", []):
        if not trusted_run(run, repository):
            continue
        artifacts = gh_json(f"repos/{repository}/actions/runs/{run['id']}/artifacts?name={ARTIFACT}")
        if not any(a.get("name") == ARTIFACT and not a.get("expired", True)
                   for a in artifacts.get("artifacts", [])):
            continue
        try:
            with tempfile.TemporaryDirectory(prefix="solidcore-ssi-report-") as temporary:
                candidate = Path(temporary)
                subprocess.run(["gh", "run", "download", str(run["id"]), "--repo", repository,
                                "--name", ARTIFACT, "--dir", str(candidate)], check=True,
                               capture_output=True, text=True)
                validate_bundle(candidate, repository, root, run)
                directory.mkdir(parents=True, exist_ok=True)
                for name in (REPORT, PROVENANCE):
                    shutil.copyfile(candidate / name, directory / name)
            print(f"fetch-ssi-report: {ARTIFACT} from main run {run['id']} -> {directory}")
            return
        except (ReportError, subprocess.CalledProcessError) as error:
            print(f"fetch-ssi-report: skipping run {run['id']}: {error}")
    raise ReportError("no valid SSI report artifact from main; no committed-report fallback")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("record", "validate", "fetch"))
    parser.add_argument("directory", nargs="?", type=Path, default=ROOT / "reports/ssi")
    parser.add_argument("--repository", default=os.environ.get("GITHUB_REPOSITORY", "tangency-labs/solidcore"))
    args = parser.parse_args()
    try:
        if args.command == "fetch":
            fetch(args.directory, args.repository)
        elif args.command == "record":
            record(args.directory)
        else:
            validate_bundle(args.directory, args.repository)
    except (ValueError, OSError, KeyError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"ssi-report: {error}\n")


if __name__ == "__main__":
    main()
