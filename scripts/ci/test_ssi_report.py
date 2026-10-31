"""SSI publication contract: completeness, truthful summaries and Actions provenance."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import ssi_report as s

REPOSITORY = "tangency-labs/solidcore"
SOURCE = "a" * 40


def report():
    cases = [dict(id=f"case{i:03}", **{e: {"verdict": "PASS"} for e in s.ENGINES},
                  parity=dict(pass_=True, mismatches=0)) for i in range(s.CASE_COUNT)]
    for case in cases:
        case["parity"]["pass"] = case["parity"].pop("pass_")
    return dict(corpus=dict(repo="https://github.com/jyviko/ssi-corpus", commit=s.corpus_pin(), glob=None),
                approx_tol="1e-6", acc_tol=10.0 * 1e-6, samples_per_curve=4001, cases=cases,
                summary={**{e: {"PASS": s.CASE_COUNT} for e in s.ENGINES},
                         "parity": dict(compared=119, identical=119, mismatched=0, same_verdict=119)})


def provenance():
    return dict(schema_version=1, repository=REPOSITORY, source_sha=SOURCE,
                upstream_commit="b" * 40, upstream_version="1.0.1", corpus_commit=s.corpus_pin(),
                ref="refs/heads/main", event="schedule", run_id=101,
                run_url=f"https://github.com/{REPOSITORY}/actions/runs/101",
                generated_at="2026-10-28T03:23:00+00:00")


def run():
    return dict(id=101, head_sha=SOURCE, event="schedule", head_branch="main", path=s.WORKFLOW,
                status="completed", conclusion="success", repository={"full_name": REPOSITORY},
                head_repository={"full_name": REPOSITORY})


def write_bundle(directory, data=None, metadata=None):
    directory.mkdir(parents=True, exist_ok=True)
    (directory / s.REPORT).write_text(json.dumps(data if data is not None else report()))
    (directory / s.PROVENANCE).write_text(json.dumps(metadata if metadata is not None else provenance()))


class ReportTests(unittest.TestCase):
    def test_full_report_and_bound_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            write_bundle(directory)
            data, metadata = s.validate_bundle(directory, REPOSITORY, run=run())
            self.assertEqual(len(data["cases"]), 119)
            self.assertEqual(metadata["source_sha"], SOURCE)

    def test_missing_and_incomplete_reports_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(s.ReportError):
                s.validate_bundle(Path(tmp), REPOSITORY)
        changes = [lambda r: r["cases"].pop(),
                   lambda r: r["cases"].append(copy.deepcopy(r["cases"][0])),
                   lambda r: r["cases"][1].update(id=r["cases"][0]["id"]),
                   lambda r: r["cases"][0].pop("solidcore-ssi"),
                   lambda r: r["cases"][0]["rust"].update(verdict="?"),
                   lambda r: r["cases"][0]["rust"].update(verdict=[]),
                   lambda r: r["summary"].pop("solidcore-ssi"),
                   lambda r: r["summary"]["parity"].update(same_verdict=0),
                   lambda r: r["summary"]["cpp"].update(PASS=118),
                   lambda r: r["corpus"].update(glob="c1_*"),
                   lambda r: r["corpus"].update(commit="c" * 40),
                   lambda r: r.update(approx_tol="default"),
                   lambda r: r.update(samples_per_curve=100)]
        for change in changes:
            data = report()
            change(data)
            with self.subTest(change=change), self.assertRaises(s.ReportError):
                s.validate_report(data, s.corpus_pin())

    def test_geometry_failures_and_parity_mismatches_remain_report_data(self):
        data = report()
        for engine in s.ENGINES:
            data["cases"][0][engine]["verdict"] = "INCOMPLETE"
            data["cases"][1][engine]["verdict"] = "ERROR"
            data["summary"][engine] = {"PASS": 117, "INCOMPLETE": 1, "ERROR": 1}
        s.validate_report(data, s.corpus_pin())
        data["cases"][2]["rust"]["verdict"] = "WRONG TOPOLOGY"
        data["cases"][2]["parity"].update({"pass": False, "mismatches": 1})
        data["summary"]["rust"] = {"PASS": 116, "INCOMPLETE": 1, "ERROR": 1, "WRONG TOPOLOGY": 1}
        data["summary"]["parity"].update(identical=118, mismatched=1, same_verdict=118)
        s.validate_report(data, s.corpus_pin())

    def test_crashes_and_deferrals_are_counted_without_inventing_parity(self):
        data = report()
        data["cases"][0]["rust"] = {"verdict": "DEFERRED", "deferred": {"unit": "a.b"}}
        data["cases"][0].pop("parity")
        for engine in ("cpp", "rust"):
            data["cases"][1][engine] = {"verdict": "CRASH"}
        data["cases"][1].pop("parity")
        data["summary"]["cpp"] = {"PASS": 118, "CRASH": 1}
        data["summary"]["rust"] = {"PASS": 117, "CRASH": 1, "DEFERRED": 1}
        data["summary"]["parity"] = dict(compared=117, identical=117, mismatched=0, same_verdict=118)
        s.validate_report(data, s.corpus_pin())

    def test_provenance_must_identify_the_downloaded_main_run(self):
        for key, value in [("repository", "someone/solidcore"), ("source_sha", "c" * 40),
                           ("corpus_commit", "d" * 40), ("ref", "refs/pull/1/merge"),
                           ("event", "pull_request"), ("event", []), ("generated_at", "not a date"),
                           ("generated_at", []),
                           ("generated_at", "2026-10-28T03:23:00"),
                           ("run_id", 102), ("run_url", "https://example.com")]:
            metadata = provenance()
            metadata[key] = value
            with tempfile.TemporaryDirectory() as tmp:
                directory = Path(tmp)
                write_bundle(directory, metadata=metadata)
                with self.subTest(key=key, value=value), self.assertRaises(s.ReportError):
                    s.validate_bundle(directory, REPOSITORY, run=run())

    def test_only_this_workflow_on_this_repositorys_main_is_trusted(self):
        self.assertTrue(s.trusted_run(run(), REPOSITORY))
        self.assertFalse(s.trusted_run(None, REPOSITORY))
        for change in [dict(path=".github/workflows/ci.yml"), dict(event="pull_request"),
                       dict(head_branch="other"), dict(conclusion="failure"),
                       dict(status="in_progress"), dict(event=[]), dict(repository={"full_name": "fork/solidcore"}),
                       dict(repository=None),
                       dict(head_repository={"full_name": "fork/solidcore"})]:
            with self.subTest(change=change):
                self.assertFalse(s.trusted_run({**run(), **change}, REPOSITORY))

    def test_fetch_skips_expired_artifacts_and_binds_download_to_run_sha(self):
        expired, corrupt, good = run(), {**run(), "id": 102}, {**run(), "id": 103}
        calls = []

        def api(path):
            if "workflows/" in path:
                return {"workflow_runs": [expired, corrupt, good]}
            return {"artifacts": [dict(name=s.ARTIFACT, expired="/101/" in path)]}

        def download(args, **_kwargs):
            calls.append(args)
            metadata = provenance()
            metadata.update(run_id=int(args[3]), run_url=f"https://github.com/{REPOSITORY}/actions/runs/{args[3]}")
            if args[3] == "102":
                metadata["source_sha"] = "c" * 40
            write_bundle(Path(args[-1]), metadata=metadata)

        with tempfile.TemporaryDirectory() as tmp, patch.object(s, "gh_json", side_effect=api), \
                patch.object(s.subprocess, "run", side_effect=download):
            directory = Path(tmp) / "reports"
            s.fetch(directory, REPOSITORY)
            _, metadata = s.validate_bundle(directory, REPOSITORY, run=good)
            self.assertEqual(metadata["run_id"], 103)
            self.assertEqual([call[3] for call in calls], ["102", "103"])

    def test_missing_artifacts_never_fall_back_to_committed_report(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(s, "gh_json", return_value={"workflow_runs": []}):
            directory = Path(tmp) / "ssi"
            (Path(tmp) / s.REPORT).write_text(json.dumps(report()))
            with self.assertRaises(s.ReportError):
                s.fetch(directory, REPOSITORY)
            self.assertFalse(directory.exists())

    def test_record_checks_actual_revisions_and_keeps_report_deterministic(self):
        pin = dict(line.split("=", 1) for line in (s.ROOT / "scripts/upstream-pin.env").read_text().splitlines()
                   if line.startswith("UPSTREAM_"))
        environment = dict(GITHUB_SHA=SOURCE, GITHUB_REPOSITORY=REPOSITORY, GITHUB_RUN_ID="101",
                           GITHUB_REF="refs/heads/main", GITHUB_EVENT_NAME="workflow_dispatch")
        with tempfile.TemporaryDirectory() as tmp, patch.dict(s.os.environ, environment):
            directory = Path(tmp)
            write_bundle(directory)
            before = (directory / s.REPORT).read_bytes()
            with patch.object(s, "git", side_effect=[SOURCE, "", pin["UPSTREAM_COMMIT"], s.corpus_pin()]):
                s.record(directory)
            _, metadata = s.validate_bundle(directory, REPOSITORY)
            self.assertEqual(metadata["upstream_commit"], pin["UPSTREAM_COMMIT"])
            self.assertEqual(metadata["event"], "workflow_dispatch")
            self.assertEqual((directory / s.REPORT).read_bytes(), before)
            with patch.object(s, "git", side_effect=[SOURCE, "crates/changed.rs"]), self.assertRaises(s.ReportError):
                s.record(directory)


if __name__ == "__main__":
    unittest.main()
