"""SSI page tests exercise measured data, independent engines and absent artifacts."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import ssi

REPOSITORY = "tangency-labs/solidcore"
PIN = "b" * 40


def report():
    return dict(corpus={"repo": "https://github.com/jyviko/ssi-corpus", "commit": "c" * 40},
                approx_tol="1e-6", acc_tol=1e-5, samples_per_curve=4001,
                summary={**{e: {"PASS": 1} for e in ssi.ssi_report.ENGINES},
                         "parity": {"compared": 1, "identical": 1, "mismatched": 0}},
                cases=[dict(id="c5_cyl_quilt_seam", **{e: {"verdict": "PASS"} for e in ssi.ssi_report.ENGINES},
                            surfaces=["cylinder", "cylinder"], parity={"pass": True, "mismatches": 0})])


def provenance():
    return dict(source_sha="a" * 40, upstream_commit=PIN, upstream_version="1.0.1",
                generated_at="2026-10-28T03:23:00+00:00",
                run_url=f"https://github.com/{REPOSITORY}/actions/runs/101")


class SsiSiteTests(unittest.TestCase):
    def test_absent_report_is_unavailable_instead_of_zero_scores(self):
        for render in (ssi.body, ssi.overview):
            text = render(None, None, REPOSITORY, PIN)
            self.assertIn("No SSI report available", text)
            self.assertNotIn("0 compared", text)
            self.assertNotIn("<table", text)

    def test_four_engines_keep_their_own_scores_and_only_cpp_rust_parity(self):
        data = report()
        data["cases"][0]["ssi"] = {"verdict": "INCOMPLETE", "flags": ["missing branch"]}
        data["summary"]["ssi"] = {"INCOMPLETE": 1}
        text = ssi.body(data, provenance(), REPOSITORY, PIN)
        for label in ("C++", "Rust", "ssi", "solidcore-ssi"):
            self.assertIn(f'>{label}</th>', text)
        self.assertIn('data-ssi="INCOMPLETE"', text)
        self.assertIn('name="solidcore-ssi"', text)
        self.assertIn("ssi: missing branch", text)
        self.assertIn("C++/Rust parity", text)
        self.assertIn("<td>agree</td>", text)
        self.assertNotIn("waits on", text)

    def test_renderer_can_show_two_engines_without_inventing_other_scores(self):
        data = report()
        for engine in ("ssi", "solidcore-ssi"):
            data["summary"].pop(engine)
            data["cases"][0].pop(engine)
        text = ssi.body(data, provenance(), REPOSITORY, PIN)
        self.assertNotIn('name="solidcore-ssi"', text)
        self.assertIn('name="cpp"', text)
        self.assertIn('name="rust"', text)

    def test_measured_revision_and_upstream_are_distinct_from_site_revision(self):
        text = ssi.body(report(), provenance(), REPOSITORY, "d" * 40)
        self.assertIn("2026-10-28 03:23 UTC", text)
        self.assertIn("/commit/" + "a" * 40, text)
        self.assertIn("SMLib v1.0.1", text)
        self.assertIn("different SMLib pin", text)
        self.assertIn('href="reports/ssi-corpus.json"', text)
        self.assertIn("/actions/runs/101", text)

    def test_genuine_deferral_and_engine_text_are_escaped(self):
        data = report()
        data["cases"][0]["rust"] = dict(verdict="DEFERRED",
            deferred={"unit": "intersect.surface_drop_lib_p3", "symbol": "<DropCurve>"})
        data["cases"][0]["ssi"] = dict(verdict="ERROR", detail="<script>bad</script>")
        data["cases"][0].pop("parity")
        text = ssi.body(data, provenance(), REPOSITORY, PIN)
        self.assertIn("waits on", text)
        self.assertIn("units/intersect.surface_drop_lib_p3.html", text)
        self.assertIn("&lt;DropCurve&gt;", text)
        self.assertIn("&lt;script&gt;bad&lt;/script&gt;", text)
        self.assertNotIn("<script>", text)

    def test_loader_does_not_read_a_legacy_or_invalid_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "reports").mkdir()
            (root / "reports/ssi-corpus.json").write_text('{"cases": []}')
            with patch.object(ssi.ssi_report, "corpus_pin", return_value="c" * 40):
                self.assertEqual(ssi.load(root, REPOSITORY), (None, None))
                (root / "reports/ssi").mkdir()
                (root / "reports/ssi/ssi-corpus.json").write_text("invalid JSON")
                self.assertEqual(ssi.load(root, REPOSITORY), (None, None))


if __name__ == "__main__":
    unittest.main()
