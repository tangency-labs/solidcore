"""parked_upstream.py: unit selection, heap qualification and the tracking issue.

No network, no build: reports are data, as in the workflow's `report` job.
"""
import io
import json
from contextlib import redirect_stdout
from pathlib import Path
import tempfile
import unittest

import parked_upstream as p


def write(root, rel, text):
    path = Path(root) / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


PARKED_ND = """[[mismatch]]
unit = "a.b"
case = "parked"
why = "C++ iterates a pointer-keyed hash map"
cpp_defect = "ND-1"
"""

PARKED_OTHER = """[[mismatch]]
unit = "c.d"
case = "plain"
why = "C++ reads freed memory"
cpp_defect = "UB-2"
"""

PLAIN = """[[mismatch]]
unit = "e.f"
case = "unparked"
why = "Rust diverges"
"""


def qualified(case, pads=p.REQUIRED_PADS, observed=None, padded=None, divergence=None):
    """A report's `heap_qualification` block with every pad observed by default."""
    obs = {
        "case": case,
        "observed": len(pads) if observed is None else observed,
        "padded": len(pads) if padded is None else padded,
    }
    if divergence is not None:
        obs["divergence"] = divergence
    return {"pads": list(pads), "cases": [obs]}


def adoption_rows(body):
    """Extract adoption-table data without pinning prose or cell spacing."""
    return {
        tuple(cell.strip().strip("`") for cell in line.split("|")[1:-1])
        for line in body.splitlines()
        if line.startswith("|")
    }


class ParkedUpstreamTest(unittest.TestCase):
    def run_main(self, *argv):
        out = io.StringIO()
        with redirect_stdout(out):
            code = p.main(list(argv))
        return code, out.getvalue()

    def issue_args(self, root, reports, out):
        return ("issue", "--root", root, "--report-dir", str(reports),
                "--upstream-ref", "v1.2", "--upstream-sha", "abc",
                "--title-out", str(Path(out) / "title"),
                "--body-out", str(Path(out) / "body"))

    def test_units_lists_only_units_with_a_cpp_defect(self):
        with tempfile.TemporaryDirectory() as root:
            write(root, p.KNOWN_DIR + "/a.b.toml", PARKED_ND)
            write(root, p.KNOWN_DIR + "/e.f.toml", PLAIN)
            code, out = self.run_main("units", "--root", root)
            self.assertEqual(code, 0)
            self.assertEqual(out.split(), ["a.b"])

    def test_non_nd_units_do_not_request_heap_padding(self):
        with tempfile.TemporaryDirectory() as root:
            write(root, p.KNOWN_DIR + "/a.b.toml", PARKED_ND)
            write(root, p.KNOWN_DIR + "/c.d.toml", PARKED_OTHER)
            self.assertEqual(
                self.run_main("pads", "--root", root, "--unit", "c.d"), (0, ""))

    def test_nd_pass_is_adoptable_only_with_complete_heap_qualification(self):
        with tempfile.TemporaryDirectory() as root, tempfile.TemporaryDirectory() as out:
            write(root, p.KNOWN_DIR + "/a.b.toml", PARKED_ND)
            write(root, p.KNOWN_DIR + "/c.d.toml", PARKED_OTHER)
            reports = Path(out) / "reports"
            write(reports, "a.b.json", json.dumps(
                {"unit": "a.b", "cases": [{"name": "parked", "outcome": "pass"}]}))
            write(reports, "c.d.json", json.dumps(
                {"unit": "c.d", "cases": [{"name": "plain", "outcome": "pass"}]}))

            # A non-ND defect keeps the plain pass rule; the ND pass alone proves nothing.
            found, missing, unqualified = p.agreeing(reports, p.parked_cases(Path(root)))
            self.assertEqual((found, missing), ([("c.d", "plain", "UB-2")], []))
            self.assertEqual(unqualified, [("a.b", "parked", "ND-1")])

            report = json.loads((reports / "a.b.json").read_text())
            report["heap_qualification"] = qualified("parked")
            write(reports, "a.b.json", json.dumps(report))
            found, _, unqualified = p.agreeing(reports, p.parked_cases(Path(root)))
            self.assertEqual(found, [("a.b", "parked", "ND-1"), ("c.d", "plain", "UB-2")])
            self.assertEqual(unqualified, [])

            # Evidence from a different unit cannot qualify a same-named case.
            report["unit"] = "other.unit"
            write(reports, "a.b.json", json.dumps(report))
            found, _, unqualified = p.agreeing(reports, p.parked_cases(Path(root)))
            self.assertEqual(found, [("c.d", "plain", "UB-2")])
            self.assertEqual(unqualified, [("a.b", "parked", "ND-1")])

    def test_incomplete_or_wrong_heap_qualification_is_not_adopted(self):
        case = "parked"
        bad = {
            "no pads block": None,
            "subset of pads": qualified(case, pads=p.REQUIRED_PADS[:-1]),
            "extra pad": qualified(case, pads=p.REQUIRED_PADS + (144,)),
            "reordered pads": qualified(case, pads=tuple(reversed(p.REQUIRED_PADS))),
            "rerun stopped early": qualified(case, observed=3),
            "in-process reruns": qualified(case, padded=0),
            "address-dependent divergence": qualified(case, divergence="pad=1: v: 1 vs 2"),
            "empty divergence record": qualified(case, divergence=""),
            "nonintegral observation count": qualified(case, observed=10.0),
            "boolean pad": qualified(case, pads=(True,) + p.REQUIRED_PADS[1:]),
            "observation for another case": qualified("other"),
            "no observation": {"pads": list(p.REQUIRED_PADS), "cases": []},
            "integer observations": {"pads": list(p.REQUIRED_PADS), "cases": 5},
            "boolean observations": {"pads": list(p.REQUIRED_PADS), "cases": True},
            "float observations": {"pads": list(p.REQUIRED_PADS), "cases": 2.5},
            "null observations": {"pads": list(p.REQUIRED_PADS), "cases": None},
            "object observations": {"pads": list(p.REQUIRED_PADS), "cases": {}},
            "string observations": {"pads": list(p.REQUIRED_PADS), "cases": "parked"},
            "block not an object": [list(p.REQUIRED_PADS)],
        }
        with tempfile.TemporaryDirectory() as root, tempfile.TemporaryDirectory() as out:
            write(root, p.KNOWN_DIR + "/a.b.toml", PARKED_ND)
            reports = Path(out) / "reports"
            for why, block in bad.items():
                report = {"unit": "a.b", "cases": [{"name": case, "outcome": "pass"}]}
                if block is not None:
                    report["heap_qualification"] = block
                write(reports, "a.b.json", json.dumps(report))
                with self.subTest(why):
                    self.assertFalse(p.heap_qualified(report, case), why)
                    found, _, unqualified = p.agreeing(
                        reports, p.parked_cases(Path(root)))
                    self.assertEqual(found, [], why)
                    self.assertEqual(unqualified, [("a.b", case, "ND-1")], why)

    def test_non_nd_pass_requires_matching_report_unit(self):
        with tempfile.TemporaryDirectory() as root, tempfile.TemporaryDirectory() as out:
            write(root, p.KNOWN_DIR + "/c.d.toml", PARKED_OTHER)
            reports = Path(out) / "reports"
            parked = p.parked_cases(Path(root))
            for report_unit in ("other.unit", None):
                report = {"cases": [{"name": "plain", "outcome": "pass"}]}
                if report_unit is not None:
                    report["unit"] = report_unit
                write(reports, "c.d.json", json.dumps(report))
                with self.subTest(report_unit=report_unit):
                    self.assertEqual(
                        p.agreeing(reports, parked),
                        ([], [], [("c.d", "plain", "UB-2")]))
                    self.assertEqual(
                        self.run_main(*self.issue_args(root, reports, out))[0], 0)
                    self.assertNotIn(
                        ("c.d", "plain", "UB-2"),
                        adoption_rows((Path(out) / "body").read_text()))

            # Correctly identified non-ND passes still need no heap qualification.
            write(reports, "c.d.json", json.dumps(
                {"unit": "c.d", "cases": [{"name": "plain", "outcome": "pass"}]}))
            self.assertEqual(
                p.agreeing(reports, parked), ([("c.d", "plain", "UB-2")], [], []))
            self.assertEqual(self.run_main(*self.issue_args(root, reports, out))[0], 0)
            self.assertIn(
                ("c.d", "plain", "UB-2"),
                adoption_rows((Path(out) / "body").read_text()))

    def test_report_with_nothing_selected_adopts_nothing_without_error(self):
        # `--run-parked` on a libc where the parked entry does not apply (the glibc-only
        # sweep entries on a libc++ host) selects no case: that is not evidence either.
        with tempfile.TemporaryDirectory() as root, tempfile.TemporaryDirectory() as out:
            write(root, p.KNOWN_DIR + "/a.b.toml", PARKED_ND)
            reports = Path(out) / "reports"
            write(reports, "a.b.json", json.dumps({"unit": "a.b", "cases": []}))
            self.assertEqual(
                p.agreeing(reports, p.parked_cases(Path(root))), ([], [], []))
            self.assertEqual(self.run_main(*self.issue_args(root, reports, out))[0], 0)
            self.assertNotIn(
                ("a.b", "parked", "ND-1"),
                adoption_rows((Path(out) / "body").read_text()))

    def test_issue_lists_only_qualified_cases_and_keeps_missing_reports_an_error(self):
        with tempfile.TemporaryDirectory() as root, tempfile.TemporaryDirectory() as out:
            write(root, p.KNOWN_DIR + "/a.b.toml", PARKED_ND)
            write(root, p.KNOWN_DIR + "/c.d.toml", PARKED_OTHER)
            reports = Path(out) / "reports"
            nd = {"unit": "a.b", "cases": [{"name": "parked", "outcome": "pass"}]}
            nd["heap_qualification"] = qualified("parked")
            write(reports, "a.b.json", json.dumps(nd))
            write(reports, "c.d.json", json.dumps(
                {"unit": "c.d", "cases": [{"name": "plain", "outcome": "pass"}]}))
            body = Path(out) / "body"
            args = self.issue_args(root, reports, out)
            self.assertEqual(self.run_main(*args)[0], 0)
            rows = adoption_rows(body.read_text())
            self.assertIn(("a.b", "parked", "ND-1"), rows)
            self.assertIn(("c.d", "plain", "UB-2"), rows)

            # The same ND pass, now without qualification: not in the adopt list.
            del nd["heap_qualification"]
            write(reports, "a.b.json", json.dumps(nd))
            self.assertEqual(self.run_main(*args)[0], 0)
            rows = adoption_rows(body.read_text())
            self.assertNotIn(("a.b", "parked", "ND-1"), rows)
            self.assertIn(("c.d", "plain", "UB-2"), rows)

            # No adoptable case at all: the issue must not prompt adoption.
            write(reports, "c.d.json", json.dumps(
                {"unit": "c.d", "cases": [{"name": "plain", "outcome": "cpp_defect"}]}))
            self.assertEqual(self.run_main(*args)[0], 0)
            rows = adoption_rows(body.read_text())
            self.assertNotIn(("a.b", "parked", "ND-1"), rows)
            self.assertNotIn(("c.d", "plain", "UB-2"), rows)

            # A missing unit report stays an error, not silence.
            (reports / "c.d.json").unlink()
            self.assertEqual(self.run_main(*args)[0], 2)


if __name__ == "__main__":
    unittest.main()
