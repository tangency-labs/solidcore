"""changes.sh: which pull-request diffs skip the build, test and difftest jobs (code=false)."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

SCRIPT = Path(__file__).with_name("changes.sh")

KNOWN = (
    '[[mismatch]]\nunit = "u"\ncase = "a"\nwhy = "w"\ncpp_defect = "ND-7"\n\n'
    '[[mismatch]]\nunit = "u"\ncase = "b"\nwhy = "w"\ncpp_defect = "WR-12-2"\n'
)


class Changes(unittest.TestCase):
    def setUp(self):
        self.dir = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.dir)
        self.git("init", "-q")
        self.write("crates/difftest/known_mismatches/u.toml", KNOWN)
        self.write("plan/upstream-issues/ND-7--parked.md", "#### ND-7. t\n- **Confirmed:** yes\n")
        self.write("plan/upstream-issues/WR-12-2--parked.md", "#### WR-12-2. t\n")
        self.write("plan/upstream-issues/WR-12--other.md", "#### WR-12. t\n")
        self.write("plan/upstream-issues/WR-8--other.md", "#### WR-8. t\n")
        self.write("src/lib.rs", "")
        self.base = self.commit()

    def git(self, *args):
        return subprocess.run(
            ["git", "-c", "user.name=t", "-c", "user.email=t@t", *args],
            cwd=self.dir,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()

    def write(self, rel, text):
        path = self.dir / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)

    def commit(self):
        self.git("add", "-A")
        self.git("commit", "-q", "--allow-empty", "-m", "c")
        return self.git("rev-parse", "HEAD")

    def code(self):
        self.commit()
        out = subprocess.run(
            ["bash", str(SCRIPT), "pull_request", self.base],
            cwd=self.dir,
            check=True,
            capture_output=True,
            text=True,
        ).stdout
        return dict(line.split("=", 1) for line in out.splitlines())["code"]

    def test_code_change_runs_the_jobs(self):
        self.write("src/lib.rs", "fn f() {}\n")
        self.assertEqual(self.code(), "true")

    def test_findings_and_unreferenced_records_skip_them(self):
        self.write("plan/findings/u--x.md", "x\n")
        self.write("plan/upstream-issues/WR-8--other.md", "#### WR-8. changed\n")
        self.write("plan/upstream-issues/WR-9--new.md", "#### WR-9. t\n")
        self.assertEqual(self.code(), "false")

    def test_an_addendum_of_a_cpp_defect_record_skips_them(self):
        # difftest reads the record itself, never its addenda.
        self.write("plan/upstream-issues/ND-7--addendum-2026-10-29-more.md", "Updates: ND-7\n")
        self.assertEqual(self.code(), "false")

    def test_editing_a_record_a_cpp_defect_names_runs_them(self):
        self.write("plan/upstream-issues/ND-7--parked.md", "#### ND-7. t\n- **Confirmed:** no\n")
        self.assertEqual(self.code(), "true")

    def test_renaming_or_deleting_that_record_runs_them(self):
        self.git("mv", "plan/upstream-issues/ND-7--parked.md", "plan/upstream-issues/ND-7--x.md")
        self.assertEqual(self.code(), "true")
        self.base = self.git("rev-parse", "HEAD")
        self.git("rm", "-q", "plan/upstream-issues/ND-7--x.md")
        self.assertEqual(self.code(), "true")

    def test_a_second_record_with_that_id_runs_them(self):
        self.write("plan/upstream-issues/ND-7--second.md", "#### ND-7. t\n")
        self.assertEqual(self.code(), "true")

    def test_ids_are_matched_whole(self):
        # WR-12-2 is named; WR-12 is not.
        self.write("plan/upstream-issues/WR-12--other.md", "#### WR-12. changed\n")
        self.assertEqual(self.code(), "false")
        self.base = self.git("rev-parse", "HEAD")
        self.write("plan/upstream-issues/WR-12-2--parked.md", "#### WR-12-2. changed\n")
        self.assertEqual(self.code(), "true")


if __name__ == "__main__":
    unittest.main()
