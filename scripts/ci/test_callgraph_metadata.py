"""Check callgraph source selection without libclang, CMake or a C++ compiler."""

import contextlib
import importlib.util
import io
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch


CALLGRAPH = Path(__file__).resolve().parents[2] / "harness/callgraph"
TREES = ("source/SMLib", "source/SM_API", "tests/SM_API_test", "tests/prog_test")
SOURCES = (
    "source/SMLib/src/SmSecond.cpp",
    "source/SMLib/src/SmFirst.cpp",
    "source/SM_API/src/SmApi.cpp",
    "tests/SM_API_test/src/SM_API_test.cpp",
    "tests/SM_API_test/src/TestAllowed.cpp",
    "tests/prog_test/src/test.cpp",
)


class CallgraphMetadataTests(unittest.TestCase):
    """Drive the real extractor's metadata boundary with a recording parser."""

    def setUp(self):
        """Make a cleared-source fixture under a path containing spaces."""
        self.temp = tempfile.TemporaryDirectory(prefix="callgraph metadata ")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        for tree in TREES:
            (self.root / "upstream" / tree / "inc").mkdir(parents=True)
            (self.root / "upstream" / tree / "src").mkdir()
        for source in SOURCES:
            self.write_source(source)
        with patch.object(sys, "path", [str(CALLGRAPH), *sys.path]):
            spec = importlib.util.spec_from_file_location("callgraph_extract_test", CALLGRAPH / "extract.py")
            self.extract = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(self.extract)
        self.seen = []

    def write_source(self, name):
        """Write a harmless source fixture, never a real upstream file."""
        path = self.root / "upstream" / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("void fixture() {}\n")
        return path

    def run_extract(self, *args):
        """Record analysis commands in place of libclang and return any refusal."""
        seen = self.seen

        class RecordingPool:
            def __init__(self, jobs):
                self.jobs = jobs

            def __enter__(self):
                return self

            def __exit__(self, *exc):
                return False

            def imap_unordered(self, work, entries):
                for entry in entries:
                    seen.append(entry)
                    yield entry["file"], "metadata only"

        upstream = self.root / "upstream"
        allowed = [str(upstream / tree) + os.sep for tree in TREES]
        with patch.multiple(self.extract, ROOT=str(self.root), UPSTREAM=str(upstream), ALLOWED=allowed,
                            Pool=RecordingPool), patch.object(sys, "argv", ["extract.py", "--jobs", "1", *args]), \
                patch.object(self.extract, "CCDB", str(self.root / "bin/cpp/compile_commands.json"), create=True), \
                contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            try:
                self.extract.main()
            except (Exception, SystemExit) as exc:
                return str(exc)
        return None

    def test_direct_extraction_needs_no_retired_build_directory(self):
        """Fresh extraction reaches the parser with no old compilation database."""
        self.assertIsNone(self.run_extract())
        self.assertEqual({Path(e["file"]).relative_to(self.root / "upstream").as_posix() for e in self.seen},
                         set(SOURCES))
        self.assertFalse((self.root / "bin/cpp").exists())

    def test_only_four_cleared_nonrecursive_source_trees_are_selected(self):
        """Uncleared, nested, non-C++ and example-dependent sources stay out."""
        for name in ("source/SMLib_extra/src/no.cpp", "unrelated/no.cpp", "source/SMLib/src/nested/no.cpp",
                     "source/SMLib/src/no.c", "tests/SM_API_test/src/TestSmBooleanBoundsExample.cpp",
                     "tests/SM_API_test/src/TestSmSTEPDerivativesExample.cpp"):
            self.write_source(name)
        self.assertIsNone(self.run_extract())
        self.assertEqual(len(self.seen), len(SOURCES))

    def test_release_definitions_and_target_specific_definitions(self):
        """Carry the old parser configuration, not parallel/debug code paths."""
        self.assertIsNone(self.run_extract())
        for entry in self.seen:
            args = entry["arguments"]
            for flag in ("-std=c++17", "-DNDEBUG", "-DNOMINMAX", "-DSM_GRAPHICS_CALLBACKS=1",
                         "-D_GLIBCXX_USE_CXX11_ABI=1"):
                self.assertIn(flag, args)
            self.assertEqual("-DSM_API_EXPORTS" in args, "/source/SM_API/src/" in entry["file"])
            self.assertEqual("-DSM_API_TEST_STANDALONE" in args, Path(entry["file"]).name == "SM_API_test.cpp")
            self.assertFalse(any("SM_USE_TBB" in arg or "SM_USE_GLOBAL_CACHE" in arg for arg in args))

    def test_include_precedence_and_paths_with_spaces(self):
        """Keep each unescaped include operand intact and in target order."""
        self.assertIsNone(self.run_extract())
        for entry in self.seen:
            source = Path(entry["file"]).relative_to(self.root / "upstream")
            owner = source.parts[:2]
            trees = ["source/SMLib"]
            if owner != ("source", "SMLib") and owner != ("tests", "prog_test"):
                trees.insert(0, "source/SM_API")
            if owner[0] == "tests":
                trees.insert(0, "/".join(owner))
            expected = [str(self.root / "upstream" / tree / "inc") for tree in trees]
            args = entry["arguments"]
            self.assertEqual([args[i + 1] for i, a in enumerate(args) if a == "-I"], expected)
            self.assertEqual(entry["directory"], str(self.root))
            self.assertIn(entry["file"], args)
            parsed = self.extract.args_for(entry)
            self.assertTrue(all(path in parsed for path in expected))

    def test_selection_is_repeatable_and_filter_is_preserved(self):
        """Equal-size inputs have stable order and the existing filter still works."""
        self.assertIsNone(self.run_extract())
        first = list(self.seen)
        self.seen.clear()
        self.assertIsNone(self.run_extract())
        self.assertEqual(self.seen, first)
        self.assertEqual([e["file"] for e in self.seen], sorted(e["file"] for e in self.seen))
        self.seen.clear()
        self.assertIsNone(self.run_extract("--filter", "SmFirst.cpp"))
        self.assertEqual(len(self.seen), 1)
        self.assertTrue(self.seen[0]["file"].endswith("/SmFirst.cpp"))

    def test_missing_source_tree_refuses_before_parser(self):
        """An incomplete checkout cannot silently produce a partial graph."""
        (self.root / "upstream/tests/prog_test/src/test.cpp").unlink()
        (self.root / "upstream/tests/prog_test/src").rmdir()
        error = self.run_extract()
        self.assertIsNotNone(error)
        self.assertIn("source directory", error)
        self.assertEqual(self.seen, [])

    def test_empty_source_tree_refuses_before_parser(self):
        """An empty source group is distinguished from a valid empty filter."""
        (self.root / "upstream/tests/prog_test/src/test.cpp").unlink()
        error = self.run_extract()
        self.assertIsNotNone(error)
        self.assertIn("no C++ sources", error)
        self.assertEqual(self.seen, [])


if __name__ == "__main__":
    unittest.main()
