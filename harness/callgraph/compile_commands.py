#!/usr/bin/env python3
"""Describe licence-cleared C++ inputs for analysis, without a reference build.

The release definitions and include order come from the former reference CMake
configuration. Commands use Clang's JSON compilation database arguments format;
this module never executes them. Direct extraction consumes entries in memory.
Run this file to print the same metadata as JSON for inspection.
"""

import argparse
import json
from pathlib import Path


TREES = ("source/SMLib", "source/SM_API", "tests/SM_API_test", "tests/prog_test")
DEFINES = ("NDEBUG", "NOMINMAX", "SM_GRAPHICS_CALLBACKS=1", "_GLIBCXX_USE_CXX11_ABI=1")
# These tests include uncleared Examples sources; the former build excluded them too.
EXAMPLE_TESTS = {"TestSmBooleanBoundsExample.cpp", "TestSmSTEPDerivativesExample.cpp"}


def entries(root):
    """Return deterministic analysis commands for the four cleared source groups.

    A missing or empty group refuses before analysis rather than yielding a
    seemingly complete graph from an incomplete checkout. Upstream itself may
    be a symlink, but source-group directories and source files may not be.
    """
    root = Path(root).absolute()
    upstream = root / "upstream"
    commands = []
    for tree in TREES:
        directory = upstream / tree / "src"
        if (upstream / tree).is_symlink() or directory.is_symlink() or not directory.is_dir():
            raise ValueError(f"{directory}: missing or symlinked source directory")
        sources = []
        for path in sorted(directory.glob("*.cpp")):
            if tree == "tests/SM_API_test" and path.name in EXAMPLE_TESTS:
                continue
            if path.is_symlink():
                raise ValueError(f"{path}: source symlinks are not analysed")
            if path.is_file():
                sources.append(path)
        if not sources:
            raise ValueError(f"{directory}: no C++ sources")
        includes = ["source/SMLib"]
        defines = list(DEFINES)
        if tree in ("source/SM_API", "tests/SM_API_test"):
            includes.insert(0, "source/SM_API")
        if tree.startswith("tests/"):
            includes.insert(0, tree)
        if tree == "source/SM_API":
            defines.append("SM_API_EXPORTS")
        for source in sources:
            args = ["clang++", "-std=c++17", *["-D" + d for d in defines]]
            # Keep the standalone test root. The old database also listed this
            # TU without its main; both wrote the same per-TU output filename.
            if tree == "tests/SM_API_test" and source.name == "SM_API_test.cpp":
                args.append("-DSM_API_TEST_STANDALONE")
            for include in includes:
                args.extend(["-I", str(upstream / include / "inc")])
            args.extend(["-fsyntax-only", str(source)])
            commands.append({"directory": str(root), "file": str(source), "arguments": args})
    return sorted(commands, key=lambda entry: entry["file"])


def main():
    """Print inspectable metadata; do not invoke a compiler or write build files."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    try:
        commands = entries(args.root)
    except ValueError as exc:
        parser.error(str(exc))
    print(json.dumps(commands, indent=2))


if __name__ == "__main__":
    main()
