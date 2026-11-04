#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Check #2381's split allocation against an already-built SMLib fork.

Run the normal runner independently first:
  nice -n 10 bin/difftest/difftest --unit blend.fillet_executive_p3 \
    --filter issue-2349 --digest-dir reports/difftest/issue-2349-fork-rust --digest-debug
Then invoke this script with --fork-source, --fork-lib-dir and --rust-debug.
The fork directory must be the clean a18bbe523b9c source; its existing
libsmlib_api.a / libsmlib_core.a are reused, never rebuilt. No shared pin changes.

All ten complete native records are retained in the report. The successful
comparison claimed here is ONLY split allocation at solver_geometry/corner_edges:
complete prefix topology, vertex state and both split vertices' five values.
It does not certify the rest of the redesigned construction. --pinned-shim can
add a full-record diagnostic against the old reference; any differences fail
that additional comparison rather than being hidden or given wider tolerances.
"""

import argparse
import ctypes
import hashlib
import json
from pathlib import Path
import re
import struct
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
FIXTURE = ROOT / 'fixtures/blend.fillet_executive_p3/issue_2349_closed_angle.json'
FORK = 'a18bbe523b9c7819f1a9da8c40a54c15f9936125'
CODES = {'corners': 2, 'first_topology': 3, 'first_vertex_geometry': 4,
         'solver_geometry': 7, 'corner_edges': 8}


def native(shim, index):
    """One call per process: upstream context state cannot leak between cases."""
    case = json.loads(FIXTURE.read_text())['cases'][index]['inputs']['real_solid']
    library = ctypes.CDLL(str(shim.resolve()))
    path = str((ROOT / case['source']['path']['$repo_file']).resolve()).encode()
    wire = [case['radius'], case['xsect'], case['continuity'],
            case['thumbweight'], len(case['edges']), *case['edges']]
    checkpoint = case.get('checkpoint', 'full')
    prefix = checkpoint != 'full'
    wire = [1.0, 1.0, 1.0, *wire] if prefix else [1.0, *([0.0] * 12), *wire]
    fn = (library.sc_blend_fillet_executive_p3_real_prefix_file if prefix else
          library.sc_blend_fillet_executive_p3_real_source)
    fn.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_double), ctypes.c_int64]
    if prefix:
        fn.argtypes += [ctypes.c_int32]
    fn.argtypes += [ctypes.POINTER(ctypes.c_double), ctypes.c_int64,
                   ctypes.POINTER(ctypes.c_int64)]
    fn.restype = ctypes.c_int32
    data = (ctypes.c_double * len(wire))(*wire)
    output = (ctypes.c_double * 1000000)()
    length = ctypes.c_int64()
    args = [path, data, len(wire)]
    if prefix:
        args += [CODES[checkpoint]]
    rc = fn(*args, output, len(output), ctypes.byref(length))
    if rc:
        raise RuntimeError(f'case {index}: native shim returned {rc}')
    return list(output[:length.value])


def allocation(values, checkpoint):
    cursor = 0

    def section():
        nonlocal cursor
        count = int(values[cursor])
        result = values[cursor + 1:cursor + 1 + count]
        cursor += 1 + count
        return result

    section()  # selected edge, retained in the complete raw record
    cursor += CODES[checkpoint] - 1  # statuses (6 or 7)
    section()  # before_direct_topology
    state = section()
    section()  # solver/corner parameters
    section()  # corner flags
    vertex_state = section()
    geometry = section()
    return checked_allocation(state, vertex_state, geometry)


def checked_allocation(state, vertex_state, geometry):
    # This operand has one corner, two original corner vertices and two split
    # vertices. Snapshot records seven integers per corner vertex, including its
    # owning-list label. -2 means present in the corner, absent from owning list.
    assert state[2] == 1 and state[7] == 4
    assert state[9 + 2 * 7 + 1] == -2 and state[9 + 3 * 7 + 1] == -2
    processed = {}
    cursor = 0
    for i in range(0, len(vertex_state), 4):
        corner, vertex, _, ready = vertex_state[i:i + 4]
        if ready:
            assert corner == 0
            processed[int(vertex)] = geometry[cursor:cursor + 5]
            cursor += 5
    assert cursor == len(geometry)
    split = [processed[2], processed[3]]
    assert all(len(v) == 5 and v[3] == 0.0 for v in split)
    return state, vertex_state, split


def rust_allocation(line):
    def field(key):
        match = re.search(r'"' + key + r'": (\[[^\]]*\])', line)
        if not match:
            raise ValueError(f'missing Rust field {key}')
        return json.loads(match.group(1))

    return checked_allocation(field('prefix_state'), field('all_corner_vertex_state'),
                              field('all_corner_vertices'))


def bits(value):
    if isinstance(value, list) or isinstance(value, tuple):
        return [bits(v) for v in value]
    return struct.pack('!d', float(value)).hex()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fork-source', type=Path, required=True)
    parser.add_argument('--fork-lib-dir', type=Path, required=True)
    parser.add_argument('--rust-debug', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--pinned-shim', type=Path)
    args = parser.parse_args()
    revision = subprocess.check_output(['git', '-C', str(args.fork_source),
                                        'rev-parse', 'HEAD'], text=True).strip()
    assert revision == FORK, revision
    assert not subprocess.check_output(['git', '-C', str(args.fork_source),
                                       'status', '--porcelain'], text=True).strip()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    source = (ROOT / 'harness/cpp/shim/blend_fillet_executive_p3.cpp').read_text()
    unrelated = '#include "blend_fillet_executive_p3_two_offset_chamfer.inc"'
    assert source.count(unrelated) == 1
    # The fork deleted the old chamfer FSG. Only that unrelated export is left
    # out; real_prefix_body, snapshot and every closed-round observation remain.
    translation_unit = args.out_dir / 'executive.cpp'
    translation_unit.write_text('#include <SmAttribute.h>\n' + source.replace(unrelated, ''))
    shim = (args.out_dir / 'fork-shim.so').resolve()
    libraries = [args.fork_lib_dir / f'libsmlib_{name}.a' for name in ['api', 'core']]
    command = ['nice', '-n', '10', 'g++', '-std=c++17', '-O2', '-ffp-contract=off',
               '-fno-fast-math', '-fPIC', '-shared', '-DNDEBUG', '-DNOMINMAX', '-DSM_STATIC',
               '-I', str(args.fork_source / 'source/SM_API/inc'),
               '-I', str(args.fork_source / 'source/SMLib/inc'),
               '-I', str(ROOT / 'harness/cpp/shim'), str(translation_unit),
               str(ROOT / 'harness/cpp/shim/sc_shim.cpp'), '-Wl,--no-undefined',
               *map(str, libraries), '-lpthread', '-o', str(shim)]
    subprocess.run(command, check=True)
    debug = dict(line.split('\t', 1) for line in args.rust_debug.read_text().splitlines())
    records = []
    allocation_passes = 0
    full_differences = 0
    for i, case in enumerate(json.loads(FIXTURE.read_text())['cases']):
        def run(shim_path):
            result = subprocess.run([sys.executable, __file__, '--worker',
                                     str(shim_path.resolve()), str(i)],
                                    check=True, capture_output=True, text=True)
            return json.loads(result.stdout)
        values = run(shim)
        record = {'name': case['name'], 'fork': values}
        checkpoint = case['inputs']['real_solid'].get('checkpoint', 'full')
        if checkpoint in ['solver_geometry', 'corner_edges']:
            fork = allocation(values, checkpoint)
            rust = rust_allocation(debug[case['name']])
            assert bits(fork) == bits(rust), case['name']
            allocation_passes += 1
            print(f'PASS allocation {case["name"]}: topology, vertex state, two split vertices bit-identical')
        if args.pinned_shim:
            pinned = run(args.pinned_shim)
            equal = bits(pinned) == bits(values)
            full_differences += not equal
            record.update(pinned=pinned, complete_record_equal=equal)
            print(f'{"PASS" if equal else "DIFF"} complete record {case["name"]}')
        records.append(record)
    assert allocation_passes == 2
    report = {'fork_commit': revision, 'compile_command': command,
              'libraries': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in libraries},
              'allocation_passes': allocation_passes, 'records': records}
    (args.out_dir / 'checkpoints.json').write_text(json.dumps(report, indent=2))
    print('2/2 allocation checkpoints match independent Rust outputs bit for bit')
    return bool(full_differences)


if __name__ == '__main__':
    if len(sys.argv) == 4 and sys.argv[1] == '--worker':
        print(json.dumps(native(Path(sys.argv[2]), int(sys.argv[3])), allow_nan=False))
    else:
        sys.exit(main())
