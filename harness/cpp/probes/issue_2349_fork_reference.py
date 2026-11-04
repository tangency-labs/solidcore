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

Decode saved complete records against independent Rust output (no C++ rebuild):
  python3 -B harness/cpp/probes/issue_2349_fork_reference.py --analyze-report \
    reports/difftest/issue-2349-fork-reference/checkpoints.json \
    reports/difftest/issue-2349-fork-rust/blend.fillet_executive_p3.debug \
    reports/difftest/issue-2349-named-differences.json
This mode exits unsuccessfully on any Rust/fork disagreement. The separately
observed fillet_point_step is outside the shim's record and is not certified here.
"""

import argparse
import ctypes
import hashlib
import json
import math
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


def decode_record(values, checkpoint):
    """Decode every observation using the existing Rust runner's wire layout."""
    cursor = 0
    fields = {}

    def take(count):
        nonlocal cursor
        if count < 0 or cursor + count > len(values):
            raise ValueError('short native record')
        result = values[cursor:cursor + count]
        cursor += count
        return result

    def section(name):
        count, = take(1)
        if not math.isfinite(count) or count < 0 or int(count) != count:
            raise ValueError('invalid section length')
        fields[name] = take(int(count))

    section('selected_edge_before')
    if checkpoint == 'full':
        fields['api_status_vertex_edge_face_counts'] = take(4)
        section('result_topology')
        section('result_geometry')
        fields['mass_status'] = take(1)
        fields['mass_properties'] = take(12 if fields['mass_status'] == [0] else 0)
    else:
        code = CODES[checkpoint]
        fields['prefix_statuses'] = take({2: 3, 3: 3, 4: 4, 7: 6, 8: 7}[code])
        for name in ['before_direct_topology', 'prefix_state',
                     'actual_solver_corner_parameters']:
            section(name)
        if code == 4:
            for name in ['calculated_vertex_state', 'calculated_vertices']:
                section(name)
        if code >= 5:
            for name in ['corner_calculation_flags', 'all_corner_vertex_state',
                         'all_corner_vertices']:
                section(name)
        if code >= 4:
            for phase in ['before', 'after']:
                section(f'offsets_{phase}_state')
                section(f'offsets_{phase}_distances')
        if code >= 6:
            section('corner_adjustment_flags')
        if code >= 7:
            section('generated_geometry_state')
            section('generated_geometry_samples')
        if code >= 8:
            section('corner_edge_geometry_state')
            section('corner_edge_geometry_samples')
    if cursor != len(values):
        raise ValueError('trailing native observations')
    return fields


def differences(left, right):
    """Exact discrete fields, relative 1e-9 values, no absolute floor."""
    result = []
    for field in sorted(left.keys() | right.keys()):
        a, b = left.get(field), right.get(field)
        if a is None or b is None or len(a) != len(b):
            result.append({'field': field, 'index': None, 'left': a, 'right': b})
            continue
        exact = (field.endswith(('_state', '_flags', '_statuses', '_topology'))
                 or field in ['prefix_state', 'before_direct_topology', 'mass_status',
                              'api_status_vertex_edge_face_counts'])
        for index, (x, y) in enumerate(zip(a, b)):
            equal = (math.isfinite(x) and math.isfinite(y)
                     and (x == y if exact else abs(x-y) <= 1e-9 * max(abs(x), abs(y))))
            if not equal:
                result.append({'field': field, 'index': index, 'left': x, 'right': y})
    return result


def analyze_report(report_path, rust_debug, output):
    """Analyze saved native records without rebuilding or rerunning C++."""
    report = json.loads(Path(report_path).read_text())
    if report['fork_commit'] != FORK:
        raise ValueError('unexpected fork revision in saved report')
    debug = dict(line.split('\t', 1) for line in Path(rust_debug).read_text().splitlines())
    cases = {case['name']: case for case in json.loads(FIXTURE.read_text())['cases']}
    summary = []
    for record in report['records']:
        name = record['name']
        checkpoint = cases[name]['inputs']['real_solid'].get('checkpoint', 'full')
        fork = decode_record(record['fork'], checkpoint)
        # The runner's separate fillet_point_step observation is not in this wire
        # record; every wire field is compared, with no other fields dropped.
        rust = {key: json.loads(value) for key, value in re.findall(
            r'"([^"]+)": (\[[^\]]*\])', debug[name]) if key != 'fillet_point_step'}
        delta = differences(rust, fork)
        entry = {'name': name, 'rust_vs_fork': delta}
        if 'pinned' in record:
            pinned = decode_record(record['pinned'], checkpoint)
            entry['pinned_vs_fork'] = differences(pinned, fork)
            entry['rust_vs_pinned'] = differences(rust, pinned)
        summary.append(entry)
        print(f'{name}: Rust/fork {len(delta)} differing observations')
        for d in delta[:3]:
            print(f'  {d["field"]}[{d["index"]}]: {d["left"]} vs {d["right"]}')
    Path(output).write_text(json.dumps(summary, indent=2))
    return any(entry['rust_vs_fork'] for entry in summary)


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
    def digest(path):
        with path.open('rb') as source:
            return hashlib.file_digest(source, 'sha256').hexdigest()

    report = {'fork_commit': revision, 'compile_command': command,
              'libraries': {str(p): digest(p) for p in libraries},
              'allocation_passes': allocation_passes, 'records': records}
    (args.out_dir / 'checkpoints.json').write_text(json.dumps(report, indent=2))
    print('2/2 allocation checkpoints match independent Rust outputs bit for bit')
    return bool(full_differences)


if __name__ == '__main__':
    if len(sys.argv) == 5 and sys.argv[1] == '--analyze-report':
        sys.exit(analyze_report(*sys.argv[2:]))
    elif len(sys.argv) == 4 and sys.argv[1] == '--worker':
        print(json.dumps(native(Path(sys.argv[2]), int(sys.argv[3])), allow_nan=False))
    else:
        sys.exit(main())
