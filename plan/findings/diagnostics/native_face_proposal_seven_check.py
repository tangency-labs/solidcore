"""Audit frozen native face outputs; optionally rerun the standalone native evaluator.

This checks fields, admissibility, ideal closed forms and native scalar identity.
It does not replace the kernel's unchanged own-enclosure/Bernstein certificate.
No golden file is written.
"""
import argparse
import json
import math
import subprocess
from decimal import Decimal, localcontext
from pathlib import Path


def digest(values):
    total = absolute = weighted = 0.0
    for i, value in enumerate(values):
        total += value
        absolute += abs(value)
        weighted += (i + 1) * value
    return [total, absolute, weighted]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--records', type=Path, default=Path(__file__).with_name('native_face_proposal_seven_records.jsonl'))
    parser.add_argument('--native-exe', type=Path)
    args = parser.parse_args()
    rows = [json.loads(line) for line in args.records.read_text(encoding='utf-8-sig').splitlines()]
    assert len(rows) == 7 and len({row['case'] for row in rows}) == 7
    for row in rows:
        before, after = row['before_static'], row['after_static']
        assert len(before) == len(after) == 19
        assert before[:11] == after[:11]
        assert after[5:11] == [1, 1, 2, 2, 0, -2]
        assert before[14:] == [1, 1, 0, 0, 0]
        assert after[14:] == [0, 0, 0, 0, 0]
        inp = row['input']
        assert inp['first'][0] == [.2, .3, 2] and inp['entry'] == 5
        assert inp['operation'] == 15 and inp['request'] == 0
        kind = inp['faces'][0]['kind']
        with localcontext() as ctx:
            ctx.prec = 80
            d = Decimal
            if kind == 'cylinder':
                ideal = (d('.7')**2 + d('.8')**2).sqrt()
                assert 0 <= after[12] <= math.pi and 0 <= after[13] <= 1
            elif kind == 'cone':
                ideal = (d('.3')**2 + (d('.3') + d('1.5') / d(3).sqrt())**2).sqrt()
            elif kind == 'sphere':
                ideal = d('.5') + (d('.2')**2 + d('.3')**2).sqrt() if 'band' in row['case'] else (d('.7')**2 + d('.3')**2).sqrt()
            else:
                assert kind == 'torus'
                ideal = d('.25') + (d('.2')**2 + d('1.3')**2).sqrt()
            # An independent ideal-geometry diagnostic, not the acceptance slack.
            assert abs(Decimal.from_float(after[11]) - ideal) <= Decimal.from_float(2 * math.ulp(after[11]))
        if kind != 'cylinder':
            assert all(0 <= t <= 1 for t in after[12:14])
        if kind == 'torus':
            assert inp['faces'][0]['center'] == [0, 0, -2]
            assert inp['faces'][0]['x'] == [0, 1, 0] and inp['faces'][0]['y'] == [0, 0, 1]
            assert inp['faces'][0]['major_radius'] == 1 and inp['faces'][0]['minor_radius'] == .25
            assert inp['faces'][0]['minor'] == [0, 360]
            for side in ('before', 'after'):
                values = row[side + '_static']
                native = row[side + '_native']
                if args.native_exe:
                    fresh = json.loads(subprocess.check_output([str(args.native_exe.resolve()), str(inp['faces'][0]['major'][1]), str(values[12]), str(values[13])]))
                    assert fresh == native
                assert native['domain'] == [0, 0, 1, 1]
                assert native['value'] == values[11]
                x, y, _ = native['point']
                assert math.sqrt((x - .2)**2 + (y - .3)**2) == values[11]
        else:
            assert before[11:14] == after[11:14]
        expected = row['expected_at_0a']['values']['static']
        assert expected['n'] == 19
        assert expected['blocks'][0] != digest(after)
        print(row['case'], 'value=', after[11], 'uv=', after[12:14], 'delta=', [b-a for a,b in zip(digest(before),digest(after))])
    print('PASS: 7 typed native solutions; 4 torus evaluations; no expected output written')


if __name__ == '__main__':
    main()
