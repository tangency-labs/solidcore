"""Original Decimal power-series oracle for two nonrational, degree-one-v cases.

Definition: tensor B-spline basis, then S+d*(Su cross Sv)/|Su cross Sv|.
Basis definition reference (not an implementation dependency):
https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.BSpline.html
Only input JSON is consumed. No Rust/C++ output is used to calculate the oracle.
"""
from decimal import Decimal as D, localcontext
from functools import lru_cache
from math import factorial
from pathlib import Path
import json
import argparse
from fractions import Fraction
ORDER = 4

def dec(x):
    return D.from_float(x) if isinstance(x, float) else D(x)

def const(x):
    return [dec(x)] + [D(0)] * ORDER

def add(a, b):
    return [x + y for (x, y) in zip(a, b)]

def scale(a, s):
    return [x * s for x in a]

def sub(a, b):
    return add(a, scale(b, D(-1)))

def mul(a, b):
    return [sum((a[j] * b[i - j] for j in range(i + 1)), D(0)) for i in range(ORDER + 1)]

def deriv(a):
    return [a[i + 1] * (i + 1) for i in range(ORDER)] + [D(0)]

def sqrt(a):
    b = [a[0].sqrt()] + [D(0)] * ORDER
    for i in range(1, ORDER + 1):
        b[i] = (a[i] - sum((b[j] * b[i - j] for j in range(1, i)), D(0))) / (2 * b[0])
    return b

def div(a, b):
    c = [D(0)] * (ORDER + 1)
    for i in range(ORDER + 1):
        c[i] = (a[i] - sum((b[j] * c[i - j] for j in range(1, i + 1)), D(0))) / b[0]
    return c

def cross(a, b):
    return [sub(mul(a[1], b[2]), mul(a[2], b[1])), sub(mul(a[2], b[0]), mul(a[0], b[2])), sub(mul(a[0], b[1]), mul(a[1], b[0]))]

def evaluate(inputs, precision):
    with localcontext() as ctx:
        ctx.prec = precision
        off = inputs['off']
        b = off['base']['bspline']
        step = inputs['steps'][0]
        assert 1 <= b['p'] <= 3 and b['q'] == 1 and (step['hv'] == 0) and step['upper']
        assert 0 <= step['hu'] <= 3
        assert all((p[3] == -float.fromhex('0x1.fffffffffffffp+1023') for row in b['net'] for p in row))
        U = [dec(x) for x in b['u']]
        V = [dec(x) for x in b['v']]
        (u, v) = map(dec, step['uv'])
        P = [[[dec(c) for c in p[:3]] for p in row] for row in b['net']]
        assert U[b['p']] < u < U[len(P)] and V[1] < v < V[len(P[0])]
        j = next((j for j in range(1, len(V) - 2) if V[j] <= v < V[j + 1]))
        alpha = (v - V[j]) / (V[j + 1] - V[j])
        points = [[(1 - alpha) * row[j - 1][c] + alpha * row[j][c] for c in range(3)] for row in P]
        vders = [[(row[j][c] - row[j - 1][c]) / (V[j + 1] - V[j]) for c in range(3)] for row in P]
        x = [u, D(1)] + [D(0)] * (ORDER - 1)

        @lru_cache(None)
        def basis(i, p):
            if p == 0:
                return tuple(const(int(U[i] <= u < U[i + 1])))
            out = const(0)
            if U[i + p] != U[i]:
                out = add(out, scale(mul(sub(x, const(U[i])), basis(i, p - 1)), 1 / (U[i + p] - U[i])))
            if U[i + p + 1] != U[i + 1]:
                out = add(out, scale(mul(sub(const(U[i + p + 1]), x), basis(i + 1, p - 1)), 1 / (U[i + p + 1] - U[i + 1])))
            return tuple(out)
        B = [basis(i, b['p']) for i in range(len(P))]
        S = [sum_series([scale(B[i], points[i][c]) for i in range(len(P))]) for c in range(3)]
        Sv = [sum_series([scale(B[i], vders[i][c]) for i in range(len(P))]) for c in range(3)]
        n = cross([deriv(s) for s in S], Sv)
        length = sqrt(sum_series([mul(c, c) for c in n]))
        out = [add(S[c], scale(div(n[c], length), dec(off['d']))) for c in range(3)]
        return {'precision': precision, 'normal_length': str(length[0]), 'derivatives': {str(i): [str(out[c][i] * factorial(i)) for c in range(3)] for i in range(step['hu'] + 1)}}

def sum_series(items):
    out = const(0)
    for a in items:
        out = add(out, a)
    return out

def analytic_checks():
    sentinel = -float.fromhex('0x1.fffffffffffffp+1023')

    def inp(p, knots, xy):
        return {'off': {'d': 0.75, 'base': {'bspline': {'p': p, 'q': 1, 'u': knots, 'v': [0, 0, 1, 1], 'net': [[[x, y, z, sentinel] for z in [0, 1]] for (x, y) in xy]}}}, 'steps': [{'hu': 2, 'hv': 0, 'upper': True, 'uv': [0.25, 0.5]}]}
    plane = evaluate(inp(1, [0, 0, 1, 1], [(0, 0), (1, 0)]), 140)['derivatives']
    assert list(map(D, plane['0'])) == [D('0.25'), D('-0.75'), D('0.5')]
    assert list(map(D, plane['1'])) == [D(1), D(0), D(0)]
    assert list(map(D, plane['2'])) == [D(0)] * 3
    curve = evaluate(inp(2, [0, 0, 0, 1, 1, 1], [(0, 0), (0.5, 0), (1, 1)]), 140)['derivatives']
    with localcontext() as ctx:
        ctx.prec = 140
        u = D('0.25')
        v = D('0.5')
        d = D('0.75')
        g = 1 + 4 * u * u
        root = g.sqrt()
        expected = {0: [u + d * 2 * u / root, u * u - d / root, v], 1: [1 + d * 2 / (g * root), 2 * u + d * 4 * u / (g * root), D(0)], 2: [-24 * d * u / (g * g * root), 2 + d * (4 - 32 * u * u) / (g * g * root), D(0)]}
        for (order, values) in expected.items():
            assert max((abs(D(a) - b) for (a, b) in zip(curve[str(order)], values))) < D('1e-125')
    return {'plane': 'exact', 'parabolic_cylinder': 'point,first,second derivatives agree with closed form to absolute1e-125 at140digits'}

def exact_mixed_derivative(inputs):
    """Compute the bilinear patch's mixed derivative as exact binary-float fractions."""
    b = inputs['off']['base']['bspline']
    if b['p'] != 1 or b['q'] != 1:
        return None
    (u, v) = map(Fraction, inputs['steps'][0]['uv'])
    (U, V) = [list(map(Fraction, b[k])) for k in ['u', 'v']]
    i = next((i for i in range(1, len(U) - 2) if U[i] <= u < U[i + 1]))
    j = next((j for j in range(1, len(V) - 2) if V[j] <= v < V[j + 1]))
    points = [[[Fraction(x) for x in p[:3]] for p in row] for row in b['net']]
    values = [(points[i][j][k] - points[i][j - 1][k] - (points[i - 1][j][k] - points[i - 1][j - 1][k])) / ((U[i + 1] - U[i]) * (V[j + 1] - V[j])) for k in range(3)]
    return {'base_duv': [str(x) for x in values], 'exactly_zero': all((x == 0 for x in values))}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=['gen0138', 'gen0956'])
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {'method': __doc__, 'analytic_checks': analytic_checks(), 'scope': 'First request only of two nonrational degree-one-v surfaces at interior parameters; no native outputs used; precision convergence is not interval arithmetic or a universal certificate.', 'cases': []}
    for path in sorted(Path(__file__).with_name('offset_grid_cases').glob('*.json')):
        case = json.loads(path.read_text())
        if args.case and case['case'] != args.case:
            continue
        result['cases'].append({'case': case['case'], 'source': case['source'], 'evaluations': [evaluate(case['inputs'], p) for p in [80, 140]], 'exact_mixed_derivative': exact_mixed_derivative(case['inputs'])})
    text = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(text)
    else:
        print(text, end='')
if __name__ == '__main__':
    main()
