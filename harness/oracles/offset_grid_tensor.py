"""Input-only Decimal Taylor diagnostic for interior nonrational offset surfaces.

Original implementation of the tensor B-spline definition and S+d*(Su x Sv)/|Su x Sv|.
Basis definition: https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.BSpline.html
No engine output enters evaluation. Precision convergence is not an interval certificate.
Scope: degrees 1..3, first evaluation request, strictly inside a knot span, orders <=3.
"""
from decimal import Decimal as D, localcontext
from functools import lru_cache
from math import factorial
from pathlib import Path
from fractions import Fraction as F
import json

N = 4
INDICES = sorted(((i, j) for i in range(N + 1) for j in range(N + 1)), key=lambda t: (sum(t), t))
ZERO = (0, 0)

def dec(x):
    return D.from_float(x) if isinstance(x, float) else D(x)

def constant(x):
    return {k: dec(x) if k == ZERO else D(0) for k in INDICES}

def add(a, b):
    return {k: a[k] + b[k] for k in INDICES}

def scale(a, s):
    return {k: a[k] * s for k in INDICES}

def subtract(a, b):
    return add(a, scale(b, D(-1)))

def multiply(a, b):
    return {(i, j): sum((a[k, l] * b[i-k, j-l]
                         for k in range(i+1) for l in range(j+1)), D(0)) for i, j in INDICES}

def derivative(a, axis):
    return {(i, j): ((i+1) * a[i+1, j] if i < N else D(0)) if axis == 0
            else ((j+1) * a[i, j+1] if j < N else D(0)) for i, j in INDICES}

def divide(a, b):
    c = constant(0)
    assert b[ZERO] != 0
    for i, j in INDICES:
        others = sum((b[k, l] * c[i-k, j-l] for k in range(i+1) for l in range(j+1)
                      if (k, l) != ZERO), D(0))
        c[i, j] = (a[i, j] - others) / b[ZERO]
    return c

def square_root(a):
    b = constant(0)
    assert a[ZERO] > 0, 'the definition has no regular normal here'
    b[ZERO] = a[ZERO].sqrt()
    for i, j in INDICES[1:]:
        others = sum((b[k, l] * b[i-k, j-l] for k in range(i+1) for l in range(j+1)
                      if (k, l) not in (ZERO, (i, j))), D(0))
        b[i, j] = (a[i, j] - others) / (2 * b[ZERO])
    return b

def cross(a, b):
    return [subtract(multiply(a[1], b[2]), multiply(a[2], b[1])),
            subtract(multiply(a[2], b[0]), multiply(a[0], b[2])),
            subtract(multiply(a[0], b[1]), multiply(a[1], b[0]))]

def basis_coefficients(knots, degree, parameter, count):
    x = [parameter, D(1)] + [D(0)] * (N - 1)
    def mul(a, b):
        return [sum((a[k] * b[i-k] for k in range(i+1)), D(0)) for i in range(N+1)]
    @lru_cache(None)
    def basis(i, p):
        if p == 0:
            return (D(int(knots[i] <= parameter < knots[i+1])),) + (D(0),) * N
        result = [D(0)] * (N+1)
        for offset, start, end, sign in [(0, i, i+p, 1), (1, i+1, i+p+1, -1)]:
            width = knots[end] - knots[start]
            if width == 0:
                continue
            linear = x.copy()
            linear[0] -= knots[start] if sign == 1 else knots[end]
            term = mul(linear, basis(i+offset, p-1))
            result = [a + b * sign / width for a, b in zip(result, term)]
        return tuple(result)
    return [basis(i, degree) for i in range(count)]

def evaluate(inputs, precision):
    with localcontext() as context:
        context.prec = precision
        off = inputs['off']; b = off['base']['bspline']; step = inputs['steps'][0]
        assert 1 <= b['p'] <= 3 and 1 <= b['q'] <= 3
        assert 0 <= step['hu'] <= 3 and 0 <= step['hv'] <= 3
        assert step['upper'] or step['hu'] + step['hv'] <= 3
        assert all(p[3] == -float.fromhex('0x1.fffffffffffffp+1023') for row in b['net'] for p in row)
        U, V = [dec(x) for x in b['u']], [dec(x) for x in b['v']]
        u, v = map(dec, step['uv']); nu, nv = len(b['net']), len(b['net'][0])
        assert U[b['p']] < u < U[nu] and V[b['q']] < v < V[nv]
        assert u not in U and v not in V, 'one-sided knot derivatives excluded'
        Bu = basis_coefficients(U, b['p'], u, nu); Bv = basis_coefficients(V, b['q'], v, nv)
        points = [[[dec(x) for x in p[:3]] for p in row] for row in b['net']]
        surface = [{(i, j): sum((points[k][l][c] * Bu[k][i] * Bv[l][j]
                                for k in range(nu) for l in range(nv)), D(0))
                    for i, j in INDICES} for c in range(3)]
        normal = cross([derivative(s, 0) for s in surface], [derivative(s, 1) for s in surface])
        length = square_root(add(add(multiply(normal[0], normal[0]), multiply(normal[1], normal[1])), multiply(normal[2], normal[2])))
        result = [add(s, scale(divide(n, length), dec(off['d']))) for s, n in zip(surface, normal)]
        wanted = [(i, j) for i in range(step['hu']+1) for j in range(step['hv']+1)
                  if not step['upper'] or i+j <= max(step['hu'], step['hv'])]
        return {'precision': precision, 'normal_length': str(length[ZERO]),
                'derivatives': {f'{i},{j}': [str(r[i, j] * factorial(i) * factorial(j)) for r in result] for i, j in wanted}}

def analytic_checks():
    sentinel = -float.fromhex('0x1.fffffffffffffp+1023')
    inp = {'off': {'d': 0.75, 'base': {'bspline': {
        'p': 2, 'q': 2, 'u': [0,0,0,1,1,1], 'v': [0,0,0,1,1,1],
        'net': [[[i/2,j/2,int(i==2)+int(j==2),sentinel] for j in range(3)] for i in range(3)]}}},
        'steps': [{'hu': 3, 'hv': 3, 'upper': True, 'uv': [0.25,0.5]}]}
    result = evaluate(inp, 140)['derivatives']
    with localcontext() as ctx:
        ctx.prec = 140
        u,v,d = D('0.25'),D('0.5'),D('0.75');g=1+4*u*u+4*v*v;root=g.sqrt()
        expected = {'0,0': [u-2*d*u/root,v-2*d*v/root,u*u+v*v+d/root],
                    '1,1': [d*(8*v/(g*root)-96*u*u*v/(g*g*root)),
                            d*(8*u/(g*root)-96*u*v*v/(g*g*root)),48*d*u*v/(g*g*root)],
                    '3,0': [d*(24/(g*root)-576*u*u/(g*g*root)+1920*u**4/(g*g*g*root)),
                            d*(-288*u*v/(g*g*root)+1920*u**3*v/(g*g*g*root)),
                            d*(144*u/(g*g*root)-960*u**3/(g*g*g*root))]}
        for key, values in expected.items():
            assert max(abs(D(a)-b) for a,b in zip(result[key], values)) < D('1e-125')
    return {'paraboloid': 'point, mixed uv, third u derivative agree with independent closed forms to absolute 1e-125 at 140 digits'}

def bilinear_closed_form(inputs, precision=140):
    """Exact rational Su, Sv and Suv, followed by the closed form of offset Duu.

    A bilinear surface has n=Su x Sv affine in u, so n_uu=0. Writing
    A=n.n_u, B=n_u.n_u, L=|n| gives (n/L)_uu=-(2*n_u*A+n*B)/L^3+3*n*A^2/L^5.
    This calculation has no spline-basis recurrence or Taylor-series operations.
    """
    b=inputs['off']['base']['bspline']
    if b['p'] != 1 or b['q'] != 1:
        return None
    u,v=map(F,inputs['steps'][0]['uv']); U,V=[list(map(F,b[k])) for k in ['u','v']]
    i=next(i for i in range(1,len(U)-2) if U[i]<=u<U[i+1])
    j=next(j for j in range(1,len(V)-2) if V[j]<=v<V[j+1])
    p=[[[F(c) for c in point[:3]] for point in row] for row in b['net']]
    du,dv=U[i+1]-U[i],V[j+1]-V[j]; alpha,beta=(u-U[i])/du,(v-V[j])/dv
    delta=[p[i][j][k]-p[i][j-1][k]-p[i-1][j][k]+p[i-1][j-1][k] for k in range(3)]
    su=[(p[i][j-1][k]-p[i-1][j-1][k]+beta*delta[k])/du for k in range(3)]
    sv=[(p[i-1][j][k]-p[i-1][j-1][k]+alpha*delta[k])/dv for k in range(3)]
    suv=[x/(du*dv) for x in delta]
    def vec_cross(a,b):
        return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
    def dot(a,b):
        return sum((x*y for x,y in zip(a,b)),F(0))
    n=vec_cross(su,sv); nu=vec_cross(su,suv); L2=dot(n,n); A=dot(n,nu); B=dot(nu,nu)
    assert L2 > 0
    with localcontext() as ctx:
        ctx.prec=precision
        cv=lambda x:D(x.numerator)/D(x.denominator)
        L=cv(L2).sqrt(); distance=dec(inputs['off']['d'])
        duu=[distance*(-cv(2*nu[k]*A+n[k]*B)/(L*L*L)+3*cv(n[k]*A*A)/(L*L*L*L*L)) for k in range(3)]
        return {'base_duv':[str(x) for x in suv], 'offset_duu':[str(x) for x in duu],
                'normal_length':str(L), 'exact_affine':all(x==0 for x in suv)}

def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {'method': __doc__, 'analytic_checks': analytic_checks(), 'cases': []}
    for path in sorted(Path(__file__).with_name('offset_grid_tensor_cases').glob('*.json')):
        case = json.loads(path.read_text())
        if args.case and case['case'] != args.case:
            continue
        evaluations = [evaluate(case['inputs'], precision) for precision in [80, 140]]
        closed = bilinear_closed_form(case['inputs'])
        if closed and '2,0' in evaluations[-1]['derivatives']:
            with localcontext() as context:
                context.prec = 160
                # Compare two independently formulated calculations; exact zeros are
                # reported by the rational closed form, not accepted through a floor.
                for series, expected in zip(evaluations[-1]['derivatives']['2,0'], closed['offset_duu']):
                    if D(expected):
                        assert abs(D(series)-D(expected))/abs(D(expected)) < D('1e-100')
        result['cases'].append({'case': case['case'], 'source': case['source'],
                                'evaluations': evaluations, 'bilinear_closed_form': closed})
    if not result['cases']:
        parser.error('no captured input matches the requested case')
    text = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(text)
    else:
        print(text, end='')

if __name__ == '__main__':
    main()
