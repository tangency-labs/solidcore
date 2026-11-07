#!/usr/bin/env python3
"""Independent SciPy witness check for #3060 JSONL audits; search is NOT a proof."""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import differential_evolution, minimize


def curve(spec):
    kind, data = next(iter(spec.items()))
    if kind == 'circle':
        kind, data = 'ellipse', {'c': data['center'], 'x': data['x'], 'y': data['y'],
            'radii': [data['radius'], data['radius']],
            'domain': [data['a0'], data['a1']], 'inside_out': False}
    if kind == 'line':
        origin, direction = np.array(data['p']), np.array(data['v'])
        direction = direction / np.linalg.norm(direction) * data['scale']
        lo, hi = data['domain']
        def evaluate(t):
            step = lo + hi - t if data['inside_out'] else t
            return origin + step * direction
        return evaluate, (lo, hi), [lo, hi]
    if kind == 'hermite':
        a, b, end, tangent = (np.array(data[k]) for k in ('p1','d1','p2','d2'))
        c, d = -3*a-2*b+3*end-tangent, 2*a+b-2*end+tangent
        return lambda t: a+t*(b+t*(c+t*d)), (0.,1.), [0.,1.]
    if kind == 'ellipse':
        # Independent homogeneous quadratic-arc construction, at most 90deg/arc.
        x, y = np.array(data['x']), np.array(data['y'])
        x = x / np.linalg.norm(x)
        z = np.cross(x,y); z = z / np.linalg.norm(z)
        y = np.cross(z,x); y = y / np.linalg.norm(y)
        lo, hi = data['domain']
        if data['inside_out']:
            lo, hi, y = 360-hi, 360-lo, -y
        count = int(math.ceil((hi-lo)/90.))
        angles = np.linspace(math.radians(lo), math.radians(hi), count+1)
        c, rx, ry = np.array(data['c']), *data['radii']
        cp = []
        for i in range(count):
            a, b = angles[i:i+2]
            if i == 0: cp.append([*(c+rx*math.cos(a)*x+ry*math.sin(a)*y),1.])
            w = math.cos((b-a)/2)
            mid = (a+b)/2
            cp.append([*(w*c+rx*math.cos(mid)*x+ry*math.sin(mid)*y),w])
            cp.append([*(c+rx*math.cos(b)*x+ry*math.sin(b)*y),1.])
        if hi-lo == 360.: cp[-1] = cp[0].copy()
        knots = [0.]*3 + [v for i in range(1,count) for v in [i/count]*2] + [1.]*3
        data = {'cpts':cp,'u':knots,'p':2}
    elif kind != 'bspline':
        raise ValueError(f'unsupported kind {kind}')
    cp = np.array(data['cpts'], dtype=float)
    cp[:,2] = np.where(cp[:,2] == np.finfo(float).max, 0., cp[:,2])
    cp[:,3] = np.where(cp[:,3] == -np.finfo(float).max, 1., cp[:,3])
    knots = np.array(data['u']); lo, hi = knots[0], knots[-1]
    spline = BSpline(knots, cp, data['p'], extrapolate=False)
    def evaluate(t):
        h = cp[0] if t == lo else cp[-1] if t == hi else spline(t)
        return h[:3] / h[3]
    return evaluate, (lo,hi), sorted(set(knots.tolist()))


def audit(row):
    if row['unit'] != 'curve.curve_p9':
        solution = (row.get('solutions') or [None])[0]
        result = {'case':row['case'], 'analytic_maximum':math.sqrt(2.5),
            'reported':solution['value'] if solution else None,
            'difference':solution['value']-math.sqrt(2.5) if solution else None,
            'proof':'projected rectangle and line endpoint convexity'}
        if solution and len(solution['parameters']) == 3:
            t,u,v = solution['parameters']
            if u == 0.0:
                # Angular parameter zero is the placement X direction in either
                # angular or rational coordinates; height is [0,1].
                a = np.array([-1+2*t,0.,0.]); b = np.array([-.5+v,.5,0.])
                result['measured'] = float(np.linalg.norm(a-b))
                result['admissible'] = 0.<=t<=1. and 0.<=v<=1.
        return result
    i = row['input']; a, da, ka = curve(i['c1']); b, db, kb = curve(i['c2'] or i['c1'])
    intervals = [(max(da[0],i['ivl1'][0]),min(da[1],i['ivl1'][1])),
        (max(db[0],i['ivl2'][0]),min(db[1],i['ivl2'][1]))]
    if any(lo > hi for lo,hi in intervals):
        return {'case':row['case'],'empty_natural_intersection':True}
    normal = np.array(i['vecs'][0])
    def metric(st):
        return float(np.linalg.norm(np.cross(np.cross(normal,a(st[0])-b(st[1])),normal)))
    samples = [sorted(set([*np.linspace(lo,hi,25),*(x for x in knots if lo<=x<=hi)]))
        for (lo,hi),knots in zip(intervals,[ka,kb])]
    probes = [(metric([s,t]),[s,t]) for s in samples[0] for t in samples[1]]
    best, where = max(probes,key=lambda v:v[0])
    result = differential_evolution(lambda st:-metric(st),intervals,seed=3060,
        maxiter=150,popsize=12,tol=1e-11,polish=True,workers=1)
    if -result.fun > best: best,where = -float(result.fun),result.x.tolist()
    # Also polish the strongest grid seeds, including parameter-domain edges.
    for _, seed in sorted(probes,key=lambda v:v[0],reverse=True)[:4]:
        result = minimize(lambda st:-metric(st),seed,bounds=intervals,method='Nelder-Mead',
            options={'xatol':1e-12,'fatol':1e-12,'maxiter':1000})
        if -result.fun > best: best,where = -float(result.fun),result.x.tolist()
    checked=[]
    for solution in row.get('solutions') or []:
        params=solution['parameters']; measured=metric(params)
        checked.append({'parameters':params,'admissible':all(lo<=p<=hi for p,(lo,hi) in zip(params,intervals)),
            'reported':solution['value'],'measured':measured,'scalar_error':measured-solution['value'],
            'larger_witness_gap':best-solution['value'],
            'disproved_beyond_slack':best-solution['value']>i['tol']+1e-10})
    return {'case':row['case'],'checked':checked,'search_maximum':best,'search_parameters':where,
        'target':i['target'],'tolerance':i['tol'],'search_is_certificate':False}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('audit_jsonl');args=p.parse_args()
    for line in Path(args.audit_jsonl).read_text(encoding='utf-8-sig').splitlines():
        if not line.strip(): continue
        row=json.loads(line)
        try: result=audit(row)
        except Exception as exc: result={'case':row['case'],'audit_error':str(exc)}
        print(json.dumps(result,sort_keys=True),flush=True)


if __name__=='__main__': main()
