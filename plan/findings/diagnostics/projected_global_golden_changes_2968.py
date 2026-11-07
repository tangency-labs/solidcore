"""Reproduce the 55 current point/curve golden explanations and before/after plots.

NumPy polynomial roots independently enumerate stationary candidates on each
rational span. Endpoints complete the compact-domain extrema oracle; plots are
diagnostics, not substitutes for the stationary-polynomial calculation.
"""
import json
from fractions import Fraction
from decimal import Decimal, localcontext
from pathlib import Path
import numpy as np
from numpy.polynomial import Polynomial as P
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def spans(inp):
    """Expand every clipped span into homogeneous coordinate polynomials."""
    kind,c=next(iter(inp['curve'].items()))
    lo,hi=inp['ivl']
    if kind=='hermite':
        lo,hi=max(lo,0),min(hi,1)
        a,b,d,e=map(lambda key:np.array(c[key],float),['p1','p2','d1','d2'])
        coeff=np.array([a,d,-3*a+3*b-2*d-e,2*a-2*b+d+e])
        polys=[P(coeff[:,i])(P([lo,hi-lo])) for i in range(3)]
        if c['dim']==2: polys[2]=P([0])
        return [(lo,hi,polys,P([1]))]
    if kind=='line':
        lo,hi=max(lo,c['domain'][0]),min(hi,c['domain'][1])
        v=np.array(c['v'],float);v=v/np.linalg.norm(v)*c['scale']
        a=np.array(c['p'],float)
        if c['inside_out']:
            a=a+sum(c['domain'])*v;v=-v
        return [(lo,hi,[P([a[i]+lo*v[i],(hi-lo)*v[i]]) for i in range(3)],P([1]))]
    cp=np.array(c['cpts'],float);u=c['u'];degree=c['p']
    # NOZ is the absent-coordinate sentinel. Dim metadata alone does not erase
    # an explicitly supplied Z coordinate (gen0975 supplies ordinary Z values).
    cp[cp[:,2]>1e300,2]=0
    if np.all(cp[:,3]<-1e300):cp[:,3]=1
    result=[]
    for a,b in zip(u,u[1:]):
        left,right=max(a,lo),min(b,hi)
        if right<=left:continue
        t=P([left,right-left]);mid=(left+right)/2
        cache={}
        def basis(i,p):
            if (i,p) in cache:return cache[i,p]
            if p==0:r=P([float(u[i]<=mid<u[i+1])])
            else:
                r=P([0])
                d=u[i+p]-u[i]
                if d:r=r+(t-u[i])/d*basis(i,p-1)
                d=u[i+p+1]-u[i+1]
                if d:r=r+(u[i+p+1]-t)/d*basis(i+1,p-1)
            cache[i,p]=r;return r
        coords=[sum((basis(i,degree)*cp[i,j] for i in range(len(cp))),P([0])) for j in range(4)]
        result.append((left,right,coords[:3],coords[3]))
    return result

def stored_domain(case):
    """Read the exact stored interval that the solver cache covers."""
    kind,c=next(iter(case["curve"].items()))
    if kind=="hermite": return 0.,1.
    if kind=="line": return tuple(c["domain"])
    return c["u"][c["p"]],c["u"][-c["p"]-1]


def unpack(raw):
    """Decode the runner point-count and six-value solution records."""
    assert len(raw)==1+int(raw[0])*6
    return [(raw[6+6*i],raw[5+6*i]) for i in range(int(raw[0]))]


def oracle(case):
    """Enumerate distance extrema from span endpoints and polynomial roots."""
    normal=np.array(case["vecs"][0]); point=np.array(case["p"])
    matrix=np.dot(normal,normal)*np.eye(3)-np.outer(normal,normal)
    pieces=spans(case); candidates=[]; plotted=[]
    def metric(piece,s):
        a,b,coords,w=piece
        return float(np.linalg.norm(matrix@(point-np.array([c(s) for c in coords])/w(s))))
    for piece in pieces:
        a,b,coords,w=piece
        q=[sum((matrix[i,j]*(point[j]*w-coords[j]) for j in range(3)),P([0])) for i in range(3)]
        numerator=sum((p*p for p in q),P([0]))
        derivative=numerator.deriv()*w-2*numerator*w.deriv()
        roots=derivative.roots()
        params=[0.,1.]
        for root in roots:
            if abs(root.imag)<1e-7 and 0<root.real<1:
                params.append(float(root.real))
        for x in params:candidates.append((metric(piece,x),a+x*(b-a)))
        xs=np.linspace(0,1,301)
        plotted.append((a+xs*(b-a),np.array([metric(piece,x) for x in xs])))
    best=(min if case["sop"]==14 else max)(candidates)
    def evaluate(t):
        piece=next(p for p in pieces if p[0]<=t<=p[1])
        return metric(piece,(t-piece[0])/(piece[1]-piece[0]))
    return best,evaluate,plotted,candidates


def exact_line_refresh(record):
    """Check refreshed points against the exact stored affine definition.

    This independent rational diagnostic does not call the private production
    interval certificate. Decimal square roots only display exact rational
    squared-distance results; the golden-comparison check uses exact squares.
    """
    r=record["line_refresh"]
    def cross(a,b):
        return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
    def dot(a,b):
        return sum(x*y for x,y in zip(a,b))
    def dec(x):
        return Decimal(x.numerator)/Decimal(x.denominator)
    f=Fraction
    n=list(map(f,r["n"]));q=list(map(f,r["q"]))
    anchor=list(map(f,r["anchor"]));v=list(map(f,r["direction"]))
    nat=list(map(f,r["natural"]));step=list(map(f,r["step"]))
    assert not r["inside_out"]
    slope=(step[1]-step[0])/(nat[1]-nat[0])
    intercept=step[0]-slope*nat[0];scale=f(r["scale"])
    a=[aa+scale*intercept*vv-qq for aa,vv,qq in zip(anchor,v,q)]
    d=[scale*slope*vv for vv in v]
    x=cross(cross(n,a),n);y=cross(cross(n,d),n)
    dd=dot(y,y);xd=dot(x,y);t=-xd/dd
    minimum_squared=dot(x,x)-xd*xd/dd
    ct=f(record["after"][6]);reported=f(record["after"][5])
    requested=list(map(f,record["inputs"]["ivl"]))
    assert nat[0]<t<nat[1] and requested[0]<t<requested[1]
    assert nat[0]<ct<nat[1] and requested[0]<ct<requested[1]
    at=[xx+ct*yy for xx,yy in zip(x,y)]
    actual_squared=dot(at,at)
    relative=f(1,10**9)
    # This stricter symmetric band suffices for the unchanged relative rule.
    assert (reported*(1-relative))**2 <= actual_squared <= (reported*(1+relative))**2
    with localcontext() as context:
        context.prec=90
        return {"exact_minimum_parameter":str(dec(t)),
                "exact_minimum_distance":str(dec(minimum_squared).sqrt()),
                "reported_minus_exact_metric":str(dec(reported)-dec(actual_squared).sqrt()),
                "candidate_objective_gap":str(dec(actual_squared).sqrt()-dec(minimum_squared).sqrt()),
                "parameter_admissible":True,
                "unchanged_relative_golden_comparison_passes":True}


def main():
    """Validate all saved records and reproduce seven comparison figures."""
    here=Path(__file__).resolve().parent
    data=json.loads((here/"projected_global_golden_changes_2968.json").read_text())
    report=[]; figures=[]
    for index,(name,record) in enumerate(data["cases"].items()):
        if index%8==0:
            fig,axes=plt.subplots(4,2,figsize=(12,13),layout="constrained")
            figures.append(fig)
        ax=axes.flat[index%8]
        case=record["inputs"]; domain=stored_domain(case)
        op="min" if case["sop"]==14 else "max"
        status=record["new"]["ints"]["status"][0]
        item={"case":name,"operation":op,"requested":case["ivl"],"stored_domain":domain,"tolerance":case["tol"]}
        outside=case["ivl"][0]<domain[0] or case["ivl"][1]>domain[1]
        if status:
            assert status==1021 and outside
            item["reason"]="Requested interval is not covered by the stored-domain cache; status 1021 refuses an uncertified partial or empty search."
            ax.axvspan(*domain,color="#247ab5",alpha=.25,label="Stored domain")
            ax.axvspan(*case["ivl"],color="#db8417",alpha=.3,label="Requested interval")
            for t,value in unpack(record["before"]):ax.scatter(t,0,color="#d34a4a",marker="x",label="Old answer parameter")
            ax.set_yticks([]);ax.set_title(f"{name}: {op}, status 0 -> 1021")
        else:
            assert not outside
            best,evaluate,plotted,candidates=oracle(case)
            item["oracle_distance"],item["oracle_parameter"]=best
            item["returned"]=[]
            for t,value in unpack(record["after"]):
                measured=evaluate(t)
                item["returned"].append({"parameter":t,"reported":value,"evaluated":measured,"gap_to_best":abs(measured-best[0])})
                assert abs(measured-value)<=1e-10*max(1.,abs(value)),(name,t,value,measured)
                assert abs(measured-best[0])<=case["tol"]+1e-10*max(1.,best[0]),(name,measured,best)
            assert item["returned"]
            item["reason"]="All returned parameters evaluate to the reported projected distance and lie within requested retention slack of the independent span-root/end-point extremum."
            for ts,ds in plotted:ax.plot(ts,ds,color="#247ab5",lw=1.6)
            ax.axhline(best[0],color="#555",ls=":",lw=1,label="Independent best")
            for state,color,marker in [("before","#d34a4a","x"),("after","#21985d","o")]:
                vals=unpack(record[state])
                if vals:ax.scatter(*zip(*vals),color=color,marker=marker,s=28,label=state)
                else:ax.plot([],[],color=color,marker=marker,ls="",label=state+": no solution")
            ax.set_title(f"{name}: {op}; best {best[0]:.9g}; slack {case['tol']:g}")
            ax.set_ylabel("Projected distance")
        ax.set_xlabel("Curve parameter");ax.legend(fontsize=7)
        if "line_refresh" in record:
            item["exact_line_refresh"]=exact_line_refresh(record)
        report.append(item)
    for unused in range((len(report)-1)%8+1,8):axes.flat[unused].set_visible(False)
    for i,fig in enumerate(figures,1):
        fig.suptitle("Point/curve golden changes: old red, current green; blue domain/objective",fontsize=12)
        fig.savefig(here.parent/"figures"/f"projected-global-goldens-2968-{i}.png",dpi=140,bbox_inches="tight")
        plt.close(fig)
    (here/"projected_global_golden_changes_2968_report.json").write_text(json.dumps(report,indent=2)+"\n")
    print(f"Validated {len(report)} cases; {sum('oracle_distance' in x for x in report)} numerical extrema and {sum('oracle_distance' not in x for x in report)} explicit domain refusals.")

if __name__=="__main__": main()
