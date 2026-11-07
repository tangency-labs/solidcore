"""Independent Decimal affine corner audit of the seven changed #3060 All cases."""
import argparse
import json
from pathlib import Path
from decimal import Decimal, getcontext
getcontext().prec=80
class MP:
 @staticmethod
 def mpf(x):return Decimal.from_float(x) if isinstance(x,float) else Decimal(x)
 @staticmethod
 def sqrt(x):return x.sqrt()
mp=MP()
v=lambda x:[mp.mpf(float(y)) for y in x]
add=lambda a,b:[x+y for x,y in zip(a,b)]
sub=lambda a,b:[x-y for x,y in zip(a,b)]
scale=lambda a,s:[x*s for x in a]
dot=lambda a,b:sum(x*y for x,y in zip(a,b))
def cross(a,b):return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def affine(spec,domain):
 kind,d=next(iter(spec.items()))
 if kind=='line':
  p,dr=v(d['p']),v(d['v']); dr=scale(dr,mp.mpf(d['scale'])/mp.sqrt(dot(dr,dr)))
  if d['inside_out']:p=add(p,scale(dr,sum(v(d['domain']))));dr=scale(dr,-1)
  return p,dr
 assert kind=='bspline' and d['p']==1
 knots=v(d['u']); cps=[v(cp[:3]) for cp in d['cpts']]
 spans=[j for j in range(1,len(cps)) if knots[j]<knots[j+1] and knots[j]<=domain[0] and domain[1]<=knots[j+1]]
 assert len(spans)==1
 j=spans[0]; dr=scale(sub(cps[j],cps[j-1]),1/(knots[j+1]-knots[j])); return sub(cps[j-1],scale(dr,knots[j])),dr
names={'gen0011','gen1860','gen2016','gen2050','gen2187','gen2501','gen2585'}
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('audit_jsonl')
args=parser.parse_args()
results=[]
for line in Path(args.audit_jsonl).read_text(encoding='utf-8-sig').splitlines():
 r=json.loads(line)
 if r['case'] not in names:continue
 i=r['input']; specs=[i['c1'],i['c2'] or i['c1']]; domains=[]
 for spec,requested in zip(specs,[i['ivl1'],i['ivl2']]):
  kind,d=next(iter(spec.items()));natural=d['domain'] if kind=='line' else [d['u'][0],d['u'][-1]]
  domains.append(v([max(natural[0],requested[0]),min(natural[1],requested[1])]))
 (a,da),(b,db)=[affine(spec,d) for spec,d in zip(specs,domains)]
 n=v(i['vecs'][0]); project=lambda a:cross(cross(n,a),n)
 daP,dbP=project(da),project(db)
 assert dot(daP,daP)>0 and dot(dbP,dbP)>0
 corners=[]
 for si,s in enumerate(domains[0]):
  for ti,t in enumerate(domains[1]):
   delta=project(sub(add(a,scale(da,s)),add(b,scale(db,t))))
   ds=2*dot(delta,daP)*(1 if si==0 else -1)
   dt=-2*dot(delta,dbP)*(1 if ti==0 else -1)
   corners.append({'parameters':[float(s),float(t)],'value':float(mp.sqrt(dot(delta,delta))),'inward_derivatives':[float(ds),float(dt)],'strict_maximum':ds<0 and dt<0})
 best=max(x['value'] for x in corners); expected=[c for c in corners if c['strict_maximum'] and c['value']>=best-i['tol']]
 actual=r['solutions'];match=len(expected)==len(actual) and all(any(c['parameters']==s['parameters'] and abs(c['value']-s['value'])<1e-12*max(1,c['value']) for s in actual) for c in expected)
 results.append({'case':r['case'],'corners':corners,'expected_count':len(expected),'actual_count':len(actual),'matches':match,'proof':'squared projected norm is separately strictly convex; strict inward derivative signs classify all corner maxima; unchanged slack selects answer set'})
for result in results: print(json.dumps(result,sort_keys=True))
