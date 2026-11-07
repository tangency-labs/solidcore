# Exact rational check of actual analytic-Line defining data and directed bounds.
# The square-root decimals below are display only; acceptance uses Fraction inequalities.
import json
from pathlib import Path
import sys
from fractions import Fraction as F
from decimal import Decimal,localcontext

def cross(a,b):return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def dec(q):return Decimal(q.numerator)/Decimal(q.denominator)
with localcontext() as ctx:
 ctx.prec=90
 for line in (Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('projected_trim_affine_coefficients.jsonl')).read_text(encoding='utf-8-sig').splitlines():
  r=json.loads(line);n=list(map(F,r['n']));a=list(map(F,r['anchor']));q=list(map(F,r['q']));v=list(map(F,r['direction']));scale=F(r['scale']);nat=list(map(F,r['natural']));step=list(map(F,r['step']));assert not r['inside_out']
  slope=(step[1]-step[0])/(nat[1]-nat[0]);intercept=step[0]-slope*nat[0]
  a=[x+scale*intercept*y for x,y in zip(a,v)];d=[scale*slope*y for y in v]
  x=cross(cross(n,[x-y for x,y in zip(a,q)]),n);d=cross(cross(n,d),n)
  dd=dot(d,d);xd=dot(x,d);t=-xd/dd;opt2=dot(x,x)-xd*xd/dd
  assert nat[0]<t<nat[1]
  candidates=[]
  for c in r['candidates']:
   ct=F(c['t']);assert nat[0]<=ct<=nat[1]
   obj2=dot([xx+ct*dx for xx,dx in zip(x,d)],[xx+ct*dx for xx,dx in zip(x,d)])
   lo,hi=map(F,c['own']);assert lo*lo<=opt2<=hi*hi;assert lo*lo<=obj2<=hi*hi
   centered=F(c['residual_lower']);assert 0<=centered and centered*centered<=opt2
   assert centered>=lo and c['residual_certifies']
   candidates.append({'t':c['t'],'objective_gap':str(dec(obj2).sqrt()-dec(opt2).sqrt()),'own_contains_exact_optimum':True})
  print(json.dumps({'case':r['case'],'point':r['q'],'exact_t_decimal':str(dec(t)),'exact_minimum_decimal':str(dec(opt2).sqrt()),'carrier_lower':r['carrier'][0],'candidate_own_lower':r['candidates'][0]['own'][0],'candidates':candidates}))
