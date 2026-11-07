"""Render identical-view geometry before/after for the frozen #3060 audit.

Sampling draws geometry only; the construction/audit records provide the proof.
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from check_projected_pair_witnesses_3060 import curve

ROOT=Path(__file__).resolve().parent
OUT=ROOT.parent/'figures'
def rows(name):return [json.loads(s) for s in (ROOT/name).read_text(encoding='utf-8').splitlines()]
before={r['case']:r for r in rows('projected_pair_single_before_3060.jsonl')}
cases=[]
for r in rows('projected_pair_single_after_3060.jsonl'):
 if r['unit']=='curve.curve_p9' and r['output']!=before[r['case']]['output']:
  cases.append((r,before[r['case']]['solutions']))
for r in rows('projected_pair_affine_all_after_3060.jsonl'):cases.append((r,r['before_solutions']))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9})
for page,start in enumerate(range(0,len(cases),6),1):
 group=cases[start:start+6];fig,axes=plt.subplots(len(group),2,figsize=(10,3.1*len(group)),squeeze=False)
 for (r,old),pair in zip(group,axes):
  i=r['input'];a,da,_=curve(i['c1']);b,db,_=curve(i['c2'] or i['c1']);n=np.array(i['vecs'][0],dtype=float)
  nhat=n/np.linalg.norm(n) if np.linalg.norm(n)>0 else np.array([0.,0.,1.])
  helper=np.eye(3)[np.argmin(abs(nhat))];x=np.cross(nhat,helper);x/=np.linalg.norm(x);y=np.cross(nhat,x)
  def project(p):
   p=np.cross(np.cross(n,p),n);return np.array([np.dot(p,x),np.dot(p,y)])
  domains=[(max(d[0],v[0]),min(d[1],v[1])) for d,v in zip([da,db],[i['ivl1'],i['ivl2']])]
  samples=[np.array([project(c(t)) for t in np.linspace(lo,hi,241)]) for c,(lo,hi) in zip([a,b],domains)]
  allpoints=np.vstack(samples);mins=allpoints.min(axis=0);maxs=allpoints.max(axis=0);width=max(float(max(maxs-mins)),1e-8);center=(mins+maxs)/2
  for ax,sols,label in zip(pair,[old,r['solutions']],['Before','After']):
   for pts,color in zip(samples,['#1f6fb2','#555555']):ax.plot(*pts.T,color=color,lw=2)
   for sol in sols or []:
    s,t=sol['parameters'];p,q=project(a(s)),project(b(t));ax.plot([p[0],q[0]],[p[1],q[1]],color='#111111',lw=1,ls='--');ax.scatter([p[0],q[0]],[p[1],q[1]],s=18,color='#111111',zorder=3)
   vals=(f"max {max(s['value'] for s in sols):.9g} (n={len(sols)})" if sols else 'empty')
   ax.set_title(f"{r['case']} {label}: {vals}",loc='left',fontsize=10)
   ax.set_xlim(center[0]-.6*width,center[0]+.6*width);ax.set_ylim(center[1]-.6*width,center[1]+.6*width);ax.set_aspect('equal');ax.grid(alpha=.15);ax.set_xlabel('Projected coordinate 1');ax.set_ylabel('Projected coordinate 2')
 fig.suptitle('Projected curve-pair candidates: same geometry, view and magnification',fontsize=13,y=1)
 fig.tight_layout();fig.savefig(OUT/f'projected-pair-audit-3060-{page}.png',dpi=125,bbox_inches='tight');plt.close(fig)
 print(page,','.join(r['case'] for r,_ in group))
