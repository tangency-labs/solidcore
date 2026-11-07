"""Actual projected geometry and scalar discrepancy for native torus pair maxima."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path(__file__).resolve().parent
r=json.loads((p/'torus_pair_native_audit_3109.json').read_text())
fig,axes=plt.subplots(1,2,figsize=(11,5))
a=np.linspace(0,2*np.pi,401)
# Projection of this torus is a radius-.25 capsule around y in [-1,1].
x=.25*np.cos(a); y=np.sin(a)+.25*np.sin(a)
# Support boundary: the two semicircles joined by vertical segments.
x=np.r_[.25*np.cos(np.linspace(0,np.pi,201)),.25*np.cos(np.linspace(np.pi,2*np.pi,201))]
y=np.r_[1+.25*np.sin(np.linspace(0,np.pi,201)),-1+.25*np.sin(np.linspace(np.pi,2*np.pi,201))]
q=np.array([-1.,0.]);actual=np.array(r['full_torus_native_point'][:2]);exact=np.array([.25/np.sqrt(2),-1-.25/np.sqrt(2)])
for ax,pt,title in zip(axes,[actual,exact],['Current returned native pair','Exact global reference (not a new solver output)']):
 ax.fill(x,y,color='#f0f0f0');ax.plot(np.r_[x,x[0]],np.r_[y,y[0]],color='#1f6fb2',lw=2)
 ax.plot([-1,1],[0,0],color='#1f6fb2',lw=2)
 ax.plot([q[0],pt[0]],[q[1],pt[1]],'--',color='#555555');ax.scatter([q[0],pt[0]],[q[1],pt[1]],color='black',zorder=5)
 ax.set(xlim=(-1.1,1.1),ylim=(-1.4,1.7),xlabel='Projected X',ylabel='Projected Y',title=title);ax.set_aspect('equal');ax.grid(alpha=.15)
axes[0].text(-1.04,1.43,'Reported: 1.664215486667304\nNative:    1.664213545940380',fontsize=9)
axes[1].text(-1.04,1.43,'Bound attained: sqrt(2) + .25\n1.664213562373095',fontsize=9)
fig.suptitle('Same projected view and scale: fitted scalar is not attained on the surface')
fig.tight_layout();fig.savefig(p.parent/'figures/torus-pair-attainment-3109.png',dpi=160)
