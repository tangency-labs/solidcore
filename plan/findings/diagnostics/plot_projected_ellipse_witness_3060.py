"""Draw the original ellipse All result and a larger feasible pair; no solve proof."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from check_projected_pair_witnesses_3060 import curve
spec={'ellipse':{'c':[1.9892188469923955,-3.246996295465655,-3.4858186076942532],
 'x':[1.,0.,0.],'y':[0.,1.,0.],'radii':[4.723327346620286,2.505285780586296],
 'domain':[-90.,101.59684508222117],'inside_out':False}}
evaluate,_,_=curve(spec)
points=np.array([evaluate(t) for t in np.linspace(0.,1.,401)])
fig,axes=plt.subplots(1,2,figsize=(11,5))
for ax,pair,label in zip(axes,[(0.,0.9378100257810158),(.3812096251375738,1.)],['Reported All pair','Larger admissible pair']):
 a,b=[evaluate(t) for t in pair];distance=np.linalg.norm((a-b)[:2])
 ax.plot(points[:,0],points[:,1],color='#1f6fb2',lw=2.5)
 ax.plot([a[0],b[0]],[a[1],b[1]],'--',color='#555555',lw=1.5)
 ax.scatter([a[0],b[0]],[a[1],b[1]],color='black',s=24,zorder=3)
 ax.set_title(f'{label}\nDistance {distance:.12f}',fontsize=13)
 ax.set_xlim(.2,7.5);ax.set_ylim(-6.5,.8);ax.set_aspect('equal');ax.grid(alpha=.15)
 ax.set_xlabel('X');ax.set_ylabel('Y')
fig.suptitle('gen0472: identical ellipse, Z projection, same view and magnification',fontsize=14)
fig.tight_layout()
fig.savefig(Path(__file__).resolve().parent.parent/'figures'/'projected-ellipse-witness-3060.png',dpi=150)
