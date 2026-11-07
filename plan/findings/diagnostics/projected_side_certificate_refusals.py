import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, PathPatch
from matplotlib.path import Path as MPath

root=Path(sys.argv[1] if len(sys.argv)>1 else '.')
# Coordinates from fixtures/topology.brep_p10/raw_vote_*.smb.
# Midpoints are exact quadratic evaluations, not a geometric approximation.
fig,axes=plt.subplots(1,2,figsize=(11.8,5.4),constrained_layout=True)
for ax,title in zip(axes,['Rectangle: all four vertex votes refuse','Lens: on-curve vertices; both midpoint votes refuse']):
    ax.plot([-2,-1,1,2],[0,0,1,1],color='#b42318',lw=2.3,label='Exact projected cutter')
    ax.scatter([-1.5],[-1],color='#111827',s=35,zorder=5,label='Reference (same side)')
    ax.set(xlim=(-2.1,2.1),ylim=(-1.15,1.15),xlabel='x',ylabel='y',title=title)
    ax.set_aspect('equal'); ax.grid(alpha=.18)
axes[0].add_patch(Rectangle((-.05,.39),.1,.02,facecolor='#c7e9d0',edgecolor='#157347',lw=1.8))
vertices=[(-.05,.39,0),(.05,.39,0),(.05,.41,0),(-.05,.41,0)]
axes[0].scatter([p[0] for p in vertices],[p[1] for p in vertices],color='#157347',s=20,zorder=5)
axes[0].annotate('Face entirely below cutter\nvertex and midpoint queries: 1021',xy=(0,.4),xytext=(-1.9,-.42),arrowprops={'arrowstyle':'->','color':'#157347'},fontsize=9,color='#14532d')
lens=MPath([(-1,0),(0,.1),(1,1),(0,-.1),(-1,0)],[MPath.MOVETO,MPath.CURVE3,MPath.CURVE3,MPath.CURVE3,MPath.CURVE3])
axes[1].add_patch(PathPatch(lens,facecolor='#c7e9d0',edgecolor='#157347',lw=1.8))
midpoints=[(0,.3,0),(0,.2,0)]
axes[1].scatter([p[0] for p in midpoints],[p[1] for p in midpoints],color='#157347',s=32,zorder=6)
axes[1].annotate('Actual midpoint votes\n(0, 0.3), (0, 0.2): 1021',xy=(0,.25),xytext=(-1.95,-.42),arrowprops={'arrowstyle':'->','color':'#157347'},fontsize=9,color='#14532d')
axes[0].legend(loc='upper left',fontsize=8)
fig.suptitle('Before fix: identical geometry, reversed successful trim output',fontsize=14)
fig.text(.5,-.035,'gaps 1.5e-12 / 1.9e-12: KeepPoint: 0 faces, DeletePoint: 1 face\ngap 1e-8 control: KeepPoint: 1 face, DeletePoint: 0 faces',ha='center',fontsize=11)
fig.savefig(root/'projected-side-certificate-refusals.png',dpi=170,bbox_inches='tight')
fig.savefig(root/'projected-side-certificate-refusals.svg',bbox_inches='tight')
