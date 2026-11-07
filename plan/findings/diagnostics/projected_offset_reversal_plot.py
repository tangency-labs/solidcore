"""Plot actual native gen1381 geometry with before/after metadata overlays."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
base=Path(__file__).parent
r=json.loads((base/'projected_offset_reversal_native_plot.json').read_text(encoding='utf-8-sig'))
points=r['native_samples'];markers=r['rounded_knot_points']
fig,axes=plt.subplots(1,2,figsize=(11,4.5),sharex=True,sharey=True,constrained_layout=True)
for ax,title in zip(axes,['Before: metadata query refuses (1007)','After: six reflected knot records']):
 ax.plot([p[1] for p in points],[p[2] for p in points],color='#1f6fb2',linewidth=2)
 ax.set(title=title,xlabel='native X');ax.set_aspect('equal',adjustable='box');ax.grid(color='#dddddd',linewidth=.5)
axes[0].set_ylabel('native Y')
axes[1].scatter([p[1] for p in markers],[p[2] for p in markers],color='black',s=24,zorder=3)
for j,p in enumerate(markers):axes[1].annotate(str(j),(p[1],p[2]),xytext=(5,5),textcoords='offset points',fontsize=9)
fig.suptitle('gen1381 — same native offset curve and view; metadata correction only',fontsize=13)
fig.savefig(base/'projected_offset_reversal_native_plot.png',dpi=160)
