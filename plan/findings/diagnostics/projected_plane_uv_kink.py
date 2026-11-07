"""Render the exact plane-UV regression objectives; no kernel outputs are inferred."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), sharex=True)
a, delta = 1 / 1024, 1 / 524288
for ax, maximize in zip(axes, [False, True]):
    t = np.linspace(0, 2, 1001)
    x = a * (t - 1)
    y = 1 + (-1 if maximize else 1) * delta * np.abs(t - 1)
    ax.plot(t, (x*x + y*y - 1) * 1e6, color="#305b80", linewidth=2)
    ax.scatter([1], [0], s=110, facecolors="none", edgecolors="#b44740",
               linewidths=2, label="Join omitted before")
    ax.scatter([1], [0], s=24, color="#18805e", zorder=4, label="Exact answer after")
    ax.axhline(0, color="#a8a8a8", linewidth=.6)
    ax.set(title="Projected maximum" if maximize else "Projected minimum",
           xlabel="Native parameter t")
    ax.grid(alpha=.2)
axes[0].set_ylabel("(Distance squared - 1) / 1e-6")
axes[1].legend(fontsize=8, loc="lower center")
fig.tight_layout()
fig.savefig(Path(__file__).resolve().parent.parent / "figures/projected-plane-uv-kink.png", dpi=180)
