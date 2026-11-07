"""Draw the exact dyadic near-collinear join from review 4501494576."""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-1, 1, 2001)
y = abs(x) / 1024
fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), layout="constrained")
axes[0].plot(x, y, color="#1f6fb2", lw=3)
axes[0].scatter([-1, 0, 1], [1/1024, 0, 1/1024], color="black")
axes[0].set(title="Actual kink (vertical axis magnified)", xlabel="x", ylabel="y")
axes[1].plot(x, np.sqrt(x*x+(y+1)**2), color="#1f6fb2", lw=3)
axes[1].set(title="Minimum query: (0, -1)", xlabel="x", ylabel="Projected distance")
# Stable difference of squares avoids subtracting nearly equal distances.
q = x*x + y*y - 2048*y
gap = -q / (1024 + np.sqrt(1024**2 + q))
axes[2].plot(x, gap, color="#1f6fb2", lw=3)
axes[2].set(title="Maximum query: (0, 1024)", xlabel="x", ylabel="1024 - projected distance")
for ax in axes:
    ax.axvline(0, color="#999999", ls="--", lw=1.5)
    ax.grid(alpha=.15)
fig.suptitle("A 0.112-degree kink is G1-tagged but retains a true extremum", fontsize=16)
fig.savefig(Path(__file__).with_suffix(".png"), dpi=160)
