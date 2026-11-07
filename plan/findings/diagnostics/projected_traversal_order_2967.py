"""Plot the two-band traversal regression with identical geometry and axes."""
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

fig, axes = plt.subplots(1, 2, figsize=(10, 5), layout="constrained")
for ax, title in zip(axes, ["Before: reverse face order refuses", "After: both face orders succeed"]):
    for x, radius, color in [(0, 0.5, "#2870b2"), (0.1, 0.6, "#d96c18")]:
        ax.add_patch(Circle((x, 0), radius, fill=False, color=color, lw=2))
        ax.add_patch(Circle((x, 0), radius / 2, fill=False, color=color, ls="--", alpha=0.6))
    ax.scatter([1e-6], [0], c="black", marker="x", zorder=3)
    ax.annotate("query (0.000001, 0)", (1e-6, 0), xytext=(-0.55, -0.76))
    ax.set(xlim=(-0.72, 0.88), ylim=(-0.82, 0.82), xlabel="projected x", ylabel="projected y", title=title)
    ax.set_aspect("equal")
    ax.grid(alpha=0.15)
axes[0].text(-0.65, 0.7, "Early face bound: Unsupported", fontsize=10)
axes[1].plot([1e-6, 0.7], [0, 0], color="#236b35", lw=2)
axes[1].scatter([0.7], [0], c="#236b35", zorder=3)
axes[1].text(-0.65, 0.7, "Final maximum: 0.699999", fontsize=10)
fig.suptitle("Z projection of the two spherical bands; unchanged tolerance 1e-8")
output = Path(__file__).resolve().parents[1] / "figures" / "projected-traversal-order-2967.png"
fig.savefig(output, dpi=150)
