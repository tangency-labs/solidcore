"""Render frozen public gen0472 outputs before/after at identical view and scale.

Curve samples draw geometry only; the finite construction supplies the certificate.
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from check_projected_pair_witnesses_3060 import curve

root = Path(__file__).resolve().parent
row = json.loads((root / "projected_pair_ellipse_after_3060.json").read_text(encoding="utf-8"))
evaluate, domain, _ = curve(row["input"]["c1"])
normal = np.array(row["input"]["vecs"][0], dtype=float)
assert np.array_equal(normal, [0.0, 0.0, 1.0])
points = np.array([evaluate(t) for t in np.linspace(*domain, 401)])
fig, axes = plt.subplots(1, 2, figsize=(11, 5))
for ax, solutions, label in zip(axes,
        [row["before_solutions"], row["solutions"]], ["Before", "Certified finite All"]):
    ax.plot(points[:, 0], points[:, 1], color="#1f6fb2", lw=2.5)
    for solution in solutions:
        a, b = [evaluate(t) for t in solution["parameters"]]
        ax.plot([a[0], b[0]], [a[1], b[1]], "--", color="#555555", lw=1.5)
        ax.scatter([a[0], b[0]], [a[1], b[1]], color="black", s=24, zorder=3)
    value = max(s["value"] for s in solutions)
    ax.set_title(f"{label}: {len(solutions)} ordered pairs\nDistance {value:.12f}", fontsize=13)
    ax.set_xlim(0.2, 7.5)
    ax.set_ylim(-6.5, 0.8)
    ax.set_aspect("equal")
    ax.grid(alpha=0.15)
    ax.set_xlabel("Projected X")
    ax.set_ylabel("Projected Y")
fig.suptitle("gen0472: identical ellipse, Z projection, same view and magnification", fontsize=14)
fig.tight_layout()
fig.savefig(root.parent / "figures" / "projected-ellipse-repair-3060.png", dpi=150)
