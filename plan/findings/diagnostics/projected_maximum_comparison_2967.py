"""Render the reported sphere query before and after the accuracy correction.

The reported scalars are from issue 2967 and the regression output respectively.
Curves and the maximizing point are independent analytic geometry, not sampled
solver output. The lower panels magnify distance error around the true maximum;
neither plot is the numerical certificate used by the implementation.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    """Draw both results using identical geometry, axes, and magnification."""
    query = np.array([0.2, 0.3])
    radius = 0.5
    exact = radius + np.linalg.norm(query)
    theta_max = np.arctan2(-query[1], -query[0])
    farthest = radius * np.array([np.cos(theta_max), np.sin(theta_max)])
    theta = np.linspace(0.0, 2.0 * np.pi, 1001)
    delta = np.linspace(-0.8, 0.8, 1001)
    local_theta = theta_max + np.deg2rad(delta)
    distances = np.sqrt(
        (radius * np.cos(local_theta) - query[0]) ** 2
        + (radius * np.sin(local_theta) - query[1]) ** 2
    )
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
    for column, (title, reported) in enumerate([
        ("Before: fitted spatial silhouette", 0.8605419808705617),
        ("After: surface candidates + domain bound", 0.8605551275463991),
    ]):
        ax = axes[0, column]
        ax.plot(radius * np.cos(theta), radius * np.sin(theta), color="#1f6fb2", lw=2)
        ax.plot([query[0], farthest[0]], [query[1], farthest[1]], "--", color="#999999")
        ax.scatter(*query, color="black", s=25)
        ax.scatter(*farthest, color="#555555", s=25)
        ax.annotate("Projected query (0.2, 0.3)", query, xytext=(0.0, 0.59), fontsize=10)
        ax.annotate("Analytic farthest point", farthest, xytext=(-0.64, -0.60), fontsize=10)
        ax.set(xlim=(-0.7, 0.7), ylim=(-0.7, 0.7), xlabel="x", ylabel="y", title=title)
        ax.set_aspect("equal")
        ax = axes[1, column]
        ax.axhspan(-10, 10, color="#dce9f5", label="Caller tolerance: 1e-5")
        ax.plot(delta, (distances - exact) * 1e6, color="#1f6fb2", lw=2, label="Exact equator distance")
        error = (reported - exact) * 1e6
        ax.axhline(error, color="black", linestyle="--", label="Reported maximum")
        ax.set(xlim=(-0.8, 0.8), ylim=(-22, 12), xlabel="Angle from true maximum (degrees)",
               ylabel="Distance minus exact maximum (1e-6)")
        ax.set_title(f"Reported: {reported:.16f}\nAbsolute error: {abs(reported-exact):.3g}", fontsize=11)
        ax.legend(loc="lower center", fontsize=9)
    fig.suptitle("Sphere band projected along Z: exact maximum = 0.8605551275463990", fontsize=14)
    output = Path(__file__).resolve().parents[1] / "figures" / "projected-maximum-comparison-2967.png"
    fig.savefig(output, dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
