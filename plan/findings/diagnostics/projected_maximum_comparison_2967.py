"""Render the reported sphere query before and after the accuracy correction.

The reported scalars are from issue 2967 and the regression output respectively.
Curves and the maximizing point are independent analytic geometry, not sampled
solver output. The lower panels magnify distance error around the true maximum;
neither plot is the numerical certificate used by the implementation.
"""

from pathlib import Path
import json

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
    draw_torus()
    draw_vertex_answers()


def draw_vertex_answers():
    """Draw a representative traversal record's omitted tolerance answers."""
    root = Path(__file__).resolve().parents[3]
    fixture = root / "fixtures/intersect.topology_solver_p2/static-subdivided-shell-vertices.json"
    case = next(c for c in json.loads(fixture.read_text())["cases"]
                if c["name"] == "subdivided_0_14_1_12_11_2_1")["inputs"]
    points = np.array(case["first"])[:, :2]
    query = np.array(case["point"])[:2]
    distances = np.linalg.norm(points - query, axis=1)
    limit = min(distances) + case["tolerance"]
    # The saved parent output returns object IDs 7, 8 and 3. The new output's
    # eleven object IDs equal this independently enumerated tolerance set.
    before = [7, 8, 3]
    after = np.flatnonzero(distances <= limit)
    assert len(after) == 11
    theta = np.linspace(0, 2 * np.pi, 1001)
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), layout="constrained")
    panels = [("Before: 3 returned", before), ("After: all 11 within tolerance", after)]
    # Enforce strict pairing while retaining Python 3.9 support.
    if len(axes) != len(panels):
        raise ValueError("Each comparison panel requires exactly one axis")
    for ax, (title, selected) in zip(axes, panels):
        ax.plot(query[0] + limit * np.cos(theta), query[1] + limit * np.sin(theta),
                "--", color="#999999", label="Best distance + caller tolerance")
        ax.scatter(points[:, 0], points[:, 1], marker="x", color="#999999", label="Fixture vertices")
        ax.scatter(points[selected, 0], points[selected, 1], color="#1f6fb2", s=35, label="Returned answers")
        ax.scatter(*query, color="black", marker="*", s=80, label="Query point")
        ax.set(xlim=(-3, 3.4), ylim=(-3, 3.4), xlabel="x", ylabel="y", title=title)
        ax.set_aspect("equal")
        ax.legend(loc="upper right", fontsize=8)
    fig.suptitle("Projected point/vertices minimum: unchanged tolerance 2, unchanged best 0.7615773106", fontsize=12)
    output = Path(__file__).resolve().parents[1] / "figures" / "projected-tolerance-answers-2967.png"
    fig.savefig(output, dpi=160)
    plt.close(fig)


def draw_torus():
    """Show the other changed analytic maximum against its capsule closed form."""
    query = np.array([0.2, 0.3])
    center = np.array([0.0, -1.0])
    radius = 0.25
    exact = np.linalg.norm(query - center) + radius
    peak = np.arctan2(*(center - query)[::-1])
    farthest = center + radius * np.array([np.cos(peak), np.sin(peak)])
    theta = np.linspace(0, 2 * np.pi, 1001)
    outline_x = radius * np.cos(theta)
    outline_y = np.where(np.sin(theta) >= 0, 1, -1) + radius * np.sin(theta)
    delta = np.linspace(-0.8, 0.8, 1001)
    angle = peak + np.deg2rad(delta)
    distance = np.sqrt((radius * np.cos(angle) - query[0]) ** 2
                       + (-1 + radius * np.sin(angle) - query[1]) ** 2)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
    for col, (title, reported) in enumerate([
        ("Before: fitted spatial silhouette", 1.5652926944322216),
        ("After: surface candidates + domain bound", 1.5652946437965907),
    ]):
        ax = axes[0, col]
        ax.plot(np.r_[outline_x, outline_x[0]], np.r_[outline_y, outline_y[0]], color="#1f6fb2", lw=2)
        ax.plot([query[0], farthest[0]], [query[1], farthest[1]], "--", color="#999999")
        ax.scatter(*query, color="black", s=25)
        ax.scatter(*farthest, color="#555555", s=25)
        ax.set(xlim=(-0.9, 0.9), ylim=(-1.4, 1.4), xlabel="x", ylabel="y", title=title)
        ax.set_aspect("equal")
        ax = axes[1, col]
        ax.axhspan(-10, 10, color="#dce9f5", label="Caller tolerance: 1e-5")
        ax.plot(delta, (distance - exact) * 1e6, color="#1f6fb2", lw=2, label="Exact capsule distance")
        ax.axhline((reported - exact) * 1e6, color="black", linestyle="--", label="Reported maximum")
        ax.set(xlim=(-0.8, 0.8), ylim=(-22, 12), xlabel="Angle from true maximum (degrees)",
               ylabel="Distance minus exact maximum (1e-6)")
        ax.set_title(f"Reported: {reported:.16f}\nAbsolute error: {abs(reported-exact):.3g}", fontsize=11)
        ax.legend(loc="lower center", fontsize=9)
    fig.suptitle("Torus projected along Z: exact maximum = sqrt(0.2² + 1.3²) + 0.25", fontsize=14)
    output = Path(__file__).resolve().parents[1] / "figures" / "projected-torus-comparison-2967.png"
    fig.savefig(output, dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
