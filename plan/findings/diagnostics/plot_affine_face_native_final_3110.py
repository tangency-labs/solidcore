"""Compare the historical unattained pair and the certified native repair."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    """Render actual native endpoints at one fixed projected view and scale."""
    source = Path(__file__).resolve().parent
    before = json.loads((source / "torus_pair_native_audit_3109.json").read_text())
    final = json.loads((source / "affine_face_native_final_3110.json").read_text())
    record = next(row for row in final["records"] if row["case"] == "analytic_face_torus_x_6_15")
    native = record["native"]
    angles = np.linspace(0, np.pi, 201)
    x = np.r_[.25 * np.cos(angles), .25 * np.cos(angles + np.pi)]
    y = np.r_[1 + .25 * np.sin(angles), -1 + .25 * np.sin(angles + np.pi)]
    points = [before["full_torus_native_point"][:2], native["surface_point"][:2]]
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    titles = ["Before: fitted scalar is not attained", "After: native pair and whole-domain proof"]
    for axis, point, title in zip(axes, points, titles):
        axis.fill(x, y, color="#f0f0f0")
        axis.plot(np.r_[x, x[0]], np.r_[y, y[0]], color="#1f6fb2", lw=2)
        axis.plot([-1, 1], [0, 0], color="#1f6fb2", lw=2)
        axis.plot([-1, point[0]], [0, point[1]], "--", color="#555555")
        axis.scatter([-1, point[0]], [0, point[1]], color="black", zorder=5)
        axis.set(xlim=(-1.1, 1.1), ylim=(-1.4, 1.7), xlabel="Projected X",
                 ylabel="Projected Y", title=title, aspect="equal")
        axis.grid(alpha=.15)
    axes[0].text(-1.04, 1.43, "Reported: 1.664215486667304\nNative:    1.664213545940380", fontsize=9)
    axes[1].text(-1.04, 1.43, f"Native/report: {native['reported']:.16f}\nBoth endpoint domains certified", fontsize=8)
    fig.suptitle("Actual outputs, same projected view and scale; no retention slack in the proof")
    fig.tight_layout()
    fig.savefig(source.parent / "figures/affine-face-native-final-3110.png", dpi=160)


if __name__ == "__main__":
    main()
