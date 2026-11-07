"""Draw the four point/curve domain-refusal fixtures in a shared XY view."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    """Compare recorded answers with refusals on the same analytic geometry."""
    fig, axes = plt.subplots(2, 2, figsize=(9, 9), sharex=True, sharey=True)
    angles = np.linspace(0, 2 * np.pi, 400)
    for row, (rx, ry, name, value) in enumerate([
        (2, 2, "Circle", 5.024937810560445),
        (3, 1, "Ellipse", 6.020797289396148),
    ]):
        for col in range(2):
            ax = axes[row, col]
            ax.plot(rx * np.cos(angles), ry * np.sin(angles), color="#2463a0")
            ax.scatter([0], [-7], color="#333333")
            ax.text(0, -7.7, "P = (0, -7, 0.5)", ha="center", fontsize=9)
            if col == 0:
                ax.scatter([0], [-ry], color="#d1662c")
                ax.plot([0, 0], [-7, -ry], color="#d1662c", ls="--")
                label = f"Accepted t = 0.75\ndistance {value:.9f}"
            else:
                label = "Status 1021\nuncovered requested interval"
            ax.text(0, 3.2, label, ha="center", va="center", fontsize=9)
            ax.set_title(name + (" - before" if col == 0 else " - after"))
            ax.set_aspect("equal")
            ax.set_xlim(-4, 4)
            ax.set_ylim(-8.2, 4)
            ax.grid(alpha=.15)
            ax.set_xlabel("x")
            ax.set_ylabel("y")
    fig.suptitle("Same fixture geometry; invalid parameter request", fontsize=12)
    fig.text(.5, .035,
             "XY schematic (P has z = 0.5). Native parameter: [0, 1]. Requested: [0, 360].\n"
             "Before: only the cached curve is searched. After: uncovered request refused.",
             ha="center", fontsize=10)
    fig.tight_layout(rect=(0, .08, 1, .96))
    destination = Path(__file__).resolve().parent.parent / "figures"
    fig.savefig(destination / "projected-point-curve-domain-refusal.png", dpi=140)


if __name__ == "__main__":
    main()
