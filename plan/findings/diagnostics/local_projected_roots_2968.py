"""Audit and draw the seven changed local projected-distance golden records.

Requires NumPy, SciPy and Matplotlib. Run from any directory. The adjacent JSON
contains public-generator inputs and exact old/new golden output arrays. This
independent spline/Hermite evaluator checks returned stationary roots; plots
are illustrations, not certificates of a global maximum or of nonexistence.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import BSpline
from numpy.polynomial import Polynomial


def evaluator(curve):
    """Construct the supplied rational spline or cubic Hermite definition."""
    if "hermite" in curve:
        h = curve["hermite"]
        a, da, b, db = [np.array(h[k]) for k in ("p1", "d1", "p2", "d2")]
        coeff = np.array([a, da, -3*a-2*da+3*b-db, 2*a+da-2*b+db])
        polys = [Polynomial(coeff[:, i]) for i in range(3)]
        return lambda t, d=0: np.array([p.deriv(d)(t) for p in polys])
    bs = curve["bspline"]
    cp = np.array(bs["cpts"])
    if bs["dim"] == 2:
        cp[:, 2] = 0
    if np.all(cp[:, 3] < -1e300):
        cp[:, 3] = 1
    spline = BSpline(bs["u"], cp, bs["p"])
    def evaluate(t, d=0):
        h = spline(t)
        point = h[:3]/h[3]
        if d == 0:
            return point
        dh = spline(t, nu=1)
        first = (dh[:3]-point*dh[3])/h[3]
        if d == 1:
            return first
        ddh = spline(t, nu=2)
        return (ddh[:3]-2*first*dh[3]-point*ddh[3])/h[3]
    return evaluate


def main():
    here = Path(__file__).resolve().parent
    records = json.loads((here / "local_projected_roots_2968.json").read_text())
    report = []
    for name, record in records.items():
        case = record["inputs"]
        assert case["op"] == "lpa" and case["prop"] == 21
        evaluate = evaluator(case["curve"])
        point, normal = np.array(case["vecs"])
        def project(v):
            return np.cross(np.cross(normal, v), normal)
        def facts(t):
            v = project(evaluate(t)-point)
            tangent = project(evaluate(t, 1))
            second = project(evaluate(t, 2))
            residual = np.dot(v, tangent)/(np.linalg.norm(v)*np.linalg.norm(tangent))
            curvature = np.dot(tangent, tangent)+np.dot(v, second)
            return np.linalg.norm(v), residual, curvature
        n = normal / np.linalg.norm(normal)
        seed = np.eye(3)[np.argmin(abs(n))]
        x = np.cross(n, seed)
        x /= np.linalg.norm(x)
        y = np.cross(n, x)
        ts = np.linspace(*case["ivl"], 1601)
        projected = np.array([project(evaluate(t)-point) for t in ts])
        xy = projected @ np.array([x, y]).T
        distances = np.linalg.norm(projected, axis=1)
        fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
        for column, state in enumerate(("before", "after")):
            out = record[state]
            ax = axes[0, column]
            ax.plot(*xy.T, color="#1f6fb2", lw=2)
            ax.scatter(0, 0, color="black", marker="*", s=65, label="Projected query")
            ax.set_aspect("equal", adjustable="datalim")
            ax.set(xlabel="Projected coordinate 1", ylabel="Projected coordinate 2")
            lower = axes[1, column]
            lower.plot(ts, distances, color="#1f6fb2", lw=2)
            lower.set(xlabel="Curve parameter t", ylabel="Projected distance")
            if out[0]:
                t, reported = out[2:]
                dist, residual, curvature = facts(t)
                assert np.isclose(dist, reported, rtol=2e-12, atol=1e-12), (name, dist, reported)
                assert abs(residual) < 1e-7, (name, residual)
                classification = "minimum" if curvature > 0 else "maximum"
                q = project(evaluate(t)-point) @ np.array([x, y]).T
                ax.scatter(*q, color="#555555", s=45, label="Returned stationary point")
                ax.plot([0,q[0]], [0,q[1]], "--", color="#999999")
                lower.scatter(t, reported, color="#555555", s=45)
                title = f"{state.title()}: stationary {classification}\nt={t:.9g}, distance={reported:.10g}"
                report.append(dict(case=name, state=state, t=t, distance=dist,
                                   normalized_residual=residual, half_squared_distance_second=curvature,
                                   classification=classification))
            else:
                title = f"{state.title()}: found=false\nNo point returned; roots may still exist"
            ax.set_title(title, fontsize=11)
            ax.legend(fontsize=8, loc="best")
        for row in range(2):
            limits_x = [axes[row,c].get_xlim() for c in range(2)]
            limits_y = [axes[row,c].get_ylim() for c in range(2)]
            for c in range(2):
                axes[row,c].set_xlim(min(v[0] for v in limits_x),max(v[1] for v in limits_x))
                axes[row,c].set_ylim(min(v[0] for v in limits_y),max(v[1] for v in limits_y))
        fig.suptitle(f"curve.curve_p5 / {name}: same geometry, interval and projection\nLocal Newton stationary solve; neither output certifies a global maximum", fontsize=13)
        fig.savefig(here.parent / "figures" / f"local-projected-root-{name}-2968.png", dpi=140)
        plt.close(fig)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
