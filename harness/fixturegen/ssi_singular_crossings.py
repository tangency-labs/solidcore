#!/usr/bin/env python3
"""Writes fixtures/surface.surface_p2/ssi_singular_crossings.json (issue #393, WR-77).

A reproduction corpus of surface/surface intersections whose branches cross, each a
`gssi` case (`SmSurface::GlobalSurfaceIntersect` of two `SmBSplineSurface`s over their
natural domains, approximation tolerance 1e-6 as in the ssi corpus of #17). Every
surface is an exact NURBS of a closed form, so the curves and their crossing nodes are
known exactly:

- nodes where the surface normals are parallel (the surfaces touch there and two
  branches cross; #393's defect: the marcher steps through such a node without
  ending a curve): two equal cylinders with perpendicular and with oblique crossing
  axes (Steinmetz), a torus and its bitangent plane (two Villarceau circles), the
  saddle `z = x y` and the plane `z = 0`, and a sphere and a cylinder of half its
  radius through its centre (Viviani's figure eight);
- nodes at an apex, where one or both surfaces have no normal: a double cone and a
  plane through its apex, and two double cones sharing their apex.

Rational cylinders, cones, spheres and tori are the 9-point rational circle (or the
5-point half circle) swept or revolved (`w = 1/sqrt(2)` at the corners), the way the
ssi corpus gives them; analytic surfaces would take the special-case intersectors and
never reach the marcher. Each revolved surface's circle starts 20 degrees off the
frame axes, so no node lies on a seam or a knot of the circle (where a curve would end
anyway); the saddle's domain is off centre for the same reason.

Each case's `upstream_ref` states the exact curves and nodes. Every case also runs with
the operands swapped (`-swapped`), as the intersector is not symmetric. The difftest
gate compares the pinned C++ with the default Rust kernel (parity); what the result
does at each node is measured against the closed forms in
`crates/solidcore/tests/issue_393_singular_crossing_corpus.rs`.

Run from the repo root: `python3 -I harness/fixturegen/ssi_singular_crossings.py`.
"""

import json
import math
import os

ROOT = os.getcwd()
OUT = os.path.join(ROOT, "fixtures/surface.surface_p2/ssi_singular_crossings.json")
TOL = 1.0e-6  # the ssi corpus tolerance (#17)
NOT_RATIONAL = -1.7976931348623157e308  # NOW: the control point carries no weight
W = math.sqrt(0.5)
TURN = math.radians(20.0)  # where each rational circle starts, off the frame axes

# the 9-point unit circle (x, y, w) from +x, counter-clockwise, knots at quarters
CIRCLE = [
    (1, 0, 1), (1, 1, W), (0, 1, 1), (-1, 1, W), (-1, 0, 1),
    (-1, -1, W), (0, -1, 1), (1, -1, W), (1, 0, 1),
]
CIRCLE_KNOTS = [0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1]
# the 5-point half circle (radial, axial, w) from -axis to +axis through radial 1
HALF = [(0, -1, 1), (1, -1, W), (1, 0, 1), (1, 1, W), (0, 1, 1)]
HALF_KNOTS = [0, 0, 0, 0.5, 0.5, 1, 1, 1]


def add(a, b):
    return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]


def mul(s, a):
    return [s * a[0], s * a[1], s * a[2]]


def hom(p, w):
    """The stored control point: `w x, w y, w z, w`."""
    return [p[0] * w, p[1] * w, p[2] * w, w]


def plain(p):
    return [float(p[0]), float(p[1]), float(p[2]), NOT_RATIONAL]


def srf(net, p, q, u, v):
    return {"net": net, "p": p, "q": q, "u": [float(x) for x in u], "v": [float(x) for x in v]}


def domain(s):
    return [s["u"][0], s["v"][0], s["u"][-1], s["v"][-1]]


def turned(e1, e2):
    """The frame `(e1, e2)` turned by `TURN` in its plane."""
    c, s = math.cos(TURN), math.sin(TURN)
    return add(mul(c, e1), mul(s, e2)), add(mul(-s, e1), mul(c, e2))


def patch(origin, du, dv, u0, u1, v0, v1):
    """The bilinear patch of the plane `origin + u du + v dv`, `[u0, u1] x [v0, v1]`."""
    net = [[plain(add(origin, add(mul(u, du), mul(v, dv)))) for v in (v0, v1)] for u in (u0, u1)]
    return srf(net, 1, 1, [0, 0, 1, 1], [0, 0, 1, 1])


def revolved(center, e1, e2, axis, profile, q, knots):
    """`profile` (`(radial, axial, w)` control points, degree `q`) revolved about `axis`
    through `center`: u the full turn from `e1` (turned by `TURN`) towards `e2`, v the
    profile. A negative radial runs on the far side of the axis (a double cone)."""
    e1, e2 = turned(e1, e2)
    net = []
    for cx, cy, wu in CIRCLE:
        row = []
        for rho, z, wv in profile:
            p = add(center, add(mul(rho * cx, e1), add(mul(rho * cy, e2), mul(z, axis))))
            row.append(hom(p, wu * wv))
        net.append(row)
    return srf(net, 2, q, CIRCLE_KNOTS, knots)


def cylinder(center, e1, e2, axis, r, h0, h1):
    """The cylinder `center + r (cos e1 + sin e2) + h axis`, `h` in `[h0, h1]`."""
    return revolved(center, e1, e2, axis, [(r, h0, 1), (r, h1, 1)], 1, [0, 0, 1, 1])


def double_cone(e1, e2, axis, slope, h):
    """The double cone `rho = slope |z|` about `axis` through the origin, `z` in
    `[-h, h]`: its apex is the profile's middle (v = 1/2), inside the surface."""
    return revolved([0, 0, 0], e1, e2, axis, [(-slope * h, -h, 1), (slope * h, h, 1)], 1, [0, 0, 1, 1])


def sphere(center, e1, e2, axis, r):
    return revolved(center, e1, e2, axis, [(r * a, r * b, w) for a, b, w in HALF], 2, HALF_KNOTS)


def torus(center, e1, e2, axis, big, small):
    prof = [(big + small * a, small * b, w) for a, b, w in CIRCLE]
    return revolved(center, e1, e2, axis, prof, 2, CIRCLE_KNOTS)


def saddle(x0, x1, y0, y1):
    """`z = x y` over `[x0, x1] x [y0, y1]`: bilinear, so its corners are its net."""
    net = [[plain([x, y, x * y]) for y in (y0, y1)] for x in (x0, x1)]
    return srf(net, 1, 1, [0, 0, 1, 1], [0, 0, 1, 1])


X, Y, Z = [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]
O = [0.0, 0.0, 0.0]
C60, S60 = 0.5, math.sqrt(3.0) / 2.0
B60 = [S60, 0.0, C60]  # 60 degrees off z, in the xz plane
B60_PERP = [-C60, 0.0, S60]  # B60 x Y
C30, S30 = S60, 0.5
D30 = [0.0, C30, S30]  # the bitangent plane's second direction (30 degrees off y)
T20 = math.radians(20.0)
VIVIANI_POLES = [math.sin(T20), math.cos(T20), 0.0]  # the sphere's axis, 20 degrees off y
VIVIANI_E1 = [0.0, 0.0, 1.0]
VIVIANI_E2 = [math.cos(T20), -math.sin(T20), 0.0]  # VIVIANI_POLES x Z

PAIRS = [
    (
        "steinmetz_perpendicular",
        "two unit cylinders with perpendicular crossing axes (z and x through the origin, "
        "length 4): the ellipses x = +-z, x^2 + y^2 = 1, crossing at the nodes (0, +-1, 0) "
        "where the normals are parallel",
        cylinder([0, 0, -2], X, Y, Z, 1.0, 0.0, 4.0),
        cylinder([-2, 0, 0], Y, Z, X, 1.0, 0.0, 4.0),
    ),
    (
        "steinmetz_oblique",
        "two unit cylinders whose axes cross at 60 degrees (z, and (sin 60, 0, cos 60), "
        "through the origin, length 4): two ellipses in the planes p.(a -+ b) = 0, "
        "crossing at the nodes (0, +-1, 0) where the normals are parallel",
        cylinder([0, 0, -2], X, Y, Z, 1.0, 0.0, 4.0),
        cylinder(mul(-2.0, B60), Y, B60_PERP, B60, 1.0, 0.0, 4.0),
    ),
    (
        "torus_bitangent_plane",
        "a torus (axis z, centre the origin, R 2, r 1) and its bitangent plane through the "
        "x axis, 30 degrees off the xy plane: the Villarceau circles of radius 2 about "
        "(+-1, 0, 0) in the plane, crossing at the nodes +-(0, 1.5, sqrt(3)/2) where the "
        "plane touches the torus",
        torus(O, X, Y, Z, 2.0, 1.0),
        patch(O, X, D30, -3.5, 3.5, -3.5, 3.5),
    ),
    (
        "saddle_xy_plane",
        "z = x y over [-1, 1.5] x [-0.7, 1] on the plane z = 0: the lines x = 0 and y = 0 "
        "crossing at the node (0, 0, 0) where the saddle's tangent plane is z = 0",
        saddle(-1.0, 1.5, -0.7, 1.0),
        patch(O, X, Y, -2.0, 2.0, -2.0, 2.0),
    ),
    (
        "viviani_sphere_cylinder",
        "a sphere of radius 2 about the origin (poles on (sin 20, cos 20, 0)) and a unit "
        "cylinder through its centre (axis z through (1, 0, 0), z in [-2.5, 2.5]): "
        "Viviani's figure eight (1 + cos t, sin t, 2 sin(t/2)), t in [0, 4 pi], crossing "
        "itself at the node (2, 0, 0) where the cylinder touches the sphere",
        sphere(O, VIVIANI_E1, VIVIANI_E2, VIVIANI_POLES, 2.0),
        cylinder([1, 0, -2.5], X, Y, Z, 1.0, 0.0, 5.0),
    ),
    (
        "cone_plane_apex",
        "a double cone x^2 + y^2 = z^2 (axis z, z in [-1, 1], apex inside the surface) and "
        "the plane x = z/2 through its apex: the lines y = +-(sqrt(3)/2) z crossing at the "
        "apex (0, 0, 0), where the cone has no normal",
        double_cone(X, Y, Z, 1.0, 1.0),
        patch(O, Y, [0.5, 0.0, 1.0], -1.5, 1.5, -1.5, 1.5),
    ),
    (
        "cones_shared_apex",
        "the double cones x^2 + y^2 = z^2 (axis z, z in [-1, 1]) and y^2 + z^2 = 3 x^2 "
        "(axis x, x in [-1, 1]) sharing their apex: the four lines (t, +-t, +-sqrt(2) t), "
        "|t| <= 1/sqrt(2), crossing at the apex (0, 0, 0), where neither cone has a normal",
        double_cone(X, Y, Z, 1.0, 1.0),
        double_cone(Y, Z, X, math.sqrt(3.0), 1.0),
    ),
]


def case(name, ref, a, b):
    return {
        "name": name,
        "upstream_ref": ref,
        "inputs": {
            "this": {"srf": a},
            "op": "gssi",
            "dom": domain(a),
            "other": {"srf": b},
            "other_dom": domain(b),
            "tol": TOL,
        },
    }


def main():
    cases = []
    for name, what, a, b in PAIRS:
        ref = f"issue #393 (WR-77) singular crossing corpus: {what}"
        cases.append(case(f"sx_{name}", ref, a, b))
        cases.append(case(f"sx_{name}-swapped", ref + "; operands swapped", b, a))
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump({"cases": cases}, f, indent=1)
        f.write("\n")
    print(f"{len(cases)} cases -> {os.path.relpath(OUT, ROOT)}")


if __name__ == "__main__":
    main()
