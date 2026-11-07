"""Emit a Rust audit of changed analytic-face records from two digest-debug files.

Run from the repository root, passing before.debug and after.debug. Compile the
emitted Rust against solidcore. CSV columns are case, before/after, scalar error,
UVs in natural rectangular domains, line parameter in [0,1], and the directed
query's transverse residual (zero for nondirected queries). Domain checks
use zero padding. They do not certify arbitrary trim loops, completeness,
stationarity, or global optimality, except for the explicitly checked cylinder
curve/edge maxima. The final column reports that closed-form check or 'not_checked'.
"""

import argparse
import json
import re
from pathlib import Path


def number(value):
    """Write a fixture number as a round-trippable Rust floating-point literal."""
    return repr(float(value))


def vector(value):
    """Emit a Vec3 constructor preserving the fixture's coordinate order."""
    return "Vec3::new(" + ",".join(map(number, value)) + ")"


def surface(spec):
    """Rebuild the analytic fixture using the runner's constructors and UV maps.

    Converting the same surface kind to NURBS is essential: returned parameters
    refer to that rational form, not to angles in the analytic constructor.
    """
    kind = spec["kind"]
    x, y = vector(spec["x"]), vector(spec["y"])
    if kind == "cylinder":
        origin, radius = vector(spec["origin"]), number(spec["radius"])
        lo, hi = map(number, spec["height"])
        a, b = map(number, spec["angle"])
        return f"{{ let z = {x}.cross({y}).unitize()?; let mut s = Cone::new_cylinder_bounded({origin}+z*{lo}, {x}, {y}, {radius}, {a}, {b}, {hi}-{lo}, false, false, false); s.make_nurb(&k)?; Surface::from(s) }}"
    if kind == "cone":
        origin = vector(spec["apex"])
        lo, hi = map(number, spec["height"])
        angle = number(spec["half_angle_deg"])
        return f"{{ let z = {x}.cross({y}).unitize()?; let t = ({angle}_f64 * std::f64::consts::PI / 180.0).tan(); Surface::from(Cone::create_canonical(&k, &Axis2Placement::new({origin}+z*{lo}, {x}, {y}), {lo}*t, {hi}*t, {hi}-{lo}, false, false)?) }}"
    center = vector(spec["center"])
    if kind == "sphere":
        u0, u1 = map(number, spec["lon"])
        v0, v1 = map(number, spec["lat"])
        radii = number(spec["radius"])
        constructor, options = "Sphere", "SphereOptions"
    elif kind == "torus":
        u0, u1 = map(number, spec["major"])
        v0, v1 = map(number, spec["minor"])
        radii = number(spec["major_radius"]) + "," + number(spec["minor_radius"])
        constructor, options = "Torus", "TorusOptions"
    else:
        raise ValueError(f"unsupported analytic kind: {kind}")
    return f"{{ let mut s = {constructor}::new(&k, {center}, {x}, {y}, Extent2d::from_bounds({u0},{v0},{u1},{v1}), {radii}, {options}::default())?; s.make_nurb(&k)?; Surface::from(s) }}"


def records(path):
    """Index the runner's tab-separated digest-debug records by fixture name."""
    return dict(line.split("\t", 1) for line in Path(path).read_text().splitlines())


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("before", help="digest-debug file before the change")
parser.add_argument("after", help="digest-debug file after the change")
args = parser.parse_args()
before, after = records(args.before), records(args.after)
fixtures = {}
for path in Path("fixtures/intersect.topology_solver_p2").glob("*.json"):
    for case in json.loads(path.read_text()).get("cases", []):
        fixtures[case["name"]] = case["inputs"]

print("""
use solidcore::{kernel::Kernel, math::Vec3};
use solidcore::geom::{axis2_placement::Axis2Placement, extent2d::Extent2d, vector2d::Vec2};
use solidcore::surface::{cone::Cone, sphere::{Sphere, SphereOptions}, torus::{Torus, TorusOptions}, surface::Surface};
fn in_domain(k: &Kernel, surface: &Surface, uv: Vec2) -> solidcore::Result<bool> {
    Ok(surface.natural_uv_domain_with_nurb(k)?.contains_point(uv, 0.0))
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
let k = Kernel::default();
""")
for name in sorted(before.keys() & after.keys()):
    if before[name] == after[name]:
        continue
    case = fixtures.get(name)
    if (case is None or not name.startswith("analytic_face_")
            or case.get("entry") not in (5, 6, 7, 8) or not case.get("faces")):
        continue
    for label, record in [("before", before[name]), ("after", after[name])]:
        match = re.search(r'"static": (\[[^]]*\])', record)
        values = json.loads(match.group(1))
        assert values[5] == 1, "audit expects one returned solution"
        scalar = number(values[11])
        params = values[12 : 12 + int(values[7])]
        print("{")
        print("let first = " + surface(case["faces"][0]) + ";")
        if case["entry"] == 5:
            print("let a = " + vector(case["first"][0]) + ";")
            u, v = map(number, params)
            print(f"let b = first.evaluate_point(&k, Vec2::new({u},{v}))?;")
            print(f"let uv_inside = in_domain(&k, &first, Vec2::new({u},{v}))?; let line_inside = true;")
        elif case["entry"] in (6, 7):
            t, u, v = map(number, params)
            print(f"let a = Vec3::new(-1.0+2.0*{t},0.0,0.0); let b = first.evaluate_point(&k, Vec2::new({u},{v}))?;")
            print(f"let uv_inside = in_domain(&k, &first, Vec2::new({u},{v}))?; let line_inside = (0.0..=1.0).contains(&({t}_f64));")
        elif case["entry"] == 8:
            print("let second = " + surface(case["faces"][1]) + ";")
            u, v, r, s = map(number, params)
            print(f"let a = first.evaluate_point(&k,Vec2::new({u},{v}))?; let b = second.evaluate_point(&k,Vec2::new({r},{s}))?;")
            print(f"let uv_inside = in_domain(&k, &first, Vec2::new({u},{v}))? && in_domain(&k, &second, Vec2::new({r},{s}))?; let line_inside = true;")
        else:
            raise ValueError("unsupported entry")
        print(f"let measured = ((a.x-b.x).powi(2)+(a.y-b.y).powi(2)).sqrt(); let error = (({scalar}_f64).abs()-measured).abs();")
        # The static runner projects along Z and directs along X, so its
        # directed candidates must have equal projected y coordinates.
        directed = case["operation"] in (11, 12, 13)
        print("let transverse = " + ("(a.y-b.y).abs()" if directed else "0.0_f64") + ";")
        maximality = '"not_checked"'
        if name in {
            "analytic_face_cyl_x_6_15", "analytic_face_cyl_x_7_15",
            "analytic_face_halfcyl_x_6_15", "analytic_face_halfcyl_x_7_15",
        }:
            # These fixtures project to [-.5,.5]^2 against [-1,1] x {0}.
            # Squared distance is convex in each coordinate: its maximum is
            # attained at opposite x endpoints and y=+/-.5, hence sqrt(2.5).
            spec = case["faces"][0]
            assert spec["kind"] == "cylinder" and spec["radius"] == 0.5
            assert spec["origin"] == [-0.5, 0, 0] and spec["height"] == [0, 1]
            assert spec["x"] == [0, 1, 0] and spec["y"] == [0, 0, 1]
            assert spec["angle"] in ([0, 180], [0, 360])
            print(f"let is_maximum = uv_inside && line_inside && error <= {number(case['tolerance'])} && (measured - 2.5_f64.sqrt()).abs() <= {number(case['tolerance'])};")
            maximality = 'if is_maximum { "true" } else { "false" }'
        print(f'println!("{name},{label},{{error:.17e}},{{uv_inside}},{{line_inside}},{{transverse:.17e}},{{}}", {maximality}); }}')
print("Ok(()) }")
