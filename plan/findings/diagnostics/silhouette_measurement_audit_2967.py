"""Emit a Rust audit of changed analytic-face records from two digest-debug files.

Run from the repository root, passing before.debug and after.debug. Compile the
emitted Rust against solidcore. It checks scalar/returned-point consistency,
not completeness, stationarity, trimming membership, or global optimality.
"""

import json
import re
import sys
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


before, after = map(records, sys.argv[1:3])
fixtures = {}
for path in Path("fixtures/intersect.topology_solver_p2").glob("*.json"):
    for case in json.loads(path.read_text()).get("cases", []):
        fixtures[case["name"]] = case["inputs"]

print("""
use solidcore::{kernel::Kernel, math::Vec3};
use solidcore::geom::{axis2_placement::Axis2Placement, extent2d::Extent2d, vector2d::Vec2};
use solidcore::surface::{cone::Cone, sphere::{Sphere, SphereOptions}, torus::{Torus, TorusOptions}, surface::Surface};
fn main() -> Result<(), Box<dyn std::error::Error>> {
let k = Kernel::default();
""")
for name in sorted(before.keys() & after.keys()):
    if before[name] == after[name]:
        continue
    case = fixtures[name]
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
        elif case["entry"] in (6, 7):
            t, u, v = map(number, params)
            print(f"let a = Vec3::new(-1.0+2.0*{t},0.0,0.0); let b = first.evaluate_point(&k, Vec2::new({u},{v}))?;")
        elif case["entry"] == 8:
            print("let second = " + surface(case["faces"][1]) + ";")
            u, v, r, s = map(number, params)
            print(f"let a = first.evaluate_point(&k,Vec2::new({u},{v}))?; let b = second.evaluate_point(&k,Vec2::new({r},{s}))?;")
        else:
            raise ValueError("unsupported entry")
        print(f"let measured = ((a.x-b.x).powi(2)+(a.y-b.y).powi(2)).sqrt(); let error = (({scalar}_f64).abs()-measured).abs();")
        print(f'println!("{name},{label},{{error:.17e}}"); }}')
print("Ok(()) }")
