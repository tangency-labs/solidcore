"""Verify WR-449's recorded point-to-whole-edge bound using exact arithmetic.

Run with Python 3.10+: python3 scripts/verify_wr449_boundary.py
The input doubles are observations, not interval-certified surface evaluations.
"""

from fractions import Fraction
import json
from pathlib import Path


def exact(value):
    return Fraction.from_float(float.fromhex(value))


def split(controls):
    a, b, c = controls
    ab = tuple((x + y) / 2 for x, y in zip(a, b))
    bc = tuple((x + y) / 2 for x, y in zip(b, c))
    middle = tuple((x + y) / 2 for x, y in zip(ab, bc))
    return (a, ab, middle), (middle, bc, c)


def box_distance_squared(point, controls):
    projected = [tuple(v[j] / v[3] for j in range(3)) for v in controls]
    distance = Fraction(0)
    for axis, coordinate in enumerate(point):
        low = min(v[axis] for v in projected)
        high = max(v[axis] for v in projected)
        gap = max(low - coordinate, coordinate - high, Fraction(0))
        distance += gap * gap
    return distance


def verify(data):
    point = tuple(map(exact, data["point_double_hex"]))
    controls = [tuple(map(exact, row)) for row in data["homogeneous_double_hex"]]
    knots = tuple(map(exact, data["knots"]))
    threshold = Fraction(data["whole_edge_lower_bound"])
    assert len(point) == 3 and len(controls) == 5
    assert all(len(row) == 4 and row[3] > 0 for row in controls)
    assert data["degree"] == 2 and data["knot_multiplicities"] == [3, 2, 3]
    assert len(knots) == 3 and knots[0] < knots[1] < knots[2]
    assert tuple(map(exact, data["edge_interval_double_hex"])) == (knots[0], knots[-1])
    assert threshold == Fraction(3, 50)
    assert threshold > exact(data["postmass_edge_tolerance_hex"])

    # Multiplicity=degree at the internal knot gives two complete Bezier spans.
    # Positive weights keep each rational curve inside its projected control hull.
    pending = [(controls[:3], 0), (controls[2:], 0)]
    leaves = 0
    maximum_depth = 0
    while pending:
        segment, depth = pending.pop()
        if box_distance_squared(point, segment) > threshold * threshold:
            leaves += 1
            maximum_depth = max(maximum_depth, depth)
            continue
        if depth >= 24:
            raise ValueError("Could not certify the bound within 24 subdivisions")
        left, right = split(segment)
        # Replacing a segment with both children preserves whole-domain coverage.
        pending.extend(((right, depth + 1), (left, depth + 1)))
    return leaves, maximum_depth


if __name__ == "__main__":
    source = Path(__file__).parent / "evidence" / "wr449-boundary.json"
    leaves, depth = verify(json.loads(source.read_text()))
    print(f"PASS: whole-edge distance > 0.06; {leaves} boxes, maximum depth {depth}")
