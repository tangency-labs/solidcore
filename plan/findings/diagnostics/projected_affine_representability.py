"""Exact parameter arithmetic for the retained narrow-span regression inputs.

On the central segment C(s)=(2s-1,s), squared distance from (0,y) is
5(s-(2+y)/5)^2+(1-2y)^2/5. The two binary64 parameters bracketing the
exact foot therefore include the best representable parameter globally.
Fractions retain every input bit; Decimal is used only to print square roots.
"""
from decimal import Decimal, localcontext
from fractions import Fraction as F
import json
import math

for gap, scale in [(1.5e-12, 1.0), (1.9e-12, 1.0), (1.5e-12, 1e6),
                   (1e-8, 1.0), (1e-8, 1e-3), (1e-8, 1e6)]:
    for y in [0.2, 0.3]:
        lo, hi, yy = F(scale), F((1.0 + gap) * scale), F(y)
        exact_parameter = lo + (hi-lo) * (2+yy)/5
        rounded = float(exact_parameter)
        candidates = [math.nextafter(rounded, -math.inf), rounded,
                      math.nextafter(rounded, math.inf)]
        assert F(candidates[0]) <= exact_parameter <= F(candidates[-1])
        squared = []
        for t in candidates:
            s = (F(t)-lo)/(hi-lo)
            squared.append((2*s-1)**2 + (s-yy)**2)
        exact_minimum_squared = (1-2*yy)**2/5
        best = min(squared)
        assert best >= exact_minimum_squared
        with localcontext() as context:
            context.prec = 80
            def decimal(value):
                return Decimal(value.numerator)/Decimal(value.denominator)
            excess = decimal(best).sqrt()-decimal(exact_minimum_squared).sqrt()
        print(json.dumps({"gap": gap, "parameter_scale": scale, "query_y": y,
                          "best_representable_excess": str(excess)}))

# Cross-leaf maximum: the defining binary64 normal is not the real number 1/10.
# Its large-coordinate product rounds away in the old reported scalar, although
# its exact projected distance is larger than the small, stably evaluated leaf.
normal = F(0.1)
large_cross = F(1e8)*normal-F(1e7)
large_squared = (1+normal*normal)*large_cross*large_cross
small_squared = (1+normal*normal)*F(1e-10)*F(1e-10)
assert large_cross == F(1e8)/(5*2**55)
assert large_squared > small_squared > 0
with localcontext() as context:
    context.prec = 80
    large = (Decimal(large_squared.numerator)/Decimal(large_squared.denominator)).sqrt()
    small = (Decimal(small_squared.numerator)/Decimal(small_squared.denominator)).sqrt()
print(json.dumps({"case": "cross_leaf_maximum", "true_first_value": str(large),
                  "true_second_value": str(small)}))
