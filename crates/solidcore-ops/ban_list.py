"""The disallowed-methods list of the layers that compute no geometry (docs/layers.md §4): the
float maths, and every public method of the kernel's vector, transform and matrix impl blocks
less the ones that only copy fields (KEEP). Prints the whole `disallowed-methods = [...]` array.
usage: python3 ban_list.py <solidcore crate src dir> <reason> [--floats-only]
(--floats-only for solidcore-record, which does not depend on the kernel.)"""
import re, sys
src, reason = sys.argv[1], sys.argv[2]
floats_only = "--floats-only" in sys.argv[3:]
FLOATS = ["sqrt", "cbrt", "hypot", "sin", "cos", "tan", "sin_cos", "asin", "acos", "atan", "atan2",
          "sinh", "cosh", "tanh", "asinh", "acosh", "atanh", "powi", "powf", "exp", "exp2", "exp_m1",
          "ln", "ln_1p", "log", "log2", "log10", "mul_add", "recip"]
KEEP = {"new", "try_new", "from_xy", "swap_xy", "x_axis", "y_axis", "z_axis", "translation", "translation_part",
        "from_array", "to_array", "component", "component_mut", "set_uninitialized", "is_initialized", "is_undef",
        "identity", "from_row_major", "num_rows", "num_columns", "as_slice", "at", "set_at", "row", "row_mut",
        "set_size", "clear", "matrix", "type_code", "type_string", "class_string", "is_kind_of", "read_from_db"}
def methods(path, ty):
    s = open(path).read()
    out = []
    for m in re.finditer(r'^impl ' + ty + r' \{', s, re.M):
        i, depth = m.end(), 1
        while depth:
            c = s[i]; depth += (c == '{') - (c == '}'); i += 1
        out += re.findall(r'^    pub (?:const )?fn (\w+)', s[m.end():i], re.M)
    return out
rows = []
for t in ("f64", "f32"):
    rows += [f'    {{ path = "{t}::{f}", reason = "{reason}" }},' for f in FLOATS]
if not floats_only:
    for path, ty, name in [(src + "/geom/vector3d.rs", "Vec3", "solidcore::Vec3"),
                           (src + "/geom/vector2d.rs", "Vec2", "solidcore::Vec2"),
                           (src + "/math.rs", "Transform", "solidcore::Transform"),
                           (src + "/geom/transform.rs", "Transform", "solidcore::geom::transform::Transform"),
                           (src + "/geom/matrix.rs", "Matrix", "solidcore::geom::matrix::Matrix")]:
        for f in sorted(set(methods(path, ty)) - KEEP):
            rows.append(f'    {{ path = "{name}::{f}", reason = "{reason}" }},')
print("disallowed-methods = [\n" + "\n".join(rows) + "\n]")
