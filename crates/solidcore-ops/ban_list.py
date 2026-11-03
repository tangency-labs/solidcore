"""Every public method of the kernel's Vec3/Point3, Vec2/Point2 and Transform impl blocks, as
clippy disallowed-methods entries, less the ones that only copy fields (KEEP).
usage: python3 ban_list.py <solidcore crate src dir> <reason>"""
import re, sys
src, reason = sys.argv[1], sys.argv[2]
KEEP = {"new", "from_xy", "swap_xy", "x_axis", "y_axis", "z_axis", "translation", "translation_part", "from_array", "to_array", "component", "component_mut", "set_uninitialized", "is_initialized", "is_undef", "identity"}
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
for path, ty, name in [(src + "/geom/vector3d.rs", "Vec3", "solidcore::Vec3"), (src + "/geom/vector2d.rs", "Vec2", "solidcore::Vec2"), (src + "/math.rs", "Transform", "solidcore::Transform")]:
    for f in sorted(set(methods(path, ty)) - KEEP):
        rows.append(f'    {{ path = "{name}::{f}", reason = "{reason}" }},')
print("\n".join(rows))
