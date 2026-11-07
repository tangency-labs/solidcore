//! Audit the fourteen frozen native pair records against their actual definitions.
//! The residual module is an import-only copy of the revision being audited.
#[path = "projected_residual_audit_copy.rs"]
mod residual;
use serde_json::{json, Value};
use solidcore::blend::fillet_section::{
    certificate::nurbs_partials_over_piece,
    interval::Ivl,
    station::{as_b_spline, Nurbs},
};
use solidcore::curve::curve::Curve;
use solidcore::geom::{axis2_placement::Axis2Placement, extent2d::Extent2d};
use solidcore::math::{Vec2, Vec3};
use solidcore::nurbs::{
    b_spline_curve::BSplineCurve,
    crv_basic::Curve as NlCurve,
    math::{CPoint, NOW},
};
use solidcore::surface::{
    cone::Cone,
    sphere::{Sphere, SphereOptions},
    surface::Surface,
    torus::{Torus, TorusOptions},
};
use solidcore::Kernel;
use std::io::{self, Write};
type AuditResult<T> = Result<T, Box<dyn std::error::Error>>;
// These frozen canonical nets are continuous. The singleton uses the native
// right-hand span, except at the last natural endpoint; no two-sided hull.
fn knot_cuts(
    knots: &[f64],
    p: usize,
    n: usize,
    lo: f64,
    hi: f64,
    _continuity: usize,
) -> solidcore::Result<Vec<(f64, f64, usize)>> {
    assert!(lo >= knots[p] && hi <= knots[n] && lo <= hi);
    if lo == hi {
        let span = if lo == knots[n] {
            n - 1
        } else {
            (p..n)
                .find(|&i| knots[i] <= lo && lo < knots[i + 1])
                .unwrap()
        };
        return Ok(vec![(lo, hi, span)]);
    }
    Ok((p..n)
        .filter_map(|i| {
            let a = lo.max(knots[i]);
            let b = hi.min(knots[i + 1]);
            (a < b).then_some((a, b, i))
        })
        .collect())
}
fn number(v: &Value) -> AuditResult<f64> {
    Ok(v.as_f64().ok_or("missing number")?)
}
fn vector(v: &Value) -> AuditResult<Vec3> {
    let a = v.as_array().ok_or("missing vector")?;
    if a.len() != 3 {
        return Err("wrong vector length".into());
    }
    Ok(Vec3::new(number(&a[0])?, number(&a[1])?, number(&a[2])?))
}
fn surface(k: &Kernel, s: &Value) -> AuditResult<Surface> {
    let x = vector(&s["x"])?;
    let y = vector(&s["y"])?;
    let z = x.cross(y);
    // Matches the adapter's unit_cross for the frozen orthonormal placements.
    if z.length_squared() != 1.0 {
        return Err("audit requires frozen unit placement axes".into());
    }
    Ok(match s["kind"].as_str().ok_or("surface kind")? {
        "cylinder" => {
            let h0 = number(&s["height"][0])?;
            let h1 = number(&s["height"][1])?;
            let mut a = Cone::new_cylinder_bounded(
                vector(&s["origin"])? + z * h0,
                x,
                y,
                number(&s["radius"])?,
                number(&s["angle"][0])?,
                number(&s["angle"][1])?,
                h1 - h0,
                false,
                false,
                false,
            );
            a.make_nurb(k)?;
            Surface::from(a)
        }
        "cone" => {
            let h0 = number(&s["height"][0])?;
            let h1 = number(&s["height"][1])?;
            let tangent = (number(&s["half_angle_deg"])? * std::f64::consts::PI / 180.0).tan();
            Surface::from(Cone::create_canonical(
                k,
                &Axis2Placement::new(vector(&s["apex"])? + z * h0, x, y),
                h0 * tangent,
                h1 * tangent,
                h1 - h0,
                false,
                false,
            )?)
        }
        "sphere" => {
            let mut a = Sphere::new(
                k,
                vector(&s["center"])?,
                x,
                y,
                Extent2d::from_bounds(
                    number(&s["lon"][0])?,
                    number(&s["lat"][0])?,
                    number(&s["lon"][1])?,
                    number(&s["lat"][1])?,
                ),
                number(&s["radius"])?,
                SphereOptions::default(),
            )?;
            a.make_nurb(k)?;
            Surface::from(a)
        }
        "torus" => {
            let mut a = Torus::new(
                k,
                vector(&s["center"])?,
                x,
                y,
                Extent2d::from_bounds(
                    number(&s["major"][0])?,
                    number(&s["minor"][0])?,
                    number(&s["major"][1])?,
                    number(&s["minor"][1])?,
                ),
                number(&s["major_radius"])?,
                number(&s["minor_radius"])?,
                TorusOptions::default(),
            )?;
            a.make_nurb(k)?;
            Surface::from(a)
        }
        _ => return Err("unsupported frozen surface kind".into()),
    })
}
fn endpoint(x: f64) -> [Ivl; 3] {
    // Same defining-control division as affine_point_definition_enclosure:
    // the original nonrational degree-one curve stores NOW, interpreted as 1.
    [x, 0., 0.].map(|a| Ivl::from(a) / Ivl::from(1.))
}
fn cross(a: [Ivl; 3], b: [Ivl; 3]) -> [Ivl; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn own_metric(net: &Nurbs, uv: Vec2, q: [Ivl; 3]) -> solidcore::Result<Ivl> {
    let u = knot_cuts(&net.u, net.p, net.nu, uv.x, uv.x, 0)?[0].2;
    let v = knot_cuts(&net.v, net.q, net.nv, uv.y, uv.y, 0)?[0].2;
    let a = nurbs_partials_over_piece(net, u, v, uv.x, uv.x, uv.y, uv.y, 0)?.at(0, 0);
    let n = [0., 0., 1.].map(Ivl::from);
    let d = cross(cross(n, [a.x - q[0], a.y - q[1], a.z - q[2]]), n);
    let square = |v: Ivl| Ivl::checked(v.mig() * v.mig(), v.mag() * v.mag());
    Ok(((square(d[0]) + square(d[1])) + square(d[2])).sqrt())
}
fn audit(row: &Value, contract: &Value) -> AuditResult<Value> {
    if contract["adapter_source_sha256"]
        != "ab07e2085d02e1025ea981d9456b54fd819e0feb42b32b0471752f7eeb83dcc8"
        || vector(&contract["curve_controls"][0])? != Vec3::new(-1.0, 0.0, 0.0)
        || vector(&contract["curve_controls"][1])? != Vec3::new(1.0, 0.0, 0.0)
        || vector(&contract["normal"])? != Vec3::Z
        || contract["curve_degree"].as_u64() != Some(1)
        || contract["curve_knots"] != json!([0.0, 0.0, 1.0, 1.0])
        || contract["curve_interval"] != json!([0.0, 1.0])
        || contract["best_answer_so_far_sq"].as_f64() != Some(0.0)
    {
        return Err("unverified static adapter geometry contract".into());
    }
    let i = &row["input"];
    if i["included"] != true
        || i["op"] != "static"
        || i["operation"] != 15
        || i["request"] != 0
        || !matches!(i["entry"].as_i64(), Some(6 | 7))
    {
        return Err("not a focal native pair".into());
    }
    if row["output"]["values"]["status"] != json!([0.0]) {
        return Ok(json!({"state":"solver_refused","checks_performed":false}));
    }
    let answers = row["solutions"]
        .as_array()
        .ok_or("missing decoded answers")?;
    if answers.len() != 1 {
        return Err("expected one retained Single answer".into());
    }
    let a = &answers[0];
    let p = a["parameters"].as_array().ok_or("missing parameters")?;
    if p.len() != 3
        || a["num_objects"].as_f64() != Some(2.0)
        || a["num_variables"].as_f64() != Some(3.0)
    {
        return Err("wrong pair layout".into());
    }
    let swapped = i["swap"].as_bool().unwrap_or(false);
    let (t, uv) = if swapped {
        (number(&p[2])?, Vec2::new(number(&p[0])?, number(&p[1])?))
    } else {
        (number(&p[0])?, Vec2::new(number(&p[1])?, number(&p[2])?))
    };
    let k = Kernel::default();
    let surface = surface(&k, &i["faces"][0])?;
    let domain = surface.natural_uv_domain_with_nurb(&k)?;
    let admissible = (0.0..=1.0).contains(&t) && domain.contains_point(uv, 0.0);
    if !admissible {
        return Ok(json!({"admissible":false,"parameters":[t,uv.x,uv.y]}));
    }
    let curve = Curve::from(BSplineCurve::from_nurb(
        NlCurve::new(
            vec![CPoint::new(-1., 0., 0., NOW), CPoint::new(1., 0., 0., NOW)],
            1,
            vec![0., 0., 1., 1.],
        )?,
        3,
        false,
    ));
    let cp = curve.evaluate_point(&k, t)?;
    let sp = surface.evaluate_point(&k, uv)?;
    let measured = Vec3::Z.cross(cp - sp).cross(Vec3::Z).length();
    let reported = number(&a["value"])?;
    if t != 0.0 && t != 1.0 {
        return Err("focal endpoint audit received an interior curve parameter".into());
    }
    let net = Nurbs::of(&as_b_spline(&k, &surface)?.ok_or("missing defining net")?)?;
    let own = own_metric(&net, uv, endpoint(if t == 0.0 { -1.0 } else { 1.0 }))?;
    let left = residual::maximum_is_bounded_from_definition(
        &k,
        &net,
        &domain,
        endpoint(-1.),
        Vec3::Z,
        own.hi,
    )?;
    let right = residual::maximum_is_bounded_from_definition(
        &k,
        &net,
        &domain,
        endpoint(1.),
        Vec3::Z,
        own.hi,
    )?;
    if !measured.is_finite() || measured <= 0.0 || !own.lo.is_finite() || own.lo <= 0.0 {
        return Err("expected a finite positive focal metric".into());
    }
    // Strictly below the rigorous native witness lower bound, not merely
    // below a rounded scalar evaluation.
    let reduced = f64::from_bits(own.lo.to_bits() - 1);
    let reduced_proved = residual::maximum_is_bounded_from_definition(
        &k,
        &net,
        &domain,
        endpoint(if t == 0.0 { -1.0 } else { 1.0 }),
        Vec3::Z,
        reduced,
    )?;
    Ok(json!({"admissible":true,"parameters":[t,uv.x,uv.y],
        "curve_point":[cp.x,cp.y,cp.z],"surface_point":[sp.x,sp.y,sp.z],
        "native_metric":measured,"reported":reported,"report_error":reported-measured,
        "own_enclosure":[own.lo,own.hi],"report_in_own_enclosure":own.lo<=reported && reported<=own.hi,
        "native_metric_in_own_enclosure":own.lo<=measured && measured<=own.hi,
        "left_domain_proved":left,"right_domain_proved":right,
        "reduced_threshold":reduced,"reduced_threshold_basis":"own_enclosure_lower","reduced_threshold_refused":!reduced_proved}))
}
fn main() -> AuditResult<()> {
    let file = std::env::args()
        .nth(1)
        .ok_or("usage: native-audit OUTPUTS.jsonl")?;
    let contract_file = std::env::args()
        .nth(2)
        .ok_or("missing adapter-contract.json")?;
    let contract: Value = serde_json::from_str(&std::fs::read_to_string(contract_file)?)?;
    let mut out = io::BufWriter::new(io::stdout().lock());
    for line in std::fs::read_to_string(file)?.lines() {
        let row: Value = serde_json::from_str(line)?;
        if row["focal_native_pair"] != true {
            continue;
        }
        let result = match audit(&row, &contract) {
            Ok(v) => json!({"case":row["case"],"native":v}),
            Err(e) => json!({"case":row["case"],"native_audit_error":e.to_string()}),
        };
        writeln!(out, "{result}")?;
        out.flush()?;
    }
    Ok(())
}
