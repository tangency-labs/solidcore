//! Verify proposed #3060 witnesses through the actual public Rust definitions.
//! Link against solidcore and serde_json rlibs of the audited build.
//! Input: projected_pair_unchanged_all_witnesses_3060.jsonl. No optimization,
//! golden writing or solver exclusions. Constructors mirror curve_curve_p9.
use serde_json::{json, Value};
use solidcore::curve::curve::Curve;
use solidcore::curve::ellipse::{Ellipse, EllipseOptions};
use solidcore::curve::hermite_curve::HermiteCurve;
use solidcore::curve::line::{BoundedLineOptions, Line};
use solidcore::geom::axis2_placement::Axis2Placement;
use solidcore::geom::extent1d::Extent1d;
use solidcore::kernel::Kernel;
use solidcore::math::Vec3;
use solidcore::nurbs::b_spline_curve::{BSplineCurve, NurbCircleParam};
use solidcore::nurbs::nurbs_crv::allocate_nurb_curve;
use std::error::Error;

type Result<T> = std::result::Result<T, Box<dyn Error>>;
fn number(v: &Value) -> Result<f64> {
    v.as_f64().ok_or_else(|| "missing number".into())
}
fn vector(v: &Value) -> Result<Vec3> {
    Ok(Vec3::new(number(&v[0])?, number(&v[1])?, number(&v[2])?))
}
fn interval(v: &Value) -> Result<Extent1d> {
    Ok(Extent1d::from_min_max(number(&v[0])?, number(&v[1])?))
}
fn build(k: &Kernel, spec: &Value) -> Result<Curve> {
    let map = spec.as_object().ok_or("curve must be an object")?;
    let (kind, d) = map.iter().next().ok_or("empty curve")?;
    Ok(match kind.as_str() {
        "line" => Curve::from(Line::bounded(
            vector(&d["p"])?,
            vector(&d["v"])?,
            &interval(&d["domain"])?,
            &BoundedLineOptions {
                scale: d["scale"].as_f64().unwrap_or(1.0),
                inside_out: d["inside_out"].as_bool().unwrap_or(false),
                ..Default::default()
            },
        )),
        "ellipse" => Curve::from(Ellipse::new(
            k,
            vector(&d["c"])?,
            vector(&d["x"])?,
            vector(&d["y"])?,
            interval(&d["domain"])?,
            number(&d["radii"][0])?,
            number(&d["radii"][1])?,
            &EllipseOptions {
                dim: 3,
                nurb: None,
                inside_out: d["inside_out"].as_bool().unwrap_or(false),
            },
        )?),
        "hermite" => Curve::from(HermiteCurve::new(
            vector(&d["p1"])?,
            vector(&d["d1"])?,
            vector(&d["p2"])?,
            vector(&d["d2"])?,
            d["dim"].as_u64().ok_or("missing dim")?,
        )),
        "bspline" => {
            let points = d["cpts"].as_array().ok_or("missing controls")?;
            let knots = d["u"].as_array().ok_or("missing knots")?;
            let degree = usize::try_from(d["p"].as_u64().ok_or("missing degree")?)?;
            let mut net = allocate_nurb_curve(
                points.len().checked_sub(1).ok_or("no controls")?,
                degree,
                knots.len().checked_sub(1).ok_or("no knots")?,
            )?;
            for (dest, src) in net.cpts_mut().iter_mut().zip(points) {
                dest.x = number(&src[0])?;
                dest.y = number(&src[1])?;
                dest.z = number(&src[2])?;
                dest.w = number(&src[3])?;
            }
            for (dest, src) in net.knots_mut().iter_mut().zip(knots) {
                *dest = number(src)?;
            }
            let mut curve =
                BSplineCurve::from_nurb(net, d["dim"].as_u64().ok_or("missing dim")?, false);
            curve.set_out_of_bounds_enabled(d["oob"].as_bool().unwrap_or(false));
            Curve::from(curve)
        }
        "circle" => Curve::from(BSplineCurve::create_circle_segment(
            3,
            &Axis2Placement::new(vector(&d["center"])?, vector(&d["x"])?, vector(&d["y"])?),
            number(&d["radius"])?,
            number(&d["a0"])?,
            number(&d["a1"])?,
            NurbCircleParam::Quadratic,
        )?),
        _ => return Err(format!("unsupported curve kind {kind}").into()),
    })
}
fn audit(row: &Value) -> Result<Value> {
    let k = Kernel::default();
    let i = &row["input"];
    if i["op"] != "gcs" || i["sop"] != 15 || i["requested"] != 1 {
        return Err("not original All projected maximum".into());
    }
    let mut a = build(&k, &i["c1"])?;
    if let Some(t) = i.get("transform1").filter(|v| !v.is_null()) {
        a.transform(
            &k,
            &Axis2Placement::new(vector(&t[0])?, vector(&t[1])?, vector(&t[2])?),
            None,
        )?;
    }
    let other = if i["c2"].is_null() {
        None
    } else {
        Some(build(&k, &i["c2"])?)
    };
    let b = other.as_ref().unwrap_or(&a);
    let normal = vector(&i["vecs"][0])?;
    let parameters = &row["witness"]["search_parameters"];
    let s = number(&parameters[0])?;
    let t = number(&parameters[1])?;
    let requested = [interval(&i["ivl1"])?, interval(&i["ivl2"])?];
    let natural = [a.natural_interval(), b.natural_interval()];
    let mut admissible = true;
    let mut domains = Vec::new();
    for ((requested, natural), p) in requested.into_iter().zip(natural).zip([s, t]) {
        let lo = requested.min().max(natural.min());
        let hi = requested.max().min(natural.max());
        admissible &= p.is_finite() && lo <= p && p <= hi;
        domains.push([lo, hi]);
    }
    let pa = a.evaluate_point(&k, s)?;
    let pb = b.evaluate_point(&k, t)?;
    let measured = normal.cross(pa - pb).cross(normal).length();
    let mut previous = Vec::new();
    let mut best = f64::NEG_INFINITY;
    for old in row["witness"]["checked"]
        .as_array()
        .ok_or("missing original solutions")?
    {
        let os = number(&old["parameters"][0])?;
        let ot = number(&old["parameters"][1])?;
        let value = normal
            .cross(a.evaluate_point(&k, os)? - b.evaluate_point(&k, ot)?)
            .cross(normal)
            .length();
        let reported = number(&old["reported"])?;
        best = best.max(reported);
        previous.push(json!({"parameters":[os,ot],"reported":reported,"rust_measured":value,"scalar_error":value-reported}));
    }
    let tolerance = number(&i["tol"])?;
    let target = i["target"].as_f64();
    let meets_target = target.is_none_or(|target| measured >= target.abs() - tolerance);
    let gap = measured - best;
    Ok(
        json!({"case":row["case"],"admissible":admissible,"admissible_domains":domains,
  "parameters":[s,t],"rust_points":[[pa.x,pa.y,pa.z],[pb.x,pb.y,pb.z]],"rust_measured":measured,
  "previous_solutions":previous,"previous_best":best,"tolerance":tolerance,"target":target,
  "meets_incumbent_cutoff":meets_target,"larger_than_best_by":gap,
  "disproves_best_beyond_slack":admissible&&meets_target&&measured.is_finite()&&gap>tolerance}),
    )
}
fn main() -> Result<()> {
    let path = std::env::args().nth(1).ok_or("usage: check INPUT.jsonl")?;
    for line in std::fs::read_to_string(path)?
        .lines()
        .filter(|s| !s.trim().is_empty())
    {
        let row: Value = serde_json::from_str(line.trim_start_matches('\u{feff}'))?;
        match audit(&row) {
            Ok(result) => println!("{result}"),
            Err(e) => println!(
                "{}",
                json!({"case":row["case"],"audit_error":e.to_string()})
            ),
        }
    }
    Ok(())
}
