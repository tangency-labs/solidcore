use serde::{Deserialize, Serialize};
use solidcore::curve::{curve::Curve, hermite_curve::HermiteCurve, offset_curve::OffsetCurve};
use solidcore::geom::extent1d::Extent1d;
use solidcore::nurbs::{b_spline_curve::BSplineCurve, nurbs_crv::allocate_nurb_curve};
use solidcore::{math::Vec3, Kernel};
type HarnessError = String;
fn harness(s: impl Into<String>) -> String {
    s.into()
}
#[derive(Debug, Clone, Serialize, Deserialize)]
struct BSplineIn {
    dim: u64,
    #[serde(default)]
    oob: bool,
    cpts: Vec<[f64; 4]>,
    p: usize,
    u: Vec<f64>,
    /// Knots written back over the curve's own NURB after construction, undoing the
    /// constructor's knot fix-up (base block kind `2`). Omitted when false, so
    /// existing case identities are unchanged.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    raw_knots: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct HermiteIn {
    dim: u64,
    p1: [f64; 3],
    d1: [f64; 3],
    p2: [f64; 3],
    d2: [f64; 3],
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
enum BaseIn {
    Bspline(BSplineIn),
    Hermite(HermiteIn),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct SpecIn {
    #[serde(default = "three")]
    dim: u64,
    a: f64,
    b: f64,
    normal: [f64; 3],
    dist: f64,
    #[serde(default)]
    owns: bool,
    #[serde(default)]
    nrev: usize,
}

const fn three() -> u64 {
    3
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "op", rename_all = "snake_case")]
enum OpIn {
    State,
    Eval {
        t: f64,
        n: usize,
        left: bool,
        nonzero: bool,
    },
    EvalPoint {
        t: f64,
    },
    Knots {
        ivl: Option<[f64; 2]>,
        mults: bool,
    },
    Continuities {
        angle: f64,
    },
    Reverse {
        a: f64,
        b: f64,
        probes: Vec<f64>,
    },
    Trim {
        a: f64,
        b: f64,
        /// 0 `TRUE`, 1 `FALSE`, 2 `UNSURE`.
        mode: u8,
        probes: Vec<f64>,
    },
    Eq {
        #[serde(default)]
        this: bool,
        #[serde(default)]
        copy: bool,
        other: Option<SpecIn>,
    },
    Copy {
        probes: Vec<f64>,
    },
    /// The copy constructor (`new (ctx) SmOffsetCurve(*original)`), not `Copy`.
    CloneCtor {
        probes: Vec<f64>,
    },
    IsKindOf {
        #[serde(rename = "type")]
        type_code: i64,
    },
    /// `AssertValid` at `level` (0..2, 4 = `SM_LEVEL_GIVEN`), with or without a
    /// report list; the reports as `crate::util::assert_reports` writes them.
    AssertValid {
        level: u8,
        list: bool,
    },
    ReadDb {
        binary: bool,
        version: u64,
        #[serde(default)]
        text: Option<String>,
        #[serde(default)]
        bytes: Option<Vec<u8>>,
        probes: Vec<f64>,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct Inputs {
    #[serde(default)]
    base: Option<BaseIn>,
    #[serde(default)]
    spec: Option<SpecIn>,
    #[serde(default)]
    outer: Option<SpecIn>,
    #[serde(flatten)]
    op: OpIn,
}

fn v3(p: [f64; 3]) -> Vec3 {
    Vec3::from_array(p)
}

fn build_base(c: &BaseIn) -> Result<Curve, HarnessError> {
    match c {
        BaseIn::Bspline(bs) => {
            let n = bs
                .cpts
                .len()
                .checked_sub(1)
                .ok_or_else(|| harness("no points"))?;
            let m =
                bs.u.len()
                    .checked_sub(1)
                    .ok_or_else(|| harness("no knots"))?;
            let mut cur =
                allocate_nurb_curve(n, bs.p, m).map_err(|e| harness(format!("curve: {e}")))?;
            for (d, s) in cur.cpts_mut().iter_mut().zip(&bs.cpts) {
                d.x = s[0];
                d.y = s[1];
                d.z = s[2];
                d.w = s[3];
            }
            cur.knots_mut().copy_from_slice(&bs.u);
            let mut c = BSplineCurve::from_nurb(cur, bs.dim, false);
            if bs.raw_knots {
                // undo the constructor's fix-up, as the shim does
                c.get_or_create_gw_nurb()
                    .map_err(|e| harness(format!("curve: {e}")))?
                    .knots_mut()
                    .copy_from_slice(&bs.u);
            }
            c.set_out_of_bounds_enabled(bs.oob);
            Ok(Curve::from(c))
        }
        BaseIn::Hermite(h) => Ok(Curve::from(HermiteCurve::new(
            v3(h.p1),
            v3(h.d1),
            v3(h.p2),
            v3(h.d2),
            h.dim,
        ))),
    }
}

fn build_offset(base: Curve, s: &SpecIn) -> solidcore::Result<Curve> {
    let mut c = Curve::from(OffsetCurve::new(
        s.dim,
        base,
        Extent1d::from_min_max(s.a, s.b),
        v3(s.normal),
        s.dist,
        s.owns,
    )?);
    for _ in 0..s.nrev {
        let ivl = c.natural_interval();
        // (status ignored by the shim; the offset reversal cannot fail)
        let _ = c.reverse_parameterization(&Kernel::default(), &ivl);
    }
    Ok(c)
}

/// Most reversals the runner applies: the `nrev` loop is runner code, so a fuzzed count
/// must not make it run for hours.
const MAX_RUNNER_REVERSALS: usize = 4;

/// The base curve and the op's curve.
fn build(i: &Inputs) -> Result<(Curve, Curve), HarnessError> {
    let other = match &i.op {
        OpIn::Eq { other, .. } => other.as_ref(),
        _ => None,
    };
    for s in [i.spec.as_ref(), i.outer.as_ref(), other]
        .into_iter()
        .flatten()
    {
        if s.nrev > MAX_RUNNER_REVERSALS {
            return Err(harness("nrev above the runner's reversal bound"));
        }
    }
    // the curves are built again rather than cloned: a Rust clone would fix the base
    // curve's knots up, while upstream offsets the curve itself
    let base_in = i.base.as_ref().ok_or_else(|| harness("no base"))?;
    let spec = i.spec.as_ref().ok_or_else(|| harness("no spec"))?;
    let built =
        |base: Curve, spec: &SpecIn| build_offset(base, spec).map_err(|e| harness(e.to_string()));
    let inner = built(build_base(base_in)?, spec)?;
    Ok(match &i.outer {
        None => (build_base(base_in)?, inner),
        Some(o) => (built(build_base(base_in)?, spec)?, built(inner, o)?),
    })
}

fn main() {
    let k = Kernel::default();
    for line in std::fs::read_to_string(std::env::args().nth(1).unwrap())
        .unwrap()
        .lines()
    {
        let r: serde_json::Value = serde_json::from_str(line).unwrap();
        let i: Inputs = serde_json::from_value(r["input"].clone()).unwrap();
        let (base, c) = build(&i).unwrap();
        if std::env::args().nth(2).as_deref() == Some("--plot") {
            if r["case"] == "gen1381" {
                let samples: Vec<_> = (0..=256)
                    .map(|j| {
                        let t = f64::from(j) / 256.0;
                        let p = c.evaluate_point(&k, t).unwrap();
                        [t, p.x, p.y, p.z]
                    })
                    .collect();
                let markers: Vec<_> = c
                    .knots(&k, None)
                    .unwrap()
                    .0
                    .into_iter()
                    .map(|t| {
                        let p = c.evaluate_point(&k, t).unwrap();
                        [t, p.x, p.y, p.z]
                    })
                    .collect();
                println!(
                    "{}",
                    serde_json::json!({"case":r["case"],"native_samples":samples,"rounded_knot_points":markers})
                );
            }
            continue;
        }
        let native = match &i.op {
            OpIn::Knots { ivl, mults } => {
                let interval = ivl.map(|[a, b]| Extent1d::from_min_max(a, b));
                match c.knots(&k, interval.as_ref()) {
                    Ok((ks, ms)) => {
                        let mut out = vec![ks.len() as f64];
                        out.extend(ks);
                        out.push(if *mults { ms.len() as f64 } else { 0. });
                        if *mults {
                            out.extend(ms.into_iter().map(|m| m as f64));
                        }
                        serde_json::json!({"status":[0],"out":out})
                    }
                    Err(e) => serde_json::json!({"status":[1007],"out":[],"reason":e.to_string()}),
                }
            }
            OpIn::Continuities { angle } => {
                let (min, all) = c.calculate_continuities(&k, *angle).unwrap();
                let mut out = vec![min as i32 as f64, all.len() as f64];
                out.extend(all.into_iter().map(|v| v as i32 as f64));
                serde_json::json!({"status":[0],"out":out})
            }
            _ => panic!(),
        };
        let base_knots = base
            .knots(&k, None)
            .map(|(ks, ms)| serde_json::json!({"knots":ks,"multiplicities":ms}))
            .unwrap_or_else(|e| serde_json::json!({"error":e.to_string()}));
        println!(
            "{}",
            serde_json::json!({"case":r["case"],"native":native,"immediate_base_knots":base_knots})
        );
    }
}
