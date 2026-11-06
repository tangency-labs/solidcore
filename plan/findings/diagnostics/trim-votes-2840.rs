//! Replay finite project-and-trim classification without production instrumentation.
//! Link against the exact candidate solidcore rlib. No dependencies except solidcore.
//! Split returns after UV-owner binding, whereas KeepPoint classifies before that
//! binding. This replay is not an in-function execution trace; compare actual KeepPoint.
use solidcore::curve::curve::{p2::BSPLINE_CURVE_TYPE, Curve};
use solidcore::geom::attribute::ObjectKey;
use solidcore::geom::geom_utility::p1::line_closest_point;
use solidcore::geom::global_solver::{SolutionRequested, SolverOperation};
use solidcore::geom::object::{ObjectRef, Shape};
use solidcore::intersect::topology_solver::p1::{CurveRef, TopologySolver};
use solidcore::intersect::topology_solver::p2::TopologySolveOptions;
use solidcore::math::{Vec3, BIG_DOUBLE};
use solidcore::nurbs::b_spline_curve::BSplineCurve;
use solidcore::primitive::{api_curves, api_primitives, api_trim_ops};
use solidcore::query::brep::{self as qb, BoundingBoxOptions};
use solidcore::topology::brep::{p1, p10::TrimType, p2};
use solidcore::topology::brep_cache::{make_brep_cache, CacheObject, CacheSource};
use solidcore::topology::{face::p4 as face4, loopuse};
use solidcore::topology::{Brep, EdgeId};
use solidcore::{Kernel, Result, SolidError};

struct Fixture {
    name: &'static str,
    origin: Vec3,
    length: f64,
    width: f64,
    height: f64,
    direction: Vec3,
    reference: Vec3,
    cutter: Cutter,
}
#[derive(Clone, Copy)]
enum Cutter {
    Segment(Vec3, Vec3),
    Arc {
        center: Vec3,
        radius: f64,
        start_deg: f64,
        end_deg: f64,
    },
}
fn fixture(name: &str) -> Result<Fixture> {
    Ok(match name {
        // generated:seed=12099798707737135304; exact serialized difftest input.
        "gen0070" => Fixture {
            name: "gen0070",
            origin: Vec3::new(-1.2263580749873166, -1.5599815699467587, 1.7230791658768143),
            length: 12.038917547066701,
            width: 5.624041194485802,
            height: 2.2849049798121,
            direction: Vec3::new(-0.013317243772994564, 0.02179333640208525, -1.0),
            reference: Vec3::new(4.793100698546034, -1.2787795102224686, 4.007984145688914),
            cutter: Cutter::Arc {
                center: Vec3::new(5.54295238374327, -2.6847898088439193, 9.007984145688914),
                radius: 2.911016437613864,
                start_deg: 5.717811128273231,
                end_deg: 147.61560969045235,
            },
        },
        // generated:seed=12099798707737135318; exact serialized difftest input.
        "gen0084" => Fixture {
            name: "gen0084",
            origin: Vec3::new(1.918210854095832, 0.12859772202508601, -1.4894456134355734),
            length: 3.382029559475042,
            width: 11.853244820201779,
            height: 12.211865355143106,
            direction: Vec3::new(0.09905487761716145, 0.04457169745233888, -1.0),
            reference: Vec3::new(3.6092256338333533, 0.7212599630351749, 10.722419741707533),
            cutter: Cutter::Segment(
                Vec3::new(0.22719607435831102, 2.8788985533302607, 15.722419741707533),
                Vec3::new(6.991255193308396, 10.760981248163384, 15.722419741707533),
            ),
        },
        // generated:seed=12099798707737135320; exact serialized difftest input.
        "gen0086" => Fixture {
            name: "gen0086",
            origin: Vec3::new(
                -0.4582742771281505,
                -0.6952678071585972,
                -0.13803199423746282,
            ),
            length: 14.265945282705923,
            width: 10.276070401901467,
            height: 6.361937142049744,
            direction: Vec3::new(
                0.09148384575384799,
                0.08585647501698548,
                -2.2695534158749506,
            ),
            reference: Vec3::new(6.674698364224811, -0.1814642870635238, 6.223905147812281),
            cutter: Cutter::Arc {
                center: Vec3::new(6.149438126529027, -2.7504818875388906, 11.223905147812282),
                radius: 7.311995406273734,
                start_deg: 31.236694384057458,
                end_deg: 174.6285077384253,
            },
        },
        _ => {
            return Err(SolidError::invalid_input(
                "diagnostic",
                "expected gen0070, gen0084, or gen0086",
            ))
        }
    })
}
fn box_fixture(k: &mut Kernel, f: &Fixture) -> Result<Brep> {
    api_primitives::create_box(k, f.origin, f.length, f.width, f.height)
}
fn cutter(f: &Fixture) -> Result<Curve> {
    match f.cutter {
        Cutter::Segment(a, b) => {
            println!("cutter=segment start={a:?} end={b:?}");
            Ok(Curve::from(api_curves::create_line_segment(a, b)?))
        }
        Cutter::Arc {
            center,
            radius,
            start_deg,
            end_deg,
        } => {
            println!("cutter=arc center={center:?} radius={radius:.17e} start_deg={start_deg:.17e} end_deg={end_deg:.17e}");
            // Exactly the difftest runner's arc constructor, including degree inputs.
            // Do not substitute an analytic full circle or manually convert angles.
            Ok(Curve::from(api_curves::create_arc(
                center, radius, start_deg, end_deg,
            )?))
        }
    }
}

// Exact public-error equivalent of the private error::ignored helper.
fn ignored<T>(r: Result<T>) -> Result<Option<T>> {
    match r {
        Ok(v) => Ok(Some(v)),
        Err(
            e @ (SolidError::Unported { .. }
            | SolidError::Cancelled { .. }
            | SolidError::FaceTrimsUnavailable { .. }
            | SolidError::PreparationRequired { .. }),
        ) => Err(e),
        Err(_) => Ok(None),
    }
}

fn side(k: &Kernel, c: &Curve, p: Vec3, n: Vec3, tol: f64) -> Result<Option<bool>> {
    let r = p1::find_side_of_curve(k, c, p, n, tol);
    // find_side uses natural interval, ProjectedMinimize, Single, EFF_ZERO_SQRT.
    match &r {
        Ok(v) => println!(
            "side point={p:?} side_tol={tol:.17e} status=0 side={}",
            u8::from(*v)
        ),
        Err(e) => println!(
            "side point={p:?} side_tol={tol:.17e} status={:?} error={e:?}",
            e.sm_status()
        ),
    }
    ignored(r)
}

fn dump(k: &Kernel, b: &Brep, label: &str) -> Result<()> {
    println!("topology label={label} faces={:?}", p1::faces(b, None)?);
    for f in p1::faces(b, None)? {
        println!(
            "face label={label} id={f:?} tolerance={:.17e}",
            b.face(f)?.tolerance
        );
        for lu in face4::upward_loopuses(b, f, false, None)? {
            println!(
                "loop face={f:?} id={lu:?} vertices={:?}",
                loopuse::vertices(b, lu, None)?
            );
            for eu in loopuse::edgeuses(b, lu, None, None)? {
                let use_ = b.edgeuse(eu)?;
                println!(
                    "edgeuse id={eu:?} edge={:?} vertexuse={:?} orientation={:?}",
                    use_.edge, use_.vertexuse, use_.orientation
                );
            }
        }
        // Preserve traversal multiplicities: do NOT deduplicate these vertices.
        for v in face4::vertices(b, f, None, false, None)? {
            println!(
                "vertex face={f:?} id={v:?} point={:?} tolerance={:.17e}",
                b.vertex(v)?.point,
                b.vertex(v)?.tolerance
            );
        }
        for e in face4::edges(b, f, None, false, None)? {
            let e = b.edge(e)?;
            println!("edge vertices={:?} interval={:?}", e.vertices, e.interval);
        }
    }
    let _ = k;
    Ok(())
}

fn run() -> Result<()> {
    let name = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "gen0084".to_owned());
    let f = fixture(&name)?;
    let mut k = Kernel::default();
    let c = cutter(&f)?;
    let n_raw = f.direction;
    let n = n_raw.unitize()?;
    let reference = f.reference;
    println!("fixture={} origin={:?} dimensions=({:.17e},{:.17e},{:.17e}) natural={:?} raw_direction={n_raw:?} direction={n:?} reference={reference:?}", f.name, f.origin, f.length, f.width, f.height, c.natural_interval());
    println!(
        "side_query=ProjectedMinimize requested=Single objective_tol={:.17e}",
        solidcore::math::EFF_ZERO_SQRT
    );
    let original = box_fixture(&mut k, &f)?;
    // Capture identities BEFORE projection, matching api_trim_ops::project_onto.
    let proven = qb::BrepIdentity::of_brep(&original);
    let identity = proven.identity();
    let initial_tol = original.tolerance() / 4.0;
    let mut tol = initial_tol;
    let mut original_tols = Vec::new();
    for f in p1::faces(&original, None)? {
        original_tols.push(original.face(f)?.tolerance);
        tol = tol.max(original.face(f)?.tolerance);
    }
    println!("initial_tol={initial_tol:.17e} original_face_tols={original_tols:?}");
    // p10 raises tol to the maximum tolerance of projected faces as it merges.
    // For this uniform-tolerance box the projected subset cannot change that value.
    // Refuse silently assuming this if the constructor or merge changes tolerances.
    if original_tols.is_empty() || original_tols.iter().any(|v| *v != original_tols[0]) {
        return Err(SolidError::invalid_input(
            "diagnostic",
            "nonuniform input face tolerances; trace private merge tolerance",
        ));
    }
    let mut b = api_trim_ops::trim_project_parallel(
        &mut k,
        api_trim_ops::TrimObject::Brep(Box::new(original)),
        &c,
        n,
        reference,
        TrimType::Split,
    )?;
    println!("split status=0 replay_classification_tol={tol:.17e}");
    dump(&k, &b, "split")?;
    if p1::faces(&b, None)?
        .iter()
        .any(|f| b.face(*f).map_or(true, |v| v.tolerance != original_tols[0]))
    {
        return Err(SolidError::invalid_input(
            "diagnostic",
            "merge changed face tolerances; private classification tolerance requires trace",
        ));
    }
    b.set_editing_enabled(false);
    // A separate completed Split has restored doing_boolean. Re-establish the
    // classification context of project_and_trim_body for this replay.
    k.context.set_doing_boolean(true);
    let ref_side = side(&k, &c, reference, n, tol)?
        .ok_or_else(|| SolidError::invalid_input("diagnostic", "reference side failed"))?;
    let mut edges: Vec<EdgeId> = Vec::new(); // deliberately carried across faces, as p10 does
    for f in p1::faces(&b, None)? {
        println!("classify face={f:?}");
        let vertices = face4::vertices(&b, f, None, false, None)?;
        let (mut left, mut right) = (0usize, 0usize);
        for &v in &vertices {
            let rec = b.vertex(v)?;
            println!("vote_source=vertex id={v:?}");
            if let Some(s) = side(&k, &c, rec.point, n, 5.0 * (rec.tolerance + tol))? {
                if s {
                    right += 1;
                } else {
                    left += 1;
                }
            }
        }
        if right == 0 && left == 0 {
            edges = face4::edges(&b, f, None, false, None)?;
            for &e in &edges {
                let rec = b.edge(e)?;
                let curve = rec
                    .curve
                    .as_ref()
                    .ok_or_else(|| SolidError::invalid_input("diagnostic", "missing curve"))?;
                let param = rec.interval.evaluate(0.5);
                println!("vote_source=midpoint edge={e:?} param={param:.17e}");
                if let Some(s) = side(
                    &k,
                    &c,
                    curve.evaluate_point(&k, param)?,
                    n,
                    5.0 * (rec.tolerance + tol),
                )? {
                    if s {
                        right += 1;
                    } else {
                        left += 1;
                    }
                }
            }
        }
        println!("votes face={f:?} left={left} right={right}");
        let mut face_side = right >= left;
        if right == left {
            println!("tie carried_edges={edges:?}");
            let ek: Vec<_> = edges.iter().map(|&e| (identity.brep_id, e)).collect();
            let vk: Vec<_> = vertices.iter().map(|&v| (identity.brep_id, v)).collect();
            let shape = Shape::new(
                &mut k,
                &(identity.brep_id, &b),
                &[],
                &[(identity.brep_id, f)],
                &ek,
                &vk,
            )?;
            let sphere =
                qb::calculate_bounding_box(&k, &b, BoundingBoxOptions::default())?.sphere_bound();
            let center = reference + n * line_closest_point(reference, n, sphere.center);
            let start = center - n * (2.0 * sphere.radius);
            let end = center + n * (2.0 * sphere.radius);
            let ray = Curve::from(BSplineCurve::create_line_segment(3, start, end)?);
            println!("tie_query=Minimize requested=Single tol={:.17e} start={start:?} end={end:?} interval={:?}", b.tolerance(), ray.natural_interval());
            let cache = make_brep_cache(&k, CacheSource::Brep(&b), identity.keys)?;
            let brep_tol = b.tolerance();
            let objects = qb::LiveBrepObjects::new(&mut b, &identity, Some(&cache))?;
            let vectors = [Vec3::ZERO, Vec3::ZERO];
            let r = TopologySolver::shape_curve_solve(
                &k,
                &objects,
                &shape,
                CurveRef {
                    curve: &ray,
                    handle: CacheObject::Other(ObjectRef {
                        type_code: BSPLINE_CURVE_TYPE,
                        key: ObjectKey(0),
                    }),
                },
                &ray.natural_interval(),
                TopologySolveOptions {
                    operation: SolverOperation::Minimize,
                    requested: SolutionRequested::Single,
                    tolerance_3d: brep_tol,
                    best_answer_so_far_sq: BIG_DOUBLE,
                    opt_vectors: &vectors,
                },
            );
            match &r {
                Ok(s) => println!("tie_status=0 solutions={:?}", s.as_slice()),
                Err(e) => println!("tie_status={:?} error={e:?}", e.sm_status()),
            }
            let solutions = r?;
            let Some(first) = solutions.as_slice().first() else {
                println!("decision face={f:?} keep=true reason=empty_fallback");
                continue;
            };
            let point = p2::point_from_solution(&k, &b, first, 1)?;
            println!("tie_witness={point:?}");
            face_side =
                side(&k, &c, point, n, 5.0 * (b.face(f)?.tolerance + tol))?.unwrap_or(!ref_side);
        }
        println!(
            "decision face={f:?} side={} reference_side={} keep={}",
            u8::from(face_side),
            u8::from(ref_side),
            face_side == ref_side
        );
    }
    // Independently execute the original public operation, with a fresh context.
    let mut actual_k = Kernel::default();
    let mut actual = box_fixture(&mut actual_k, &f)?;
    let result = api_trim_ops::project_and_trim(&mut actual_k, &mut actual, &c, n_raw, reference);
    println!("actual_keep_point_status={result:?}");
    result?;
    dump(&actual_k, &actual, "actual_keep_point")?;
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        println!("diagnostic_status={:?} error={e:?}", e.sm_status());
        std::process::exit(1);
    }
}
