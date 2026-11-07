//! Replays both seeded elbow inputs, checks the body and writes diagnostic SMB files.
//! Link against goldens, solidcore and serde_json from the same audited build.
use goldens::{CaseContext, Rng, UnitDiff};
use solidcore::curve::curve::Curve;
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::geom::axis2_placement::Axis2Placement;
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::math::Vec3;
use solidcore::nurbs::b_spline_curve::{BSplineCurve, NurbCircleParam};
use solidcore::primitive::primitive_creation::{model::PrimitiveModel, p1::PrimitiveCreation};
use solidcore::topology::api_files::write_brep_to_file;
use solidcore::topology::brep::validation::{check_brep, CheckLevel};
use solidcore::topology::brep_data_write::FromBrepEnv;
use solidcore::topology::topology::NotifyRecords;
use solidcore::Kernel;

/// Replays the two seeded elbows, validates them and writes their diagnostic bodies.
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let out = std::path::PathBuf::from(std::env::args().nth(1).ok_or("output directory required")?);
    std::fs::create_dir_all(&out)?;
    let runner = goldens::units::primitive_primitive_creation_p2::Runner;
    for index in [92, 189] {
        let input = runner
            .generate(&mut Rng::new(CaseContext::seed_for(runner.id(), index)))
            .ok_or("missing input")?;
        assert_eq!(input["op"], 9);
        let d: Vec<f64> = input["data"]
            .as_array()
            .ok_or("data")?
            .iter()
            .map(|v| v.as_f64().ok_or("non-numeric data entry"))
            .collect::<Result<_, _>>()?;
        assert_eq!(d.len(), 27);
        assert_eq!(&d[4..6], &[1., 2.]);
        let curves: Vec<Curve> = d[6..]
            .chunks_exact(7)
            .map(|c| {
                if c[0] == 0. {
                    Curve::from(
                        BSplineCurve::create_line_segment(
                            3,
                            Vec3::new(c[1], c[2], c[3]),
                            Vec3::new(c[4], c[5], c[6]),
                        )
                        .unwrap(),
                    )
                } else {
                    Curve::from(
                        BSplineCurve::create_circle_segment(
                            3,
                            &Axis2Placement::new(Vec3::new(c[1], c[2], c[3]), Vec3::X, Vec3::Y),
                            c[4],
                            c[5],
                            c[6],
                            NurbCircleParam::Quadratic,
                        )
                        .unwrap(),
                    )
                }
            })
            .collect();
        let mut k = Kernel::default();
        let (_, mut body) = PrimitiveModel::build_new(&mut k, |_, model| {
            let id = PrimitiveCreation::create_sweep_along_planar_path(
                model,
                PrimitiveModel::MAIN,
                &curves[..1],
                &curves[1..],
                d[0],
                d[1] != 0.,
                d[2] != 0.,
                d[3] != 0.,
            )?;
            model.take_brep(id)
        })?;
        body.refresh_cached_topology();
        let check = check_brep(&mut k, &mut body, CheckLevel::Level2)?;
        let store = AttributeStore::default();
        let volume = brep_volume(
            &k,
            &mut body,
            1e-12,
            &PropertyAttributes {
                store: &store,
                key_of: &|_| ObjectKey(0),
            },
        )?;
        let exact = std::f64::consts::PI * d[10] * d[10] * (d[19] + d[24]);
        let records = NotifyRecords::new();
        let env = FromBrepEnv {
            attributes: &store,
            key_of: &|e| body.owner_identity(e),
            records: &records,
            curve_key: &|_| ObjectKey(0),
            uv_curve_key: &|_| ObjectKey(0),
            surface_key: &|_| ObjectKey(0),
        };
        std::fs::write(
            out.join(format!("gen{index:04}.smb")),
            write_brep_to_file(&k, Some(&body), &env, false)?,
        )?;
        println!(
            "{}",
            serde_json::json!({"case": format!("gen{index:04}"), "input": input, "valid": check.is_valid(), "checks": check.lines(), "volume": volume, "exact_volume": exact, "relative_volume_error": (volume-exact).abs()/exact})
        );
        assert!(check.is_valid(), "{:?}", check.lines());
        assert!((volume - exact).abs() < exact * 1e-4);
    }
    Ok(())
}
