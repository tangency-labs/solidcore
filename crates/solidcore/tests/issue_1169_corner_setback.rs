//! #1169 / WR-445: the opt-in setback correction through the public fillet session.
//! The untouched native default is covered separately by executive differential fixtures.
#![allow(clippy::unwrap_used, clippy::indexing_slicing)]

use solidcore::blend::fillet_executive::session::FilletSession;
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::mass::integrator::properties::{brep_mass_properties, PropertyAttributes};
use solidcore::math::Vec3;
use solidcore::primitive::api_primitives::create_box;
use solidcore::topology::{
    brep::{p1, p5},
    vertex,
};
use solidcore::Kernel;

const SECTIONS: [Vec3; 6] = [
    Vec3::new(52.0, 35.0, 20.0),
    Vec3::new(52.0, 40.0, 15.0),
    Vec3::new(55.0, 32.0, 20.0),
    Vec3::new(60.0, 32.0, 15.0),
    Vec3::new(55.0, 40.0, 12.0),
    Vec3::new(60.0, 35.0, 12.0),
];

#[derive(Debug)]
struct Outcome {
    counts: [usize; 3],
    vertices: Vec<Vec3>,
    bounds: [Vec3; 2],
    mass: Vec<f64>,
}

impl Outcome {
    fn bits(&self) -> Vec<u64> {
        self.vertices
            .iter()
            .chain(self.bounds.iter())
            .flat_map(|p| [p.x, p.y, p.z])
            .chain(self.mass.iter().copied())
            .map(f64::to_bits)
            .collect()
    }

    fn nearest_sections(&self) -> Vec<(usize, f64)> {
        SECTIONS
            .iter()
            .map(|wanted| {
                self.vertices
                    .iter()
                    .enumerate()
                    .map(|(i, point)| (i, point.distance(*wanted)))
                    .min_by(|a, b| a.1.total_cmp(&b.1))
                    .unwrap()
            })
            .collect()
    }
}

fn fillet(fix: bool, distance: f64) -> Outcome {
    let mut k = Kernel::default();
    k.fixes.fillet_setback_from_start_vertex = fix;
    let mut body = create_box(&mut k, Vec3::ZERO, 60.0, 40.0, 20.0).unwrap();
    let high = p1::vertices(&body, None)
        .unwrap()
        .into_iter()
        .filter(|v| vertex::point(&body, *v).unwrap() == Vec3::new(60.0, 40.0, 20.0))
        .collect::<Vec<_>>();
    assert_eq!(high.len(), 1);
    let high = high[0];
    let edges = vertex::edges(&body, high, None).unwrap();
    assert_eq!(edges.len(), 3);
    let owning = p1::edges(&body, None).unwrap();
    assert_eq!(
        edges
            .iter()
            .map(|e| owning.iter().position(|v| v == e).unwrap())
            .collect::<Vec<_>>(),
        [0, 1, 5]
    );
    let boolean_before = k.context.doing_boolean();
    {
        let mut session = FilletSession::new(&mut k, &mut body).unwrap();
        session
            .prepare_constant_radius_edges(&mut k, &edges, 5.0, 1, 1, 1.0)
            .unwrap();
        session.create_corners(&k, None).unwrap();
        let corners = session
            .fillet_corners()
            .unwrap()
            .into_iter()
            .filter(|id| session.corner(*id).unwrap().filleted_vertex() == Some(high))
            .collect::<Vec<_>>();
        assert_eq!(corners.len(), 1);
        session.set_corner_setback(corners[0], distance).unwrap();
        assert_eq!(
            session
                .corner(corners[0])
                .unwrap()
                .set_back_dist()
                .to_bits(),
            distance.to_bits()
        );
        session.do_filleting(&mut k, None).unwrap();
    }
    assert_eq!(k.context.doing_boolean(), boolean_before);
    // Match the public API lifetime: the session ends before final stitching.
    p5::stitch_and_orient(&mut k, &mut body, &p5::StitchAndOrientOptions::default()).unwrap();
    let bbox = solidcore::query::brep::calculate_tight_bounding_box(&k, &mut body).unwrap();
    let store = AttributeStore::default();
    let attributes = PropertyAttributes {
        store: &store,
        key_of: &|_| ObjectKey(0),
    };
    let mass = brep_mass_properties(&k, &mut body, 1e-8, 1.0, Vec3::ZERO, &attributes).unwrap();
    let mut values = vec![
        mass.properties.area,
        mass.properties.volume,
        mass.properties.mass,
    ];
    for p in [mass.centroid, mass.inertia, mass.products] {
        values.extend([p.x, p.y, p.z]);
    }
    let vertices = p1::vertices(&body, None)
        .unwrap()
        .into_iter()
        .map(|v| vertex::point(&body, v).unwrap())
        .collect::<Vec<_>>();
    let outcome = Outcome {
        counts: [
            vertices.len(),
            p1::edges(&body, None).unwrap().len(),
            p1::faces(&body, None).unwrap().len(),
        ],
        vertices,
        bounds: [bbox.min(), bbox.max()],
        mass: values,
    };
    assert!(outcome
        .vertices
        .iter()
        .chain(outcome.bounds.iter())
        .all(|p| [p.x, p.y, p.z].iter().all(|v| v.is_finite())));
    assert!(outcome.mass.iter().all(|v| v.is_finite()));
    outcome
}

#[test]
fn setback_correction_reaches_all_six_sections_without_changing_zero_setback() {
    let default = fillet(false, 8.0);
    let corrected = fillet(true, 8.0);
    assert_eq!(corrected.counts, [23, 36, 15]);
    for (actual, expected) in corrected
        .bounds
        .iter()
        .zip([Vec3::ZERO, Vec3::new(60.0, 40.0, 20.0)])
    {
        assert!(actual.distance(expected) < 1e-9, "{actual:?}");
    }
    let sections = corrected.nearest_sections();
    assert!(
        sections.iter().all(|(_, distance)| *distance < 1e-9),
        "{sections:?}"
    );
    let mut ids = sections.iter().map(|(id, _)| *id).collect::<Vec<_>>();
    ids.sort_unstable();
    ids.dedup();
    assert_eq!(
        ids.len(),
        6,
        "six distinct actual owning vertices are required"
    );
    // Independent convexity bounds from #1169, not a golden result from this implementation.
    assert!((47271.976..47384.893).contains(&corrected.mass[1]));
    assert!(default
        .nearest_sections()
        .iter()
        .any(|(_, distance)| *distance > 1.0));
    assert_ne!(default.bounds, corrected.bounds);
    assert_ne!(default.mass[1].to_bits(), corrected.mass[1].to_bits());
    let zero_default = fillet(false, 0.0);
    let zero_corrected = fillet(true, 0.0);
    assert_eq!(zero_default.counts, [13, 21, 10]);
    assert_eq!(zero_default.counts, zero_corrected.counts);
    assert_eq!(zero_default.bits(), zero_corrected.bits());
}
