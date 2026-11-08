//! #3094: a shell of a manifold solid reported success for a body that is no manifold
//! solid. The #1528 drilled plate (`fixtures/offset.offset_executive_p2/data/issue1528_plate.smb`,
//! a manifold solid) shelled inward by 2, open at its top face, with the convex edges
//! extended, the offset not merged (`merge_results` false) and an offset solid created
//! (the goldens case `issue1795_plate_api_shell_2_top_offset_solid_unmerged`) came back
//! with 16 faces, 40 edges and 29 vertices: Euler characteristic 5, no closed surface.
//!
//! A shell of a manifold solid is a manifold solid or a refusal: the ops contract
//! (`finish`, docs/validity-and-tolerances.md section 4) at the kernel's entry point.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers

#[path = "support/offset_fixture.rs"]
mod offset_fixture;

use offset_fixture::{api_shell, planar_faces_at};
use solidcore::error::offset::{OffsetFailure, Unsoundness};
use solidcore::error::Reason;
use solidcore::geom::attribute::AttributeStore;
use solidcore::offset::offset_executive::p2::{OffsetOptions, ShellOptions};
use solidcore::query::api_queries::brep_is_manifold_solid;
use solidcore::topology::api_files::read_brep_from_file;
use solidcore::topology::brep::p10::TopologyEnv;
use solidcore::topology::brep::validation::{check_brep, CheckLevel};
use solidcore::topology::topology::NotifyRecords;
use solidcore::{Brep, Kernel, SolidError};

fn plate(k: &mut Kernel) -> Brep {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let bytes =
        std::fs::read(root.join("fixtures/offset.offset_executive_p2/data/issue1528_plate.smb"))
            .unwrap();
    let mut attributes = AttributeStore::new();
    let mut records = NotifyRecords::new();
    let fresh = Brep::new();
    let key_of = |e| fresh.owner_identity(e);
    let mut env = TopologyEnv {
        all_attributes: &[],
        attributes: &mut attributes,
        key_of: &key_of,
        records: &mut records,
        next_key: 1 << 50,
    };
    let mut read = read_brep_from_file(k, Some(bytes), true, false, &mut env).unwrap();
    read.refresh_cached_topology();
    read
}

#[test]
fn a_shell_that_is_no_solid_is_refused() {
    let mut k = Kernel::default();
    let mut part = plate(&mut k);
    assert!(
        brep_is_manifold_solid(&part).unwrap(),
        "the input is a solid"
    );
    let open = planar_faces_at(&k, &part, 2, 27.5);
    assert_eq!(open.len(), 1, "the top face");
    let options = ShellOptions {
        offset: OffsetOptions {
            distance: -2.0,
            extend_convex_edges: true,
            merge_results: false,
        },
        create_offset_solid: true,
        merge_sheet_side_walls: false,
    };
    let mut body = match api_shell(&mut k, &mut part, options, &open) {
        Ok(body) => body,
        Err(error) => {
            let SolidError::Failed(failure) = &error else {
                panic!("not a refusal: {error}");
            };
            assert_eq!(
                failure.reason,
                Reason::Offset(OffsetFailure::Unsound {
                    what: Unsoundness::NotManifoldSolid
                })
            );
            return;
        }
    };
    let manifold = brep_is_manifold_solid(&body).unwrap();
    let lines = check_brep(&mut k, &mut body, CheckLevel::Level2)
        .unwrap()
        .lines();
    panic!(
        "reported as built: {} faces, {} edges, {} vertices, manifold solid {manifold}, \
         check_brep level 2 {} lines, the first {:?}",
        body.face_list().len(),
        body.edge_list().len(),
        body.vertex_list().len(),
        lines.len(),
        lines.first()
    );
}
