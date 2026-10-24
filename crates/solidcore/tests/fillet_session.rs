//! Actual-object lifetime controls for the production fillet session.

use solidcore::blend::fillet_corner::p3::FilletCornerKind;
use solidcore::blend::fillet_executive::p3::VertexEntry;
use solidcore::blend::fillet_executive::session::{FilletBody, FilletObject, FilletSession};
use solidcore::blend::fillet_solver::p1::FilletSurfaceGenerator;
use solidcore::geom::attribute::ids;
use solidcore::math::{Point3, PI};
use solidcore::primitive::api_primitives::create_box;
use solidcore::topology::edge;
use solidcore::{EntityRef, Kernel};

#[test]
fn actual_box_solver_stays_authoritative_through_corner_creation() {
    let mut k = Kernel::default();
    let mut target = create_box(&mut k, Point3::new(0.0, 0.0, 0.0), 10.0, 20.0, 30.0).unwrap();
    let edge = target.edges().ids().next().unwrap();
    let edgeuse = edge::blend_edgeuse(&k, &mut target, edge).unwrap().unwrap();
    let mut session = FilletSession::new(&mut k, &mut target).unwrap();
    let generator = session
        .add_surface_generator(FilletSurfaceGenerator::circular(false, 0.05, false, false))
        .unwrap();
    let solver = session
        .new_constant_radius_solver(&k, edgeuse, 0.5, 1e-4, 20.0 * PI / 180.0, 10.0 * PI / 180.0)
        .unwrap();
    session.set_solver_generator(solver, generator).unwrap();
    let actual_base = std::ptr::from_ref(session.solver(solver).unwrap().base());
    assert!(session
        .solver(solver)
        .unwrap()
        .base()
        .surface(0)
        .unwrap()
        .is_some());
    assert!(session
        .solver(solver)
        .unwrap()
        .base()
        .surface(1)
        .unwrap()
        .is_some());
    session.load_solver(&k, solver).unwrap();
    assert_eq!(
        session.executive().solvers,
        vec![FilletObject::Solver(solver)]
    );
    assert!(session
        .solver(solver)
        .unwrap()
        .base()
        .fillet_executive()
        .is_some());
    assert_eq!(session.executive().filleted_vertices.len(), 2);
    session.create_corners(&k, None).unwrap();
    assert_eq!(session.executive().corners.len(), 2);
    assert!(std::ptr::eq(
        actual_base,
        session.solver(solver).unwrap().base()
    ));
    for vertex in &session.executive().filleted_vertices {
        let Some(FilletObject::Topology(FilletBody::Target, EntityRef::Vertex(expected))) = vertex
        else {
            panic!("actual target vertex required");
        };
        let Some(VertexEntry::Corner(FilletObject::Corner(id))) =
            session.executive().vertex_corner_map.get(vertex)
        else {
            panic!("actual retained corner required");
        };
        assert_eq!(session.corner(*id).unwrap().vertex, Some(*expected));
        assert_eq!(session.corner(*id).unwrap().kind, FilletCornerKind::NByOne);
        session.target().vertex(*expected).unwrap();
    }
    let key = session.target_identity(EntityRef::Edge(edge));
    assert!(session
        .attributes()
        .find_attribute(key.key, ids::MASS_PROPERTIES)
        .is_none());
    assert_ne!(
        FilletObject::Topology(FilletBody::Target, EntityRef::Edge(edge)),
        FilletObject::Topology(FilletBody::Pseudo, EntityRef::Edge(edge))
    );
}

#[test]
fn native_api_setup_uses_owned_parameter_algorithm_and_unknown_xsect_default() {
    use solidcore::blend::fillet_executive::session::SolverObject;
    use solidcore::blend::fillet_solver::p1::FilletSurfaceGeneratorType;

    let mut k = Kernel::default();
    let mut target = create_box(&mut k, Point3::new(0.0, 0.0, 0.0), 10.0, 20.0, 30.0).unwrap();
    let edge = target.edges().ids().next().unwrap();
    let tolerance = edge::tolerance(&target, edge).unwrap();
    let mut session = FilletSession::new(&mut k, &mut target).unwrap();
    session
        .prepare_constant_radius_edges(&mut k, &[edge], 0.5, 99, 1, 1.0)
        .unwrap();
    assert!(session.executive().initialized);
    // The native SetFilletParameters reset follows the API's earlier TRUE write.
    assert!(!session.executive().options.global_merge);
    assert_eq!(session.executive().surface_generators.len(), 1);
    let [FilletObject::Solver(id)] = session.executive().solvers.as_slice() else {
        panic!("one actual solver required");
    };
    let SolverObject::ConstantRadius(solver) = session.solver(*id).unwrap() else {
        panic!("actual constant-radius constructor required");
    };
    assert_eq!(solver.fillet_radii, [0.5, 0.5]);
    assert_eq!(solver.solver.approx_tol_3d, 4.0 * tolerance);
    assert_eq!(
        solver.solver.fsg.unwrap().fillet_surface_generator_type(),
        FilletSurfaceGeneratorType::Circular
    );
    session.create_corners(&k, None).unwrap();
    assert_eq!(session.executive().corners.len(), 2);
}

#[test]
fn real_box_nx1_topology_uses_live_pseudo_objects() {
    use solidcore::blend::fillet_geom::p1::{FilletEdgeType, FilletVertexType};
    use solidcore::intersect::point_class::PointClassType;

    let mut k = Kernel::default();
    let mut target = create_box(&mut k, Point3::new(0.0, 0.0, 0.0), 10.0, 20.0, 30.0).unwrap();
    let edge = target.edges().ids().next().unwrap();
    let mut session = FilletSession::new(&mut k, &mut target).unwrap();
    session
        .prepare_constant_radius_edges(&mut k, &[edge], 0.5, 1, 1, 1.0)
        .unwrap();
    session.create_corners(&k, None).unwrap();
    let FilletObject::Corner(id) = session.executive().corners[0] else {
        panic!("actual corner required");
    };
    session.make_nx1_corner_topology(&k, id).unwrap();
    let corner = session.corner(id).unwrap();
    assert_eq!(corner.vertices.len(), 2);
    assert_eq!(corner.edges.len(), 1);
    for &vertex in &corner.vertices {
        let vertex = session.pseudo_brep().fillet_vertex(vertex).unwrap();
        assert_eq!(vertex.fillet_corner(), Some(id));
        assert_eq!(vertex.fillet_vertex_type(), FilletVertexType::RailXEdgeuse);
        assert_eq!(vertex.point_class(), PointClassType::Edgeuse);
        let Some(EntityRef::Edgeuse(original)) = vertex.point_class_object() else {
            panic!("actual classified target edge use required");
        };
        session.target().edgeuse(original).unwrap();
    }
    let rib = session.pseudo_brep().fillet_edge(corner.edges[0]).unwrap();
    assert_eq!(rib.fillet_edge_type(), FilletEdgeType::FilletXSideFace);
    session.target().face(rib.original_face().unwrap()).unwrap();
}

#[test]
fn real_box_corner_vertex_stage_retains_native_geometry_and_rail_start() {
    let mut k = Kernel::default();
    let mut target = create_box(&mut k, Point3::ZERO, 10.0, 20.0, 30.0).unwrap();
    let edge = target.edges().ids().next().unwrap();
    let mut session = FilletSession::new(&mut k, &mut target).unwrap();
    session
        .prepare_constant_radius_edges(&mut k, &[edge], 0.5, 1, 1, 1.0)
        .unwrap();
    session.create_corners(&k, None).unwrap();
    session.process_smooth_edge_sequences(&mut k).unwrap();
    session.create_corner_vertices(&mut k).unwrap();
    let mut points = Vec::new();
    for &object in &session.executive().corners {
        let FilletObject::Corner(id) = object else {
            panic!("actual corner required");
        };
        let corner = session.corner(id).unwrap();
        assert!(!corner.is_calc_as_const_rad());
        for &vertex in &corner.vertices {
            assert!(session
                .pseudo_brep()
                .fillet_vertex(vertex)
                .unwrap()
                .is_processed());
            points.push(session.pseudo_brep().vertex_point(vertex).unwrap());
        }
    }
    // Direct native box-edge component trace, in actual corner/FV order.
    assert_eq!(
        points,
        [
            Point3::new(9.5, 0.0, 30.0),
            Point3::new(10.0, 0.0, 29.5),
            Point3::new(10.0, 20.0, 29.5),
            Point3::new(9.5, 20.0, 30.0),
        ]
    );
    let rail = session.pseudo_brep().edge_list()[0];
    assert_eq!(
        session
            .pseudo_brep()
            .fillet_edge(rail)
            .unwrap()
            .fillet_edge_type(),
        solidcore::blend::fillet_geom::p1::FilletEdgeType::Rail
    );
    let first_use = session.pseudo_brep().brep().edge(rail).unwrap().edgeuses[0];
    let FilletObject::Corner(first_corner) = session.executive().corners[0] else {
        panic!("actual first corner required");
    };
    let first_vertex = session.corner(first_corner).unwrap().vertices[0];
    let vertexuse = session
        .pseudo_brep()
        .brep()
        .edgeuse(first_use)
        .unwrap()
        .vertexuse
        .unwrap();
    assert_eq!(
        session
            .pseudo_brep()
            .brep()
            .vertexuse(vertexuse)
            .unwrap()
            .vertex,
        Some(first_vertex)
    );
}

#[test]
fn real_box_post_vertex_adjustment_preserves_retained_corner_topology() {
    let mut k = Kernel::default();
    let mut target = create_box(&mut k, Point3::ZERO, 10.0, 20.0, 30.0).unwrap();
    let edge = target.edges().ids().next().unwrap();
    let mut session = FilletSession::new(&mut k, &mut target).unwrap();
    session
        .prepare_constant_radius_edges(&mut k, &[edge], 0.5, 1, 1, 1.0)
        .unwrap();
    session.create_corners(&k, None).unwrap();
    session.create_corner_vertices(&mut k).unwrap();
    let corners = session.executive().corners.clone();
    let edges = session.pseudo_brep().edge_list().to_vec();
    let mut points = Vec::new();
    for &object in &corners {
        let FilletObject::Corner(id) = object else {
            panic!("actual corner required");
        };
        for &vertex in &session.corner(id).unwrap().vertices {
            points.push((vertex, session.pseudo_brep().vertex_point(vertex).unwrap()));
        }
    }
    session.adjust_corner_topology(&mut k).unwrap();
    assert_eq!(session.executive().corners, corners);
    assert_eq!(session.pseudo_brep().edge_list(), edges);
    for (vertex, point) in points {
        assert_eq!(session.pseudo_brep().vertex_point(vertex).unwrap(), point);
        assert!(session
            .pseudo_brep()
            .fillet_vertex(vertex)
            .unwrap()
            .is_processed());
    }
}

#[test]
fn actual_box_solver_builds_surface_and_rails_in_retained_geometry() {
    let mut k = Kernel::default();
    let mut target = create_box(&mut k, Point3::ZERO, 10.0, 20.0, 30.0).unwrap();
    let edge = target.edges().ids().next().unwrap();
    let mut session = FilletSession::new(&mut k, &mut target).unwrap();
    session
        .prepare_constant_radius_edges(&mut k, &[edge], 0.5, 1, 1, 1.0)
        .unwrap();
    session.create_corners(&k, None).unwrap();
    session.process_smooth_edge_sequences(&mut k).unwrap();
    session.create_corner_vertices(&mut k).unwrap();
    session.adjust_corner_topology(&mut k).unwrap();
    let FilletObject::Solver(id) = session.executive().solvers[0] else {
        panic!("actual solver");
    };
    let actual_base = std::ptr::from_ref(session.solver(id).unwrap().base());
    let geometries = session.solver(id).unwrap().base().fillet_geoms.clone();
    assert!(geometries.iter().all(|&geometry| session
        .pseudo_brep()
        .geom(geometry)
        .unwrap()
        .fillet_surface
        .is_none()));
    session.calc_solver_geometry(&k, id).unwrap();
    let solver = session.solver(id).unwrap().base();
    assert_eq!(std::ptr::from_ref(solver), actual_base);
    assert_eq!(solver.fillet_geoms, geometries);
    assert!(solver.extend_before && solver.extend_after);
    for geometry in geometries {
        let geometry = session.pseudo_brep().geom(geometry).unwrap();
        assert_eq!(geometry.solver, Some(id));
        assert!(geometry.fillet_surface.is_some());
        let center = geometry.center_line_geometry().unwrap();
        let solidcore::curve::curve::CurveKind::Line(line) = center.kind() else {
            panic!("actual native analytic centerline must remain a line");
        };
        assert!(std::ptr::eq(
            geometry.center_line_curve().unwrap(),
            line.base()
        ));
        for rail in geometry.rails {
            let rail = rail.unwrap();
            assert!(edge::curve(session.pseudo_brep().brep(), rail)
                .unwrap()
                .is_some());
        }
    }
}

#[test]
fn smooth_prelude_requires_live_tangent_context_for_multiple_solvers() {
    let mut k = Kernel::default();
    let mut target = create_box(&mut k, Point3::ZERO, 10.0, 20.0, 30.0).unwrap();
    let edges: Vec<_> = target.edges().ids().take(2).collect();
    let mut session = FilletSession::new(&mut k, &mut target).unwrap();
    session.process_smooth_edge_sequences(&mut k).unwrap();
    session
        .prepare_constant_radius_edges(&mut k, &edges, 0.5, 1, 1, 1.0)
        .unwrap();
    assert_eq!(session.executive().solvers.len(), 2);
    let solvers = session.executive().solvers.clone();
    let error = session.process_smooth_edge_sequences(&mut k).unwrap_err();
    assert!(error.is_unported());
    assert_eq!(session.executive().solvers, solvers);
}
