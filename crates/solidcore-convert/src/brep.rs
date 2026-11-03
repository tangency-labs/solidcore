//! The record's B-rep and the kernel's B-rep data ([`solidcore::topology::brep_data::BrepData`]).

use solidcore::EntityRef;
use solidcore::geom::attribute::ObjectKey;
use solidcore::geom::object::ObjectRef;

/// The topology identity [`BrepData::from_brep`](solidcore::topology::brep_data::BrepData) and
/// `make_topology_from_data` name entities by (upstream: their addresses): the type code and the
/// entity's index.
#[must_use]
pub fn key_of(e: EntityRef) -> ObjectRef {
    let (type_code, kind, i) = match e {
        EntityRef::Brep(x) => (16011, 0u64, x.index()),
        EntityRef::Region(x) => (16012, 1, x.index()),
        EntityRef::Shell(x) => (16013, 2, x.index()),
        EntityRef::Faceuse(x) => (16022, 3, x.index()),
        EntityRef::Face(x) => (16002, 4, x.index()),
        EntityRef::Loopuse(x) => (16004, 5, x.index()),
        EntityRef::Loop(x) => (16003, 6, x.index()),
        EntityRef::Edgeuse(x) => (16005, 7, x.index()),
        EntityRef::Edge(x) => (16016, 8, x.index()),
        EntityRef::Vertexuse(x) => (16006, 9, x.index()),
        EntityRef::Vertex(x) => (16021, 10, x.index()),
    };
    ObjectRef { type_code, key: ObjectKey(((kind + 1) << 32) | u64::from(i)) }
}
