//! Record indices and the kernel's entity ids. Operations take faces, edges and vertices by their
//! index in the body's record (docs/native-brep.md, "Canonical order"); an [`Index`] says which
//! kernel entity each index names.

use solidcore::{Brep, EdgeId, FaceId, VertexId};

use crate::{Error, Result};

/// The kernel entity each record index names.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Index {
    /// `faces[i]` is record face `i`.
    pub faces: Vec<FaceId>,
    /// `edges[i]` is record edge `i`.
    pub edges: Vec<EdgeId>,
    /// `vertices[i]` is record vertex `i`.
    pub vertices: Vec<VertexId>,
}

fn at<T: Copy>(list: &[T], i: usize) -> Result<T> {
    list.get(i).copied().ok_or(Error::Argument)
}

impl Index {
    /// The kernel's own order of `b`'s faces, edges and vertices.
    #[must_use]
    pub fn of_kernel_order(b: &Brep) -> Index {
        Index { faces: b.face_list().to_vec(), edges: b.edge_list().to_vec(), vertices: b.vertex_list().to_vec() }
    }

    /// Record face `i`.
    ///
    /// # Errors
    /// [`Error::Argument`] for an index past the faces.
    pub fn face(&self, i: usize) -> Result<FaceId> {
        at(&self.faces, i)
    }

    /// Record edge `i`.
    ///
    /// # Errors
    /// [`Error::Argument`] for an index past the edges.
    pub fn edge(&self, i: usize) -> Result<EdgeId> {
        at(&self.edges, i)
    }

    /// Record vertex `i`.
    ///
    /// # Errors
    /// [`Error::Argument`] for an index past the vertices.
    pub fn vertex(&self, i: usize) -> Result<VertexId> {
        at(&self.vertices, i)
    }

    /// The record index of kernel face `f`.
    #[must_use]
    pub fn face_index(&self, f: FaceId) -> Option<usize> {
        self.faces.iter().position(|&x| x == f)
    }

    /// The record index of kernel edge `e`.
    #[must_use]
    pub fn edge_index(&self, e: EdgeId) -> Option<usize> {
        self.edges.iter().position(|&x| x == e)
    }

    /// The record index of kernel vertex `v`.
    #[must_use]
    pub fn vertex_index(&self, v: VertexId) -> Option<usize> {
        self.vertices.iter().position(|&x| x == v)
    }
}
