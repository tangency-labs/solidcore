//! A kernel body and the operations every family shares: primitives, copy and transform, counts
//! and volume, and the result contract ([`finish`]).
//!
//! Operations take a body and change it, as the kernel does; an application that keeps its
//! input calls [`Body::copy`] first.

use std::sync::{Mutex, MutexGuard, PoisonError};

use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::mass::integrator::properties::{PropertyAttributes, brep_volume};
use solidcore::primitive::api_primitives as primitives;
use solidcore::topology::brep::mirror::transform_matrix;
use solidcore::topology::brep::p2::is_manifold_solid;
use solidcore::topology::brep::p3::brep_copy;
use solidcore::{Brep, Vec3};

use crate::error::status;
use crate::session::{run, run_on, with_body};
use crate::{Error, Result};

/// A kernel body, with its own lock: the kernel's operations take `&mut Brep`, even to measure.
#[derive(Debug)]
pub struct Body(Mutex<Brep>);

impl Body {
    /// A body the kernel made.
    #[must_use]
    pub fn new(b: Brep) -> Body {
        Body(Mutex::new(b))
    }

    /// The kernel's body, locked. A body an earlier panic left unfinished is still handed out;
    /// [`with_body`] refuses it first.
    pub fn lock(&self) -> MutexGuard<'_, Brep> {
        self.0.lock().unwrap_or_else(PoisonError::into_inner)
    }

    /// Whether a panic interrupted an edit of this body.
    #[must_use]
    pub fn is_poisoned(&self) -> bool {
        self.0.is_poisoned()
    }

    /// The kernel's body, taken out.
    #[must_use]
    pub fn into_inner(self) -> Brep {
        self.0.into_inner().unwrap_or_else(PoisonError::into_inner)
    }

    /// A deep copy (`SmApiBrepCopy`): the same topology in the same order, with the attributes
    /// its entities carry.
    ///
    /// # Errors
    /// Those of the kernel's copy; [`Error::Internal`] for a body a panic left unfinished.
    pub fn copy(&self) -> Result<Body> {
        if self.is_poisoned() {
            return Err(Error::Internal("SmApiBrepCopy: the body was left unfinished by an earlier panic".into()));
        }
        run_on("SmApiBrepCopy", &[("brep", self)], |k| brep_copy(k, &self.lock())).map(Body::new)
    }

    /// Moves this body by the affine `matrix` (`SmBrep::Transform(matrix, TRUE)`): rows 0-2 the
    /// images of x, y and z, row 3 the translation (`p' = p·A + t`), the fourth column unread.
    /// Reflections mirror it, a non-uniform scale or a shear turns its geometry into NURBS,
    /// tolerances scale with it.
    ///
    /// # Errors
    /// [`Error::Argument`] for an entry that is not finite; the kernel's status.
    pub fn transform(&self, matrix: &[[f64; 4]; 4]) -> Result<()> {
        if !matrix.iter().flatten().all(|x| x.is_finite()) {
            return Err(Error::Argument);
        }
        with_body("SmBrep::Transform", self, |k, b| Ok(transform_matrix(k, b, matrix, true).map_err(status)))?
    }

    /// Face, edge and vertex counts.
    ///
    /// # Errors
    /// Those of [`with_body`].
    pub fn counts(&self) -> Result<[usize; 3]> {
        with_body("counts", self, |_, b| Ok([b.face_list().len(), b.edge_list().len(), b.vertex_list().len()]))
    }

    /// The enclosed volume (`SmApiBrepComputeVolume`) at `relative_accuracy`, which must be
    /// finite and is clamped to [1e-8, 1e-1].
    ///
    /// # Errors
    /// [`Error::Argument`] for an accuracy that is not finite; the kernel's error.
    pub fn volume(&self, relative_accuracy: f64) -> Result<f64> {
        let relative_accuracy = accuracy(relative_accuracy)?;
        with_body("brep_volume", self, |k, b| {
            let store = AttributeStore::default();
            let attrs = PropertyAttributes { store: &store, key_of: &|_| ObjectKey(0) };
            brep_volume(k, b, relative_accuracy, &attrs)
        })
    }
}

/// The accuracy the kernel's integrators take: finite, clamped to [1e-8, 1e-1].
///
/// # Errors
/// [`Error::Argument`] when it is not finite.
// `f64::clamp` is on the kernel's panic-site list (it panics on a reversed range); these bounds
// are constants in order.
#[allow(clippy::manual_clamp)]
pub fn accuracy(relative_accuracy: f64) -> Result<f64> {
    if !relative_accuracy.is_finite() {
        return Err(Error::Argument);
    }
    Ok(relative_accuracy.max(1e-8).min(1e-1))
}

/// What an operation's result must be ([`finish`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Expect {
    /// Any body with faces.
    Faces,
    /// A closed manifold solid: what a solid input must give.
    Solid,
}

/// The result contract: a result has faces ([`Error::Empty`] if not), and a solid input gives a
/// closed manifold solid ([`Error::NotSolid`] if not).
///
/// # Errors
/// The kernel's error from the manifold test.
pub fn finish(b: Brep, expect: Expect) -> solidcore::Result<Result<Brep>> {
    if b.face_list().is_empty() {
        return Ok(Err(Error::Empty));
    }
    if expect == Expect::Solid && !is_manifold_solid(&b)? {
        return Ok(Err(Error::NotSolid));
    }
    Ok(Ok(b))
}

fn point(p: [f64; 3]) -> Result<Vec3> {
    if !p.iter().all(|x| x.is_finite()) {
        return Err(Error::Argument);
    }
    Ok(Vec3::new(p[0], p[1], p[2]))
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}

/// An axis-aligned box from its lower corner, every edge length positive (`SmApiCreateBox`).
///
/// # Errors
/// [`Error::Argument`] for a number that is not finite or a length that is not positive.
pub fn make_box(origin: [f64; 3], size: [f64; 3]) -> Result<Body> {
    let o = point(origin)?;
    if !size.iter().all(|&s| positive(s)) {
        return Err(Error::Argument);
    }
    let [x, y, z] = size;
    run("create_box", |k| primitives::create_box(k, o, x, y, z)).map(Body::new)
}

/// A cylinder about +Z from its base centre, radius and height positive (`SmApiCreateCylinder`).
///
/// # Errors
/// [`Error::Argument`] for a number out of range.
pub fn make_cylinder(origin: [f64; 3], radius: f64, height: f64) -> Result<Body> {
    let o = point(origin)?;
    if !(positive(radius) && positive(height)) {
        return Err(Error::Argument);
    }
    run("create_cylinder", |k| primitives::create_cylinder(k, o, radius, height)).map(Body::new)
}

/// A sphere, radius positive (`SmApiCreateSphere`).
///
/// # Errors
/// [`Error::Argument`] for a number out of range.
pub fn make_sphere(center: [f64; 3], radius: f64) -> Result<Body> {
    let o = point(center)?;
    if !positive(radius) {
        return Err(Error::Argument);
    }
    run("create_sphere", |k| primitives::create_sphere(k, o, radius)).map(Body::new)
}

/// A cone or frustum about +Z (`SmApiCreateCone`): base radius and height positive, top radius
/// at least 0 and not the base radius.
///
/// # Errors
/// [`Error::Argument`] for a number out of range.
pub fn make_cone(origin: [f64; 3], base_radius: f64, top_radius: f64, height: f64) -> Result<Body> {
    let o = point(origin)?;
    let ok = positive(base_radius) && positive(height) && top_radius.is_finite() && top_radius >= 0.0 && top_radius != base_radius;
    if !ok {
        return Err(Error::Argument);
    }
    run("create_cone", |k| primitives::create_cone(k, o, base_radius, top_radius, height)).map(Body::new)
}

/// A ring torus about +Z, major radius > minor radius > 0 (`SmApiCreateTorus`).
///
/// # Errors
/// [`Error::Argument`] for a number out of range.
pub fn make_torus(center: [f64; 3], major_radius: f64, minor_radius: f64) -> Result<Body> {
    let o = point(center)?;
    if !(positive(minor_radius) && minor_radius < major_radius && major_radius.is_finite()) {
        return Err(Error::Argument);
    }
    run("create_torus", |k| primitives::create_torus(k, o, major_radius, minor_radius)).map(Body::new)
}
