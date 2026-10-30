//! Imprint: a tool body's intersection cut into a body's faces, and a point inside each
//! face for telling the pieces apart.
//!
//! Source: tangency-labs/smlib `capi/src/smlib_imprint.cpp` (Apache-2.0, Copyright
//! SolidCraft contributors) at fork commit 4f8f998 (replay/1.0.1).
//!
//! | Rust | C API (SMLib call) |
//! |------|-------|
//! | [`Smlib::imprint`] | `smlib_imprint` (`SmApiBooleanWithOptions`, `SM_BO_IMPRINT`) |
//! | [`Smlib::face_points`] | `smlib_face_points` (`SmFace::CalculateAnInternalPoint`) |
//!
//! Deviations: [`Smlib::face_points`] returns the points instead of filling a buffer of
//! the caller's face count, and works on a deep copy of the body (the C API calls
//! `CalculateAnInternalPoint` on the body itself, whose lazy UV trims it may add).

use solidcore::boolean::api_brep::boolean_with_options;
use solidcore::boolean::merge::MergeOptions;
use solidcore::geom::a_object::BooleanOperation;
use solidcore::topology::brep::p1 as bp1;
use solidcore::topology::face::p2::calculate_an_internal_point;
use solidcore::topology::Brep;

use crate::session::{take, Fail, Smlib};
use crate::status::{Result, Status};

impl Smlib {
    /// `smlib_imprint`: a new body, `body` with each of its faces that `tool` (a sheet or
    /// a solid) meets cut along the intersection, the pieces kept as faces of one body.
    /// Both inputs are left unchanged.
    ///
    /// # Errors
    ///
    /// `SMLIB_ERR_ARG` when `body` and `tool` are the same body; `SMLIB_ERR_EMPTY` for a
    /// result without faces; `SMLIB_ERR_NO_CUT` when the tool cuts no face; the kernel's
    /// status otherwise.
    pub fn imprint(&mut self, body: &Brep, tool: &Brep) -> Result<Brep> {
        if core::ptr::eq(body, tool) {
            return Err(Status::ERR_ARG);
        }
        self.run(
            || "smlib_imprint(body, tool)".to_string(),
            &[("body", body), ("tool", tool)],
            |api| {
                let a = api.copy(body)?;
                let t = api.copy(tool)?;
                let before = bp1::num_faces(&a);
                // Imprint the intersection on A and stop; keep the tool (dropped here).
                let options = MergeOptions::new(false, true, false, false, true);
                let r = boolean_with_options(
                    api.kernel_mut(),
                    a,
                    t,
                    BooleanOperation::Imprint,
                    options,
                )?
                .result;
                let after = bp1::num_faces(&r);
                if after == 0 {
                    return Err(Fail::Status(Status::ERR_EMPTY));
                }
                if after == before {
                    return Err(Fail::Status(Status::ERR_NO_CUT));
                }
                take(r)
            },
        )
    }

    /// `smlib_face_points`: a point inside each face of `body`, in its face order: on the
    /// face's surface and inside its trims.
    ///
    /// # Errors
    ///
    /// The kernel's status (no interior found).
    pub fn face_points(&mut self, body: &Brep) -> Result<Vec<[f64; 3]>> {
        self.run(
            || format!("smlib_face_points(body, n_faces={})", bp1::num_faces(body)),
            &[("body", body)],
            |api| {
                let mut b = api.copy(body)?;
                let faces = bp1::faces(&b, None)?;
                let mut out = Vec::new();
                for f in faces {
                    let p = calculate_an_internal_point(api.kernel(), &mut b, f)?;
                    out.push([p.x, p.y, p.z]);
                }
                Ok(out)
            },
        )
    }
}
