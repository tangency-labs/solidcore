//! Stitch: a body's faces sewn along the vertices and edges that coincide, then
//! oriented.
//!
//! Source: tangency-labs/smlib `capi/src/smlib_sheets.cpp` (`smlib_stitch`; Apache-2.0,
//! Copyright SolidCraft contributors) at fork commit b018013 (replay/1.0.1).
//!
//! | Rust | C API (SMLib call) |
//! |------|-------|
//! | [`Smlib::stitch`] | `smlib_stitch` (`SmBrep::StitchFaces`, `SmBrep::OrientTrimmedSurfaces`) |
//!
//! Deviations: the closed flag is returned with the body (the C API's `*out_closed`).

use solidcore::query::api_queries::brep_is_manifold_solid;
use solidcore::topology::brep::p2::ArenaNotifyState;
use solidcore::topology::brep::p3::orient_trimmed_surfaces;
use solidcore::topology::brep::p4::stitch_faces;
use solidcore::topology::topology_traverser::Marks;
use solidcore::topology::Brep;

use crate::session::{num, take, Smlib};
use crate::status::{Result, Status};

impl Smlib {
    /// `smlib_stitch`: a new body, `brep`'s faces sewn along the vertices and edges that
    /// coincide within `tolerance` (> 0, at most 1e3), glued once at that tolerance, then
    /// oriented: outward where every edge is glued and the faces close, else turned only to
    /// agree with their neighbours across the edges they share. The flag is `true` when
    /// the result is a closed manifold solid, `false` for a sheet.
    ///
    /// # Errors
    ///
    /// `SMLIB_ERR_ARG` for a bad tolerance; `SMLIB_ERR_EMPTY` for a result without faces;
    /// the kernel's status otherwise.
    pub fn stitch(&mut self, brep: &Brep, tolerance: f64) -> Result<(Brep, bool)> {
        let ok = tolerance.is_finite() && tolerance > 0.0 && tolerance <= 1e3;
        if !ok {
            return Err(Status::ERR_ARG);
        }
        self.run(
            || format!("smlib_stitch(brep, tolerance={})", num(tolerance)),
            &[("brep", brep)],
            |api| {
                let mut w = api.copy(brep)?;
                let lamina =
                    stitch_faces(api.kernel_mut(), &mut w, tolerance, None, false)?.lamina_edges;
                let mut state = ArenaNotifyState::new(&w);
                let mut env = state.env(None);
                let mut marks = Marks::default();
                orient_trimmed_surfaces(
                    api.kernel_mut(),
                    &mut w,
                    &mut env,
                    &mut marks,
                    true,
                    lamina != 0,
                    None,
                )?;
                let closed = lamina == 0 && brep_is_manifold_solid(&w)?;
                Ok((take(w)?, closed))
            },
        )
    }
}
