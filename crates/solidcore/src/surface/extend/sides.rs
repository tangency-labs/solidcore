//! A surface extended past several sides of its natural box in turn, each side exactly
//! ([`Surface::extended_by_distance`]): [`Surface::extended_sides`], the one kernel entry the
//! callers that extend a surface "on every side" take.
//!
//! Not upstream: upstream's `CreateExtendedSurface` extends every side by an estimated
//! distance, fits its extension by least squares and removes knots within a tolerance; this
//! is the exact path side by side.
//!
//! Each asked side, in [`BoxSide::ALL`]'s order, is extended on the box the sides before it
//! left (the result's natural box), so a later side runs along the earlier extensions. A side
//! is left as it is, and said so, where the surface closes in its direction (a full turn: an
//! extension would lap the surface) or where the side is not proved longer than a point (a
//! pole or an apex: its certified diameter's lower bound is 0), as upstream leaves a closed
//! direction and a singular side. Every other side is extended or refused with the side's
//! own typed error.
//!
//! Deviations: none (no upstream counterpart). No recursion.

use super::{Continuity, Extended};
use crate::error::{Result, SolidError};
use crate::kernel::Kernel;
use crate::nurbs::b_spline_surface::SurfParamType;
use crate::surface::side::{Across, BoxSide};
use crate::surface::surface::{Surface, SurfaceKind};

const OP: &str = "surface::extend::Surface::extended_sides";

/// What [`Surface::extended_sides`] did at one side of the source's box.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SideExtension {
    /// Not asked.
    NotAsked,
    /// Extended: the length across it ([`Extended::across`]).
    Extended(Across),
    /// Left: the surface closes in the side's direction.
    Closed,
    /// Left: the side is not proved longer than a point.
    Collapsed,
}

/// A surface extended past several sides ([`Surface::extended_sides`]).
#[derive(Debug, Clone)]
pub struct ExtendedSides {
    /// The result.
    pub surface: Surface,
    /// What was done at each side of the source's box, in [`BoxSide::ALL`]'s order.
    pub sides: [SideExtension; 4],
    /// The largest construction bound of the sides' extensions ([`Extended::extension`]):
    /// the result over the source's box is the source (bound 0), each extension within this
    /// of its construction.
    pub extension: f64,
}

impl Surface {
    /// The surface extended `distance` past each side of its natural box in `asked`, as
    /// `how` says, exactly (module docs), a side whose direction closes or that is not proved
    /// longer than a point left as it is.
    ///
    /// # Errors
    ///
    /// Those of [`Surface::extended_by_distance`] at the first side that is refused, and of
    /// the closure and diameter queries.
    pub fn extended_sides(
        &self,
        k: &Kernel,
        asked: &[BoxSide],
        how: Continuity,
        distance: f64,
    ) -> Result<ExtendedSides> {
        let mut current: Option<Extended> = None;
        // the side of the current result each side of the source's box is on
        let mut on = BoxSide::ALL;
        let mut sides = [SideExtension::NotAsked; 4];
        let mut extension = 0.0_f64;
        for ((&source_side, done), at) in BoxSide::ALL.iter().zip(sides.iter_mut()).zip(0..) {
            if !asked.contains(&source_side) {
                continue;
            }
            let side = *on.get(at).ok_or_else(|| SolidError::sm_err(OP))?;
            let surface = current.as_ref().map_or(self, |e| &e.surface);
            let dom = surface.natural_uv_domain_with_nurb(k)?;
            let param = if side.is_u() {
                SurfParamType::U
            } else {
                SurfParamType::V
            };
            // an offset surface closes and collapses where its base does
            let mut carrier = surface;
            while let SurfaceKind::OffsetSurface(o) = carrier.kind() {
                carrier = o.base_surface();
            }
            if carrier.is_closed(k, &dom, param, None, false)?.0 {
                *done = SideExtension::Closed;
                continue;
            }
            if carrier.side_diameter(k, &dom, side)?.lower <= 0.0 {
                *done = SideExtension::Collapsed;
                continue;
            }
            let grown = surface.extended_by_distance(k, &dom, side, how, distance)?;
            for s in &mut on {
                *s = grown.side_of(*s);
            }
            extension = extension.max(grown.extension);
            *done = SideExtension::Extended(grown.across);
            current = Some(grown);
        }
        let surface = match current {
            Some(e) => e.surface,
            None => self.copy(k)?,
        };
        Ok(ExtendedSides {
            surface,
            sides,
            extension,
        })
    }
}
