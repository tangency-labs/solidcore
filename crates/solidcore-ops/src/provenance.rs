//! The provenance hook: data an application attaches to faces and edges, which the kernel's
//! operations carry through their edits as SMLib's `SmAObject::Notify` does. A face or edge cut
//! in two hands it to both pieces, two merged into one hand theirs to the result, a copy carries
//! it, and where a boolean finds one operand's entity lying on the other's the survivor gets it.
//!
//! What the data is and how it joins is the application's ([`Provenance`]); ops only stores it on
//! the kernel's attribute store and calls the application's rule at each event. An application's
//! construction history is built on it.

use solidcore::geom::attribute::{
    ApplicationValue, Attribute, AttributeHandlers, AttributeValue, AttributeView, Behavior, ObjectKey, Placement, ids,
};
use solidcore::{Brep, EntityRef, SolidError};

/// An application's provenance data and its rules for the kernel's edits. Each rule returns what
/// the entity it names carries afterwards: `Some` replaces what it carried (`None` there before
/// means it carried none), `None` leaves it as it was.
pub trait Provenance: Clone + PartialEq + std::fmt::Debug + Send + Sync + 'static {
    /// The attribute id this data is stored under, above the kernel's own (`ids::USER_1`
    /// and up); one id per kind of data an application attaches.
    const ID: u64;

    /// `self` was on an entity cut in two; `child` is one piece and what it carries now.
    fn split(&self, child: Option<&Self>) -> Option<Self>;

    /// `self` was on one of two entities merged into one; `result` is what the survivor carries
    /// now (the survivor may be the other entity, or this one).
    fn merge(&self, result: Option<&Self>) -> Option<Self>;

    /// `self` was on an entity found lying on another, which survives; `onto` is what that one
    /// carries now.
    fn coincident(&self, onto: Option<&Self>) -> Option<Self>;
}

/// The lowest [`Provenance::ID`] an application may take.
pub const FIRST_ID: u64 = ids::USER_1;

/// The data as the kernel's application attribute, its handlers the application's rules.
#[derive(Clone, Debug, PartialEq)]
struct Carried<P: Provenance>(P);

fn on<P: Provenance>(view: &AttributeView<'_>, x: ObjectKey) -> Option<(solidcore::geom::attribute::AttributeHandle, P)> {
    let (h, a) = view.find_attribute(x, P::ID)?;
    a.application_value::<Carried<P>>().map(|c| (h, c.0.clone()))
}

/// What `x` carries afterwards, as placements: `new` replacing what it carries.
fn place<P: Provenance>(view: &AttributeView<'_>, x: ObjectKey, new: Option<P>) -> Vec<Placement> {
    let Some(new) = new else { return Vec::new() };
    if let Some((handle, _)) = view.find_attribute(x, P::ID) {
        return vec![Placement::Edit { handle, value: AttributeValue::Application(ApplicationValue::new(Carried(new))) }];
    }
    vec![Placement::Add { on: x, attribute: Attribute::application(P::ID, Behavior::Copy, Carried(new)) }]
}

impl<P: Provenance> AttributeHandlers for Carried<P> {
    fn split(&self, view: &AttributeView<'_>, split: ObjectKey, a: ObjectKey, b: ObjectKey) -> solidcore::Result<Vec<Placement>> {
        let mut out = Vec::new();
        for x in [a, b] {
            if x != split {
                let now = on::<P>(view, x).map(|(_, p)| p);
                out.extend(place(view, x, self.0.split(now.as_ref())));
            }
        }
        Ok(out)
    }

    fn merge(&self, view: &AttributeView<'_>, _a: ObjectKey, _b: ObjectKey, result: ObjectKey) -> solidcore::Result<Vec<Placement>> {
        let now = on::<P>(view, result).map(|(_, p)| p);
        Ok(place(view, result, self.0.merge(now.as_ref())))
    }

    fn coincident(&self, view: &AttributeView<'_>, _from: ObjectKey, to: ObjectKey) -> solidcore::Result<Vec<Placement>> {
        let now = on::<P>(view, to).map(|(_, p)| p);
        Ok(place(view, to, self.0.coincident(now.as_ref())))
    }
}

/// Attaches `data` to `e` (a face or an edge), replacing what of its kind `e` carried.
///
/// # Errors
/// The kernel's, for an entity not in `b`.
pub fn attach<P: Provenance>(b: &mut Brep, e: EntityRef, data: P) -> Result<(), SolidError> {
    let old = b.entity_attributes(e).into_iter().find(|(_, a)| a.attribute_id() == P::ID).map(|(h, _)| h);
    if let Some(h) = old {
        b.remove_attribute(e, h)?;
    }
    b.add_attribute(e, Attribute::application(P::ID, Behavior::Copy, Carried(data)))?;
    Ok(())
}

/// What of kind `P` `e` carries.
#[must_use]
pub fn read<P: Provenance>(b: &Brep, e: EntityRef) -> Option<P> {
    let (_, a) = b.entity_attributes(e).into_iter().find(|(_, a)| a.attribute_id() == P::ID)?;
    a.application_value::<Carried<P>>().map(|c| c.0.clone())
}
