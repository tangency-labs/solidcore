//! Geometry as data (docs/layers.md, "Record"): the portable B-rep record with its serde form
//! and structural validation ([`brep`]), its history and metadata records, and the native body
//! file ([`native`]).
//!
//! The record computes no geometric quantity: no square roots, trigonometry, dot or cross
//! products, normalising or solves (clippy `disallowed-methods`, this crate's `clippy.toml`). It
//! depends on no geometry or maths crate. Evaluating, measuring and checking the geometry a
//! record holds is the kernel's.
#![forbid(unsafe_code)]
#![deny(clippy::disallowed_methods)]

// The record's types are SolidCraft's wire format, carried over as it is: each type says what it
// is, its formula and its SMLib and STEP counterparts, and its fields are the names in that
// formula (SolidCraft's docs/native-brep.md), so fields and variants are not documented one by
// one.
#[allow(missing_docs)]
pub mod brep;
pub mod native;
