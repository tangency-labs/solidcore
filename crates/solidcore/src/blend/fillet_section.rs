//! The G^n cross-section of a blend (opt-in, not upstream): the section of order n
//! (contact of order n with the faces along the rolling ball's contact lines, the rails)
//! and thumbweight w between the rails.
//!
//! Not in the pinned upstream: the SolidCraft fork's `source/SMLib/src/SmFilletSection.cpp`
//! (`SmGnSectionFSG`, tangency-labs/smlib@11ea2d9, which replaces the fork's
//! `SmBlendCurveCrossSectionFSG`). Here the old blend-curve section stays the default; the
//! G^n section is being ported in slices for an opt-in `UpstreamFixes` flag (#2250).
//!
//! The section, in the plane square to the leader (the centre line) at its parameter s:
//! `V(tau)` is the point of that plane at distance `tau r_i` from face i, `F_i(tau)` its
//! foot on face i (`tau = 1` at the rails), and `T_i(u) = sum_{k<=n} c_ik (u - 1)^k` the
//! order-n Taylor polynomial of `F_i` at `tau = 1`. Then
//!
//! ```text
//! S(t) = f(t) T1(1 - p1(t)) + (1 - f(t)) T2(p2(t)),  t in [0, 1],
//! f(t) = (1-t)^(n+1) / ((1-t)^(n+1) + t^(n+1)),
//! p1(t) = t / (1 - w + w t),  p2(t) = t (1 - w) / (1 - w t).
//! ```
//!
//! Over the common denominator `W(t) = ((1-t)^(n+1) + t^(n+1)) (1 - w + w t)^n (1 - w t)^n`
//! (degree 3n + 1, positive Bernstein coefficients) S is a rational Bezier curve whose
//! control points are linear in the jets `c_ik` ([`form`]).
//!
//! | Rust                                  | Fork (`SmFilletSection.cpp`)                     |
//! |---------------------------------------|--------------------------------------------------|
//! | [`interval::Ivl`] and its operators   | `sm_Ivl`, `sm_Dn`, `sm_Up`, `sm_Whole`, `sm_Checked`, `operator+ - * /` |
//! | [`interval::Ivl::sqrt`], [`interval::Ivl::mag`], [`interval::Ivl::mid`], [`interval::Ivl::hull`] | `sm_Sqrt`, `sm_Mag`, `sm_Mid`, `sm_Hull` |
//! | [`series::GnScalar`]                  | the `T` of the templates (`double`, `sm_Ivl`), `sm_Positive`, `sm_Sqrt` |
//! | [`series::Trunc`], [`series::Ser`]    | `sm_Trunc`, `sm_Ser`                              |
//! | [`series::Ser::constant`], [`series::Ser::add`], [`series::Ser::sub`], [`series::Ser::scale`], [`series::Ser::mul`], [`series::Ser::inv_sqrt`] | `sm_SerConst`, `sm_SerAdd`, `sm_SerSub`, `sm_SerScale`, `sm_SerMul`, `sm_SerInvSqrt` |
//! | [`series::V3`], [`series::VSer`] and its methods | `sm_V3`, `sm_V3From`, `sm_VSer`, `sm_VSerConst`, `sm_VSerAdd`, `sm_VSerSub`, `sm_VSerTimes`, `sm_VSerDot`, `sm_VSerCross` |
//! | [`form::bern_mul`], [`form::bern_pow`] | `sm_BernMul`, `sm_BernPow`                       |
//! | [`form::GnSectionForm`], [`form::gn_form`] | `sm_GnSectionForm`, `sm_GnForm`              |
//! | [`form::hermite_points`], [`form::hermite_knots`] | `sm_HermitePoints`, `sm_HermiteKnots` |
//!
//! Deviations:
//! - **Bounds.** `sm_Trunc`'s orders are checked against the series' fixed size
//!   (`SM_GN_MAXS`, `SM_GN_MAXE`) when a [`series::Trunc`] is made; the fork leaves them
//!   to its callers. An empty Bernstein polynomial (which the fork never forms) multiplies
//!   to an empty one instead of reading out of bounds.
//! - **`std::min`/`std::max`.** Kept as `(b < a) ? b : a` / `(a < b) ? b : a`, not
//!   `f64::min`/`f64::max`, so a NaN operand propagates as in the fork.
//! - **Out-params.** `sm_SerInvSqrt` returns its series; `sm_HermitePoints` and
//!   `sm_HermiteKnots` return theirs.

pub mod form;
pub mod interval;
pub mod series;
