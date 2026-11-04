//! Higher screw derivatives by closed-form composition, without differencing the net.
//! Deviations: no upstream counterpart. Orders are bounded at four in either parameter.

use crate::error::{Result, SolidError};
use crate::math::{Point2, Vec3};
use crate::nurbs::b_spline_surface_n_lib::p2::EvaluateFlags;
use crate::nurbs::knots::Side;

use super::{Motion, ProceduralSweep};

const OP: &str = "surface::procedural_sweep::derivatives";

impl ProceduralSweep {
    pub(super) fn screw_derivatives(
        &self,
        uv: Point2,
        du: usize,
        dv: usize,
    ) -> Result<Vec<Vec<Vec3>>> {
        self.screw_derivatives_with_flags(uv, du, dv, EvaluateFlags::default())
    }

    pub(crate) fn screw_derivatives_with_flags(
        &self,
        uv: Point2,
        du: usize,
        dv: usize,
        flags: EvaluateFlags,
    ) -> Result<Vec<Vec<Vec3>>> {
        if du > 4 || dv > 4 {
            return Err(SolidError::invalid_input(
                OP,
                "screw derivatives are bounded at four in each parameter",
            ));
        }
        if flags.only_upper_half && du != dv {
            return Err(SolidError::invalid_input(
                OP,
                "only the upper half needs equal orders",
            ));
        }
        let Motion::Screw(motion) = &self.motion else {
            return Err(SolidError::invalid_input(
                OP,
                "higher screw derivatives require a screw motion",
            ));
        };
        let wt = self.definition_parameters(uv);
        let domain = self.profile.natural_interval();
        let w = wt.x.clamp(domain.min(), domain.max());
        let (dw, dt) = if self.map.swapped { (dv, du) } else { (du, dv) };
        let (w_side, t_side) = if self.map.swapped {
            (flags.v_side, flags.u_side)
        } else {
            (flags.u_side, flags.v_side)
        };
        let ending_w = (w_side == Side::Right) != (self.map.b[0] < 0.0);
        let ending_t = (t_side == Side::Right) != (self.map.b[1] < 0.0);
        let profile = if let Some(offset) = &self.offset {
            offset.derivatives(motion, w, dw)?
        } else {
            self.profile.evaluate(w, dw, !ending_w, false)?
        };
        let (theta, rate) = self.turn_at_side(wt.y, ending_t);
        let (Some(&first), Some(&last)) = (self.stations.first(), self.stations.last()) else {
            return Err(SolidError::invalid_input(OP, "empty screw stations"));
        };
        let t = wt.y.clamp(first, last);
        let index = self
            .stations
            .partition_point(|&s| if ending_t { s < t } else { s <= t })
            .saturating_sub(1)
            .min(self.stations.len().saturating_sub(2));
        let (Some(&a), Some(&b)) = (self.stations.get(index), self.stations.get(index + 1)) else {
            return Err(SolidError::invalid_input(OP, "missing screw span"));
        };
        let alpha = 2.0 * ((b - a) * 0.25).tan() / (b - a);
        let x = alpha * (t - 0.5 * (a + b));
        let den = 1.0 + x * x;
        let second = -4.0 * alpha * alpha * x / (den * den);
        let third = -4.0 * alpha * alpha * alpha * (1.0 - 3.0 * x * x) / (den * den * den);
        let fourth =
            48.0 * alpha * alpha * alpha * alpha * x * (1.0 - x * x) / (den * den * den * den);
        let slope = motion
            .law
            .as_ref()
            .map_or(0.0, |law| law.slope_over(first, last));
        let (shift, shift_rate) = motion.shift(theta, slope);
        let mut source = vec![vec![Vec3::ZERO; dt + 1]; dw + 1];
        for (i, row) in source.iter_mut().enumerate() {
            let p = profile
                .get(i)
                .copied()
                .ok_or_else(|| SolidError::invalid_input(OP, "missing profile derivative"))?;
            let arm = if i == 0 { p - motion.origin + shift } else { p };
            let r0 = motion.turn(arm, theta);
            let r1 = motion.axis.cross(r0)
                + if i == 0 {
                    motion.turn(shift_rate, theta)
                } else {
                    Vec3::ZERO
                };
            let r2 = motion.axis.cross(motion.axis.cross(r0))
                + if i == 0 {
                    motion.axis.cross(motion.turn(shift_rate, theta)) * 2.0
                } else {
                    Vec3::ZERO
                };
            let r3 = motion.axis.cross(motion.axis.cross(motion.axis.cross(r0)))
                + if i == 0 {
                    motion
                        .axis
                        .cross(motion.axis.cross(motion.turn(shift_rate, theta)))
                        * 3.0
                } else {
                    Vec3::ZERO
                };
            let r4 = motion
                .axis
                .cross(motion.axis.cross(motion.axis.cross(motion.axis.cross(r0))))
                + if i == 0 {
                    motion.axis.cross(
                        motion
                            .axis
                            .cross(motion.axis.cross(motion.turn(shift_rate, theta))),
                    ) * 4.0
                } else {
                    Vec3::ZERO
                };
            let theta_first = r1
                + if i == 0 {
                    motion.axis * (motion.pitch / std::f64::consts::TAU)
                } else {
                    Vec3::ZERO
                };
            for (j, value) in row.iter_mut().enumerate() {
                *value = match j {
                    0 => {
                        r0 + if i == 0 {
                            motion.origin + motion.axis * motion.lift(theta)
                        } else {
                            Vec3::ZERO
                        }
                    }
                    1 => theta_first * rate,
                    2 => r2 * (rate * rate) + theta_first * second,
                    3 => {
                        r3 * (rate * rate * rate) + r2 * (3.0 * rate * second) + theta_first * third
                    }
                    _ => {
                        r4 * rate.powi(4)
                            + r3 * (6.0 * rate * rate * second)
                            + r2 * (3.0 * second * second + 4.0 * rate * third)
                            + theta_first * fourth
                    }
                };
            }
        }
        let mut result = vec![vec![Vec3::ZERO; dv + 1]; du + 1];
        for (i, row) in result.iter_mut().enumerate() {
            for (j, value) in row.iter_mut().enumerate() {
                let (w, t) = if self.map.swapped { (j, i) } else { (i, j) };
                let derivative = source
                    .get(w)
                    .and_then(|row| row.get(t))
                    .copied()
                    .ok_or_else(|| SolidError::invalid_input(OP, "missing screw derivative"))?;
                *value = derivative * (self.map.b[0].powi(w as i32) * self.map.b[1].powi(t as i32));
                if ![value.x, value.y, value.z].iter().all(|x| x.is_finite()) {
                    return Err(SolidError::invalid_input(
                        OP,
                        "screw derivatives overflowed",
                    ));
                }
            }
        }
        Ok(result)
    }
}
