//! Interior grid lines of a face: where its Steiner points go before refinement.

use solidcore::geom::extent2d::Extent2d;
use solidcore::math::{Point2, Point3, Vec3};
use solidcore::Result;

use crate::geom::{angle, segment_distance};
use crate::options::FaceParams;

/// Most interior grid lines per UV direction of one face.
pub const MAX_GRID_LINES: usize = 2048;
/// The grid is sized for this fraction of the chord tolerance, so the triangles it
/// makes (whose diagonals sag like its cells) pass refinement's check with a margin
/// instead of failing it by rounding.
const GRID_CHORD_MARGIN: f64 = 0.7;
/// A direction whose sample intervals all need less than this many cells is straight.
const STRAIGHT: f64 = 0.05;
/// Fewest cells per direction of the sample grid that sizes a face's grid.
pub(crate) const SAMPLE_GRID: usize = 8;
/// Most cells per direction of that sample grid (a B-spline asks for one per knot span).
const MAX_SAMPLE_GRID: usize = 48;
/// Largest aspect ratio, in scaled UV, of a curved face's grid cells.
pub(crate) const MAX_CELL_ASPECT: f64 = 6.0;

/// The interior grid lines `(us, vs)` over `bbox` (strictly inside it, increasing):
/// from the chord deviation, normal turning and length between neighbouring points of
/// a sample grid, with cells of bounded aspect in scaled UV on a curved face. A flat face
/// gets none. `eval` gives the surface point and normal (if any) at a UV point.
pub(crate) fn grid_lines(
    eval: &dyn Fn(Point2) -> Result<(Point3, Option<Vec3>)>,
    bbox: &Extent2d,
    scale: (f64, f64),
    p: &FaceParams,
    samples: usize,
) -> Result<(Vec<f64>, Vec<f64>)> {
    let n = samples.clamp(SAMPLE_GRID, MAX_SAMPLE_GRID);
    let at = |i: usize, j: usize| {
        Point2::new(
            bbox.u_min() + (bbox.u_max() - bbox.u_min()) * i as f64 / n as f64,
            bbox.v_min() + (bbox.v_max() - bbox.v_min()) * j as f64 / n as f64,
        )
    };
    // sample points and normals, evaluated once
    let mut pts = Vec::with_capacity((n + 1) * (n + 1));
    for j in 0..=n {
        for i in 0..=n {
            let uv = at(i, j);
            pts.push(eval(uv)?);
        }
    }
    let sample = |i: usize, j: usize| pts.get(j * (n + 1) + i);
    let need = |a: (usize, usize), b: (usize, usize)| -> Result<f64> {
        let (Some((pa, na)), Some((pb, nb))) = (sample(a.0, a.1), sample(b.0, b.1)) else {
            return Ok(0.0);
        };
        let (ua, ub) = (at(a.0, a.1), at(b.0, b.1));
        let (pm, _) = eval(Point2::new((ua.x + ub.x) / 2.0, (ua.y + ub.y) / 2.0))?;
        let mut m = (segment_distance(pm, *pa, *pb) / (GRID_CHORD_MARGIN * p.chord)).sqrt();
        if let (Some(na), Some(nb)) = (na, nb) {
            m = m.max(angle(*na, *nb) / p.angle);
        }
        if let Some(l) = p.max_len {
            m = m.max(pa.distance(*pb) / l);
        }
        Ok(m)
    };
    // cells needed per sample interval, the worst over the rows (columns)
    let mut cu = vec![0.0_f64; n];
    let mut cv = vec![0.0_f64; n];
    for j in 0..=n {
        for (i, c) in cu.iter_mut().enumerate() {
            *c = c.max(need((i, j), (i + 1, j))?);
        }
    }
    for i in 0..=n {
        for (j, c) in cv.iter_mut().enumerate() {
            *c = c.max(need((i, j), (i, j + 1))?);
        }
    }
    let us = lines_from(&cu, bbox.u_min(), bbox.u_max());
    let vs = lines_from(&cv, bbox.v_min(), bbox.v_max());
    let lu = (bbox.u_max() - bbox.u_min()) * scale.0;
    let lv = (bbox.v_max() - bbox.v_min()) * scale.1;
    let mut us = tidy(us, bbox.u_min(), bbox.u_max());
    let mut vs = tidy(vs, bbox.v_min(), bbox.v_max());
    // cells of bounded aspect: a curved direction left coarse gets uniform lines at a
    // spacing comparable to the other's mean spacing. A straight direction (a ruling:
    // no sag, no turning) gets none, so the triangles run along the rulings, as flat as
    // their boundary chords; refinement splits any that sag after all.
    let straight = |c: &[f64]| c.iter().all(|x| *x < STRAIGHT);
    if !us.is_empty() || !vs.is_empty() {
        let hu = lu / (us.len() + 1) as f64;
        let hv = lv / (vs.len() + 1) as f64;
        let h = hu.min(hv);
        if h > 0.0 {
            if !straight(&cu) {
                us = with_uniform(
                    us,
                    bbox.u_min(),
                    bbox.u_max(),
                    (lu / (MAX_CELL_ASPECT * h)).ceil(),
                );
            }
            if !straight(&cv) {
                vs = with_uniform(
                    vs,
                    bbox.v_min(),
                    bbox.v_max(),
                    (lv / (MAX_CELL_ASPECT * h)).ceil(),
                );
            }
        }
    }
    Ok(bounded(us, vs, bbox, p.max_triangles))
}

/// Grid lines `(us, vs)` from any source made safe to insert: finite, sorted, strictly
/// inside `bbox`, without near duplicates, at most [`MAX_GRID_LINES`] per direction and
/// about `max_triangles` triangles in all (every other line dropped until they fit).
pub(crate) fn bounded(
    us: Vec<f64>,
    vs: Vec<f64>,
    bbox: &Extent2d,
    max_triangles: usize,
) -> (Vec<f64>, Vec<f64>) {
    // a huge list is thinned before it is sorted, so the work stays bounded
    let thin = |l: Vec<f64>| {
        let n = l.len().div_ceil(4 * MAX_GRID_LINES).max(1);
        if n > 1 {
            l.into_iter().step_by(n).collect()
        } else {
            l
        }
    };
    let mut us = tidy(thin(us), bbox.u_min(), bbox.u_max());
    let mut vs = tidy(thin(vs), bbox.v_min(), bbox.v_max());
    let budget = (max_triangles / 2).max(1);
    while us.len() > MAX_GRID_LINES
        || vs.len() > MAX_GRID_LINES
        || (us.len() + 1).saturating_mul(vs.len() + 1) > budget
    {
        if us.len() >= vs.len() {
            us = us.into_iter().step_by(2).collect();
        } else {
            vs = vs.into_iter().step_by(2).collect();
        }
        if us.is_empty() && vs.is_empty() {
            break;
        }
    }
    (us, vs)
}

/// Lines over `[lo, hi]` for sample intervals needing `cells[i]` cells each (a
/// direction needing less than two cells in all gets none).
fn lines_from(cells: &[f64], lo: f64, hi: f64) -> Vec<f64> {
    // the needs add up along the direction: ceil(total) cells, each holding an equal
    // share of the need, so a direction needing a third of a cell per sample interval
    // over 40 intervals gets 13 cells, and one needing 1.04 per interval is not doubled
    let total: f64 = cells.iter().sum();
    let m = total.ceil().min(MAX_GRID_LINES as f64 + 1.0);
    if m.is_nan() || m < 2.0 {
        return Vec::new();
    }
    let share = total / m;
    let n = cells.len() as f64;
    let mut out = Vec::new();
    let mut acc = 0.0;
    let mut next = share;
    for (i, c) in cells.iter().enumerate() {
        let a = lo + (hi - lo) * i as f64 / n;
        let b = lo + (hi - lo) * (i + 1) as f64 / n;
        while out.len() + 1 < m as usize && *c > 0.0 && acc + c >= next {
            out.push(a + (b - a) * (next - acc) / c);
            next += share;
        }
        acc += c;
    }
    out
}

/// `lines` sorted, strictly inside `(lo, hi)`, without near duplicates.
fn tidy(mut lines: Vec<f64>, lo: f64, hi: f64) -> Vec<f64> {
    let eps = (hi - lo) * 1e-9;
    lines.retain(|t| t.is_finite() && *t > lo + eps && *t < hi - eps);
    lines.sort_by(f64::total_cmp);
    let mut out: Vec<f64> = Vec::with_capacity(lines.len());
    for t in lines {
        if out.last().is_none_or(|l| t - l > eps) {
            out.push(t);
        }
    }
    out
}

/// `lines` plus, where two neighbours (or the ends) are farther apart than
/// `(hi - lo) / cells`, uniform lines between them.
fn with_uniform(lines: Vec<f64>, lo: f64, hi: f64, cells: f64) -> Vec<f64> {
    if cells <= 1.0 || cells.is_nan() {
        return lines;
    }
    let h = (hi - lo) / cells.min(MAX_GRID_LINES as f64);
    let mut out = Vec::with_capacity(lines.len());
    let mut prev = lo;
    for t in lines.iter().copied().chain(core::iter::once(hi)) {
        let gap = t - prev;
        if gap > h {
            let m = (gap / h).ceil() as usize;
            for kk in 1..m {
                out.push(prev + gap * kk as f64 / m as f64);
            }
        }
        if t < hi {
            out.push(t);
        }
        prev = t;
    }
    tidy(out, lo, hi)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lines_are_sorted_inside_and_bounded() {
        assert!(lines_from(&[0.2, 0.3, 0.4], 0.0, 3.0).is_empty());
        assert_eq!(lines_from(&[0.25; 40], 0.0, 1.0).len(), 9);
        let l = tidy(lines_from(&[2.0, 1.0, 3.0], 0.0, 3.0), 0.0, 3.0);
        assert_eq!(l.len(), 5);
        assert!(l.windows(2).all(|w| w[0] < w[1]) && l.iter().all(|t| *t > 0.0 && *t < 3.0));
        let u = with_uniform(vec![0.1], 0.0, 1.0, 4.0);
        assert!(u.windows(2).all(|w| w[1] - w[0] <= 0.25 + 1e-12));
        assert!(lines_from(&[1e300], 0.0, 1.0).len() <= MAX_GRID_LINES);
        // any source is bounded, sorted and inside the box
        let bbox = Extent2d::from_bounds(0.0, 0.0, 1.0, 1.0);
        let huge: Vec<f64> = (0..1_000_000).rev().map(|i| f64::from(i) / 1e6).collect();
        let (us, vs) = bounded(huge.clone(), huge, &bbox, 10_000);
        assert!(us.len() <= MAX_GRID_LINES && vs.len() <= MAX_GRID_LINES);
        assert!((us.len() + 1) * (vs.len() + 1) <= 5_000);
        assert!(us.windows(2).all(|w| w[0] < w[1]) && us.iter().all(|t| *t > 0.0 && *t < 1.0));
        let (us, _) = bounded(
            vec![f64::NAN, 2.0, -1.0, 0.5, 0.5],
            Vec::new(),
            &bbox,
            10_000,
        );
        assert_eq!(us, vec![0.5]);
    }
}
