//! Runner-local `.sms` surface decoding, also included by non-C++ regression tests.
//!
//! Deviations: reject invalid header counts and overflowing dimensions, check the
//! available token count before decoding, and grow storage only as tokens are read.
//! Valid files retain the original float parsing, point order and knot values.

use serde::{Deserialize, Serialize};
use solidcore::nurbs::math::NOW;

use crate::HarnessError;

/// An inline surface: storage plus optional indices in use (default: all).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub(super) struct InlineSrf {
    /// Control net, with the v index fastest.
    pub(super) net: Vec<Vec<[f64; 4]>>,
    /// Degree in u.
    pub(super) p: usize,
    /// Degree in v.
    pub(super) q: usize,
    /// Knots in u.
    pub(super) u: Vec<f64>,
    /// Knots in v.
    pub(super) v: Vec<f64>,
    /// Optional indices in use.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(super) used: Option<[usize; 4]>,
}

/// Reads an upstream `.sms` text file, preserving the path in decoding errors.
pub(super) fn read_sms(path: &str) -> Result<InlineSrf, HarnessError> {
    let text =
        std::fs::read_to_string(path).map_err(|e| HarnessError::new(format!("{path}: {e}")))?;
    parse_sms(&text, path)
}

/// Decodes `n m p q`, a rational flag, points, then u and v knot vectors.
pub(super) fn parse_sms(text: &str, path: &str) -> Result<InlineSrf, HarnessError> {
    let error = |message: &str| HarnessError::new(format!("{path}: {message}"));
    let mut it = text.split_whitespace();
    let mut num = || -> Result<f64, HarnessError> {
        it.next()
            .ok_or_else(|| error("truncated"))?
            .parse::<f64>()
            .map_err(|e| HarnessError::new(format!("{path}: {e}")))
    };
    let int = |value: f64, name: &str| -> Result<usize, HarnessError> {
        // The exclusive power-of-two bound avoids rounding usize::MAX up when
        // converting it to f64 (notably on 64-bit hosts).
        let upper = 2.0_f64.powi(usize::BITS as i32);
        if !value.is_finite() || value.fract() != 0.0 || value < 0.0 || value >= upper {
            return Err(HarnessError::new(format!(
                "{path}: bad {name} count {value}"
            )));
        }
        Ok(value as usize)
    };
    let n = int(num()?, "n")?;
    let m = int(num()?, "m")?;
    let p = int(num()?, "p")?;
    let q = int(num()?, "q")?;
    let rat = num()? != 0.0;
    let rows = n
        .checked_add(1)
        .ok_or_else(|| error("row count overflow"))?;
    let cols = m
        .checked_add(1)
        .ok_or_else(|| error("column count overflow"))?;
    let u_count = n
        .checked_add(p)
        .and_then(|v| v.checked_add(2))
        .ok_or_else(|| error("u knot count overflow"))?;
    let v_count = m
        .checked_add(q)
        .and_then(|v| v.checked_add(2))
        .ok_or_else(|| error("v knot count overflow"))?;
    let required = rows
        .checked_mul(cols)
        .and_then(|v| v.checked_mul(if rat { 4 } else { 3 }))
        .and_then(|v| v.checked_add(u_count))
        .and_then(|v| v.checked_add(v_count))
        .ok_or_else(|| error("required token count overflow"))?;
    // Count tokens without collecting them or allocating from header dimensions.
    // Scan independently so num keeps the original float-parsing path.
    let available = text.split_whitespace().skip(5).count();
    if required > available {
        return Err(HarnessError::new(format!(
            "{path}: truncated: need {required} data tokens, have {available}"
        )));
    }
    let mut net = Vec::new();
    for _ in 0..rows {
        let mut row = Vec::new();
        for _ in 0..cols {
            let (x, y, z) = (num()?, num()?, num()?);
            let w = if rat { num()? } else { NOW };
            row.push([x, y, z, w]);
        }
        net.push(row);
    }
    let mut u = Vec::new();
    for _ in 0..u_count {
        u.push(num()?);
    }
    let mut v = Vec::new();
    for _ in 0..v_count {
        v.push(num()?);
    }
    Ok(InlineSrf {
        net,
        p,
        q,
        u,
        v,
        used: None,
    })
}
