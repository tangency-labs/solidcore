//! Goldens: every case's expected output, committed as `goldens/<unit>.jsonl`, and the check
//! of a run against it.
//!
//! A golden records what solidcore produced when it was blessed. It proves a result is
//! unchanged, not that it is correct: correctness is proved by closed forms, geometric checks
//! and invariants. A change to a golden is made with `--bless` and justified in its pull
//! request with numbers and pictures.
//!
//! File format, one JSON object per line, sorted by case name:
//! * line 1: `{"unit": <id>, "format": 1}`;
//! * then one line per case: `{"case", "error", "deferred", "values", "ints", "topology"}`.
//!
//! A float array of at most [`BLOCK`] entries is stored in full. A longer one is stored as
//! `{"n": N, "blocks": [[sum, abs_sum, weighted_sum], ...]}`, one triple per [`BLOCK`] entries
//! (the last block may be shorter), with `weighted_sum = Σ (j + 1) x_j` over the block, so a
//! change in any block, or entries swapped within one, shows. Integer arrays follow the same
//! rule with exact sums. Non-finite floats are written as the strings `"NaN"`, `"inf"`,
//! `"-inf"`. Topology and the error status are stored in full.
//!
//! Comparison, per float: [`values_agree`] (the 1e-9 relative rule of
//! `docs/definitions-and-certificates.md`, with the unit's documented absolute floors). Per
//! block: each of the three sums agrees within `REL_TOL` of the block's absolute sum (times
//! `BLOCK` for the weighted sum), plus the floor times the block length. Integers, topology,
//! errors and deferrals agree exactly. A case in the golden that the run does not produce, or
//! a case the run produces that the golden lacks, fails.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use serde_json::{json, Map, Value};

use crate::case::CaseContext;
use crate::compare::{values_agree, AbsFloor, REL_TOL};
use crate::output::Output;
use crate::runner::{guarded, RunOptions, UnitDiff};
use crate::HarnessError;

/// Entries per block of a long array, and the longest array stored in full.
pub const BLOCK: usize = 16;

/// Version of the file format above.
pub const FORMAT: u64 = 1;

/// Mismatches recorded per case (all are counted).
const MAX_RECORDED: usize = 8;

/// The goldens directory: `crates/difftest/goldens`.
pub fn default_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("goldens")
}

/// The golden file of `unit` in `dir`.
pub fn path_of(dir: &Path, unit: &str) -> PathBuf {
    dir.join(format!("{unit}.jsonl"))
}

fn float(x: f64) -> Value {
    if x.is_nan() {
        json!("NaN")
    } else if x.is_infinite() {
        json!(if x > 0.0 { "inf" } else { "-inf" })
    } else {
        json!(x)
    }
}

fn float_of(v: &Value) -> Option<f64> {
    match v {
        Value::String(s) => match s.as_str() {
            "NaN" => Some(f64::NAN),
            "inf" => Some(f64::INFINITY),
            "-inf" => Some(f64::NEG_INFINITY),
            _ => None,
        },
        _ => v.as_f64(),
    }
}

fn floats(xs: &[f64]) -> Value {
    if xs.len() <= BLOCK {
        return Value::Array(xs.iter().map(|&x| float(x)).collect());
    }
    let blocks: Vec<Value> = xs
        .chunks(BLOCK)
        .map(|b| {
            let (mut sum, mut abs, mut weighted) = (0.0f64, 0.0f64, 0.0f64);
            for (j, &x) in b.iter().enumerate() {
                sum += x;
                abs += x.abs();
                weighted += (j as f64 + 1.0) * x;
            }
            json!([float(sum), float(abs), float(weighted)])
        })
        .collect();
    json!({"n": xs.len(), "blocks": blocks})
}

fn ints(xs: &[i64]) -> Value {
    if xs.len() <= BLOCK {
        return json!(xs);
    }
    let blocks: Vec<Value> = xs
        .chunks(BLOCK)
        .map(|b| {
            let (mut sum, mut weighted) = (0i64, 0i64);
            for (j, &x) in b.iter().enumerate() {
                sum = sum.wrapping_add(x);
                weighted = weighted.wrapping_add((j as i64 + 1).wrapping_mul(x));
            }
            json!([sum, weighted])
        })
        .collect();
    json!({"n": xs.len(), "blocks": blocks})
}

/// One case's golden line.
pub fn encode(case: &str, out: &Output) -> Value {
    let values: Map<String, Value> = out
        .values
        .iter()
        .map(|(k, v)| (k.clone(), floats(v)))
        .collect();
    let int_map: Map<String, Value> = out.ints.iter().map(|(k, v)| (k.clone(), ints(v))).collect();
    let topology = out.topology.as_ref().map(|t| {
        json!({
            "entities": t.entities,
            "adjacency": t.adjacency,
        })
    });
    let deferred = out
        .deferred
        .as_ref()
        .map(|d| json!({"unit": d.unit, "symbol": d.symbol}));
    json!({
        "case": case,
        "error": out.error,
        "deferred": deferred,
        "values": values,
        "ints": int_map,
        "topology": topology,
    })
}

/// A golden file: its cases by name.
pub type Golden = BTreeMap<String, Value>;

/// Reads `unit`'s golden from `dir`; a missing file is an empty golden.
pub fn read(dir: &Path, unit: &str) -> Result<Golden, HarnessError> {
    let path = path_of(dir, unit);
    let text = match std::fs::read_to_string(&path) {
        Ok(t) => t,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(Golden::new()),
        Err(e) => return Err(HarnessError(format!("{}: {e}", path.display()))),
    };
    let mut lines = text.lines();
    let header: Value = serde_json::from_str(lines.next().unwrap_or("{}"))
        .map_err(|e| HarnessError(format!("{}: header: {e}", path.display())))?;
    if header.get("unit").and_then(Value::as_str) != Some(unit)
        || header.get("format").and_then(Value::as_u64) != Some(FORMAT)
    {
        return Err(HarnessError(format!(
            "{}: header is not {{\"unit\": \"{unit}\", \"format\": {FORMAT}}}",
            path.display()
        )));
    }
    let mut golden = Golden::new();
    for (i, line) in lines.enumerate() {
        let v: Value = serde_json::from_str(line)
            .map_err(|e| HarnessError(format!("{}:{}: {e}", path.display(), i + 2)))?;
        let name = v
            .get("case")
            .and_then(Value::as_str)
            .ok_or_else(|| HarnessError(format!("{}:{}: no case name", path.display(), i + 2)))?
            .to_string();
        if golden.insert(name.clone(), v).is_some() {
            return Err(HarnessError(format!(
                "{}: case `{name}` appears twice",
                path.display()
            )));
        }
    }
    Ok(golden)
}

/// Writes `unit`'s golden to `dir`, one case per line, sorted.
pub fn write(dir: &Path, unit: &str, golden: &Golden) -> Result<PathBuf, HarnessError> {
    let path = path_of(dir, unit);
    let mut body = json!({"unit": unit, "format": FORMAT}).to_string();
    body.push('\n');
    for line in golden.values() {
        body.push_str(&line.to_string());
        body.push('\n');
    }
    std::fs::create_dir_all(dir)
        .and_then(|()| std::fs::write(&path, body))
        .map_err(|e| HarnessError(format!("{}: {e}", path.display())))?;
    Ok(path)
}

/// How one case differs from its golden.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Change {
    /// Case name.
    pub case: String,
    /// Differences counted (at most [`MAX_RECORDED`] are kept in `details`).
    pub count: usize,
    /// `field: golden -> now`, one per recorded difference.
    pub details: Vec<String>,
}

struct Differ<'a> {
    floors: &'a [AbsFloor],
    change: Change,
}

impl Differ<'_> {
    fn note(&mut self, field: &str, golden: &Value, now: &Value) {
        self.change.count += 1;
        if self.change.details.len() < MAX_RECORDED {
            self.change
                .details
                .push(format!("{field}: {golden} -> {now}"));
        }
    }

    fn floor(&self, key: &str) -> Option<f64> {
        self.floors.iter().find(|f| f.key == key).map(|f| f.abs)
    }

    fn exact(&mut self, field: &str, golden: &Value, now: &Value) {
        if golden != now {
            self.note(field, golden, now);
        }
    }

    fn float_array(&mut self, field: &str, key: &str, golden: &Value, now: &Value) {
        let floor = self.floor(key);
        match (golden, now) {
            (Value::Array(g), Value::Array(n)) if g.len() == n.len() => {
                for (i, (a, b)) in g.iter().zip(n).enumerate() {
                    let agree = match (float_of(a), float_of(b)) {
                        (Some(a), Some(b)) => values_agree(a, b, floor),
                        _ => false,
                    };
                    if !agree {
                        self.note(&format!("{field}[{i}]"), a, b);
                    }
                }
            }
            (Value::Object(g), Value::Object(n)) if g.get("n") == n.get("n") => {
                let empty = Vec::new();
                let gb = g.get("blocks").and_then(Value::as_array).unwrap_or(&empty);
                let nb = n.get("blocks").and_then(Value::as_array).unwrap_or(&empty);
                if gb.len() != nb.len() {
                    self.note(&format!("{field}.blocks"), golden, now);
                    return;
                }
                for (i, (a, b)) in gb.iter().zip(nb).enumerate() {
                    if !block_agrees(a, b, floor) {
                        self.note(&format!("{field}.block[{i}]"), a, b);
                    }
                }
            }
            _ => self.note(field, golden, now),
        }
    }

    fn keyed(&mut self, section: &str, golden: &Value, now: &Value, floats: bool) {
        let empty = Map::new();
        let g = golden.as_object().unwrap_or(&empty);
        let n = now.as_object().unwrap_or(&empty);
        for key in g.keys().chain(n.keys().filter(|k| !g.contains_key(*k))) {
            let field = format!("{section}.{key}");
            match (g.get(key), n.get(key)) {
                (Some(a), Some(b)) if floats => self.float_array(&field, key, a, b),
                (Some(a), Some(b)) => self.exact(&field, a, b),
                (a, b) => self.note(
                    &field,
                    a.unwrap_or(&Value::Null),
                    b.unwrap_or(&Value::Null),
                ),
            }
        }
    }
}

/// A block triple `[sum, abs_sum, weighted_sum]` agrees with another when each sum is within
/// `REL_TOL` of the larger absolute sum (`BLOCK` times that for the weighted sum), plus the
/// floor times `BLOCK` (times `BLOCK` again for the weighted sum).
fn block_agrees(golden: &Value, now: &Value, floor: Option<f64>) -> bool {
    let parse = |v: &Value| -> Option<[f64; 3]> {
        let a = v.as_array()?;
        Some([float_of(a.first()?)?, float_of(a.get(1)?)?, float_of(a.get(2)?)?])
    };
    if golden == now {
        return true;
    }
    let (Some([s0, a0, w0]), Some([s1, a1, w1])) = (parse(golden), parse(now)) else {
        // Integer blocks are `[sum, weighted]` and agree only exactly.
        return false;
    };
    if [s0, a0, w0, s1, a1, w1].iter().any(|x| !x.is_finite()) {
        return [s0, a0, w0]
            .iter()
            .zip([s1, a1, w1].iter())
            .all(|(x, y)| x == y || (x.is_nan() && y.is_nan()));
    }
    let scale = a0.max(a1);
    let pad = floor.unwrap_or(0.0) * BLOCK as f64;
    let n = BLOCK as f64;
    (s0 - s1).abs() <= REL_TOL * scale + pad
        && (a0 - a1).abs() <= REL_TOL * scale + pad
        && (w0 - w1).abs() <= REL_TOL * n * scale + pad * n
}

/// How `now` differs from `golden`; `None` when they agree.
pub fn diff(case: &str, golden: &Value, now: &Value, floors: &[AbsFloor]) -> Option<Change> {
    let mut d = Differ {
        floors,
        change: Change {
            case: case.to_string(),
            ..Change::default()
        },
    };
    for field in ["error", "deferred", "topology"] {
        d.exact(
            field,
            golden.get(field).unwrap_or(&Value::Null),
            now.get(field).unwrap_or(&Value::Null),
        );
    }
    let null = Value::Null;
    d.keyed(
        "values",
        golden.get("values").unwrap_or(&null),
        now.get("values").unwrap_or(&null),
        true,
    );
    d.keyed(
        "ints",
        golden.get("ints").unwrap_or(&null),
        now.get("ints").unwrap_or(&null),
        false,
    );
    (d.change.count > 0).then_some(d.change)
}

/// What [`check_unit`] found for one unit.
#[derive(Debug, Clone, Default)]
pub struct GoldenResult {
    /// Unit id.
    pub unit: String,
    /// Cases run.
    pub cases: usize,
    /// Cases whose output differs from the golden.
    pub changed: Vec<Change>,
    /// Cases the run produced that the golden lacks.
    pub added: Vec<String>,
    /// Cases in the golden the run did not produce.
    pub removed: Vec<String>,
    /// Cases whose Rust side failed as a harness error or panicked: never blessed.
    pub broken: Vec<(String, String)>,
    /// Building the case list failed.
    pub harness_error: Option<String>,
    /// The run's goldens, for `--bless`.
    pub now: Golden,
}

impl GoldenResult {
    /// Whether the run matches the golden exactly in its case set and to rounding in values.
    pub fn passed(&self) -> bool {
        self.harness_error.is_none()
            && self.changed.is_empty()
            && self.added.is_empty()
            && self.removed.is_empty()
            && self.broken.is_empty()
    }
}

/// Runs the Rust side of every case of `unit` (in reverse order with `reverse`, so a result
/// that depends on what ran before shows as a change) and compares each output with `golden`.
pub fn check_unit(
    unit: &dyn UnitDiff,
    ctx: &CaseContext,
    opts: &RunOptions,
    golden: &Golden,
    reverse: bool,
) -> GoldenResult {
    let mut res = GoldenResult {
        unit: unit.id().to_string(),
        ..GoldenResult::default()
    };
    let mut cases = match unit.cases(ctx) {
        Ok(c) => c,
        Err(e) => {
            res.harness_error = Some(e.0);
            return res;
        }
    };
    if let Some(f) = &opts.filter {
        cases.retain(|c| c.name.contains(f.as_str()));
    }
    if reverse {
        cases.reverse();
    }
    let floors = unit.abs_floors();
    for case in &cases {
        res.cases += 1;
        let line = match guarded(|| unit.run_rust(case), "rust") {
            Ok(out) => encode(&case.name, &out),
            Err(e) => {
                res.broken.push((case.name.clone(), e));
                continue;
            }
        };
        match golden.get(&case.name) {
            Some(g) => {
                if let Some(c) = diff(&case.name, g, &line, floors) {
                    res.changed.push(c);
                }
            }
            None => res.added.push(case.name.clone()),
        }
        if res.now.insert(case.name.clone(), line).is_some() {
            res.harness_error = Some(format!("case `{}` appears twice", case.name));
        }
    }
    if opts.filter.is_none() {
        res.removed = golden
            .keys()
            .filter(|k| !res.now.contains_key(*k))
            .cloned()
            .collect();
    }
    res
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::output::Topology;

    fn out(values: &[f64]) -> Output {
        let mut o = Output::new();
        o.values.insert("v".into(), values.to_vec());
        o.ints.insert("n".into(), vec![values.len() as i64]);
        o.topology = Some(Topology::new().count("faces", 2));
        o
    }

    #[test]
    fn short_arrays_are_stored_in_full_and_round_trip() {
        let line = encode("c", &out(&[1.0, f64::NAN, f64::INFINITY]));
        assert_eq!(line["values"]["v"], json!([1.0, "NaN", "inf"]));
        assert!(diff("c", &line, &line, &[]).is_none());
    }

    #[test]
    fn rounding_passes_and_a_real_change_fails() {
        let g = encode("c", &out(&[1.0, 2.0]));
        let rounded = encode("c", &out(&[1.0 + 1e-13, 2.0]));
        assert!(diff("c", &g, &rounded, &[]).is_none());
        let moved = encode("c", &out(&[1.0 + 1e-6, 2.0]));
        let c = diff("c", &g, &moved, &[]).unwrap();
        assert_eq!(c.count, 1);
        assert!(c.details[0].starts_with("values.v[0]"));
    }

    #[test]
    fn long_arrays_are_blocks_that_catch_one_entry_and_a_swap() {
        let base: Vec<f64> = (0..40).map(|i| f64::from(i) * 0.5 + 1.0).collect();
        let g = encode("c", &out(&base));
        assert_eq!(g["values"]["v"]["n"], json!(40));
        assert_eq!(g["values"]["v"]["blocks"].as_array().unwrap().len(), 3);
        let mut one = base.clone();
        one[37] += 1e-6;
        let c = diff("c", &g, &encode("c", &out(&one)), &[]).unwrap();
        assert!(c.details[0].starts_with("values.v.block[2]"));
        let mut swapped = base.clone();
        swapped.swap(3, 4);
        assert!(diff("c", &g, &encode("c", &out(&swapped)), &[]).is_some());
        let mut rounded = base;
        rounded[20] *= 1.0 + 1e-14;
        assert!(diff("c", &g, &encode("c", &out(&rounded)), &[]).is_none());
    }

    #[test]
    fn integers_topology_and_errors_are_exact() {
        let g = encode("c", &out(&[1.0]));
        let mut o = out(&[1.0]);
        o.error = Some("nlib_7".into());
        assert!(diff("c", &g, &encode("c", &o), &[]).is_some());
        let mut o = out(&[1.0]);
        o.topology = Some(Topology::new().count("faces", 3));
        assert!(diff("c", &g, &encode("c", &o), &[]).is_some());
        let mut o = out(&[1.0]);
        o.ints.insert("extra".into(), vec![0]);
        assert!(diff("c", &g, &encode("c", &o), &[]).is_some());
    }

    #[test]
    fn a_floor_widens_only_its_key() {
        let floors = [AbsFloor { key: "v", abs: 1e-6, why: "test" }];
        let g = encode("c", &out(&[0.0]));
        assert!(diff("c", &g, &encode("c", &out(&[5e-7])), &floors).is_none());
        assert!(diff("c", &g, &encode("c", &out(&[5e-7])), &[]).is_some());
    }

    #[test]
    fn files_round_trip_sorted() {
        let dir = std::env::temp_dir().join(format!("goldens-test-{}", std::process::id()));
        let mut g = Golden::new();
        g.insert("b".into(), encode("b", &out(&[2.0])));
        g.insert("a".into(), encode("a", &out(&[1.0])));
        let path = write(&dir, "unit.x", &g).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        let names: Vec<&str> = text.lines().skip(1).map(|l| &l[9..10]).collect();
        assert_eq!(names, ["a", "b"]);
        assert_eq!(read(&dir, "unit.x").unwrap(), g);
        assert!(read(&dir, "other.unit").unwrap().is_empty());
        let _ = std::fs::remove_dir_all(&dir);
    }
}
