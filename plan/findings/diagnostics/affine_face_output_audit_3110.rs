//! Freeze/replay the complete p2 cohort without exclusions or golden blessing.
//! Link against the exact build's existing goldens and serde_json rlibs.
//! Usage: audit ROOT freeze | audit ROOT replay INPUTS.jsonl
use goldens::{Case, CaseContext, CaseSource, UnitDiff};
use serde_json::{json, Value};
use std::{
    io::{self, Write},
    path::PathBuf,
};
/// Select the native affine pair fixtures without applying golden exclusions.
fn focal(input: &Value) -> bool {
    input["op"] == "static"
        && input["operation"] == 15
        && input["request"] == 0
        && matches!(input["entry"].as_u64(), Some(6 | 7))
}
/// Decode focal Single answers while retaining explicit non-success states.
fn decode(case: &Case, output: &goldens::Output) -> Result<Value, &'static str> {
    if !focal(&case.inputs) {
        return Ok(Value::Null);
    }
    if output.error.is_some()
        || output.rejected.is_some()
        || output.deferred.is_some()
        || output
            .ints
            .get("status")
            .is_some_and(|v| v.as_slice() != [0])
        || output
            .values
            .get("status")
            .is_some_and(|v| v.as_slice() != [0.0])
    {
        return Ok(Value::Null);
    }
    let flat = output.values.get("static").ok_or("missing static output")?;
    let mut pos = 2
        + case.inputs["first"]
            .as_array()
            .ok_or("missing first")?
            .len()
        + case.inputs["second"]
            .as_array()
            .ok_or("missing second")?
            .len();
    let count = *flat.get(pos).ok_or("missing count")? as usize;
    pos += 1;
    let mut answers = Vec::new();
    for _ in 0..count {
        let h = flat.get(pos..pos + 6).ok_or("short header")?;
        pos += 6;
        if h[0] != 1.0 {
            return Err("unexpected non-SingleValue focal answer");
        }
        let n = h[1] as usize;
        let parameters = flat.get(pos..pos + n).ok_or("short parameters")?;
        pos += n;
        let nodes = flat.get(pos..pos + 5).ok_or("short nodes")?;
        pos += 5;
        answers.push(json!({"type":h[0],"num_variables":h[1],"num_objects":h[2],
            "objects":&h[3..5],"value":h[5],"parameters":parameters,"nodes":nodes}));
    }
    if pos != flat.len() {
        return Err("unparsed tail");
    }
    Ok(json!(answers))
}
/// Freeze the real runner case list or replay exactly the supplied JSONL inputs.
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    let root = PathBuf::from(args.get(1).ok_or("missing root")?);
    let mode = args.get(2).ok_or("missing freeze/replay")?;
    let runner = goldens::units::intersect_topology_solver_p2::Runner;
    let cases = if mode == "freeze" {
        runner.cases(&CaseContext { root, seeds: 16 })?
    } else if mode == "replay" {
        std::fs::read_to_string(args.get(3).ok_or("missing frozen input file")?)?
            .lines()
            .filter(|s| !s.trim().is_empty())
            .map(|s| {
                let r: Value = serde_json::from_str(s)?;
                Ok(Case {
                    unit: r["unit"].as_str().ok_or("unit")?.into(),
                    name: r["case"].as_str().ok_or("case")?.into(),
                    inputs: r["input"].clone(),
                    source: CaseSource::Inline,
                })
            })
            .collect::<Result<Vec<_>, Box<dyn std::error::Error>>>()?
    } else {
        return Err("expected freeze/replay".into());
    };
    let mut out = io::BufWriter::new(io::stdout().lock());
    for case in cases {
        let mut row = json!({"unit":case.unit,"case":case.name,"input":case.inputs,
            "focal_native_pair":focal(&case.inputs)});
        if mode == "replay" {
            match runner.run_rust(&case) {
                Ok(v) => {
                    row["output"] = json!({"values":v.values,"ints":v.ints,"error":v.error,
                        "rejected":v.rejected,"deferred":v.deferred.as_ref().map(|d|json!({"unit":d.unit,"symbol":d.symbol}))});
                    match decode(&case, &v) {
                        Ok(a) => row["solutions"] = a,
                        Err(e) => row["decode_error"] = json!(e),
                    }
                }
                Err(e) => row["harness_error"] = json!(e.to_string()),
            }
        }
        writeln!(out, "{row}")?;
        out.flush()?;
    }
    Ok(())
}
