//! Reproducible #3060 case audit using only the public goldens API.
//! Link with the goldens and serde_json rlibs of the build being audited.
//! Usage: audit ROOT [--inputs-only] or audit ROOT --replay INPUTS.jsonl.
//! --all-inputs-only freezes all generated and fixture All projected cases.
//! Outputs one JSON record per case; never applies known-failure exclusions.
use goldens::{Case, CaseContext, CaseSource, Output, Rng, UnitDiff};
use serde_json::{json, Value};
use std::io::{self, Write};
use std::path::Path;

fn selected(root: &Path, all: bool) -> Result<Vec<Case>, Box<dyn std::error::Error>> {
    let runner = goldens::units::curve_curve_p9::Runner;
    let mut cases = Vec::new();
    for i in 0..3000 {
        let seed = CaseContext::seed_for(runner.id(), i);
        let input = runner
            .generate(&mut Rng::new(seed))
            .ok_or("generator ended early")?;
        if input.get("op").and_then(Value::as_str) == Some("gcs")
            && input.get("sop").and_then(Value::as_u64) == Some(15)
            && input.get("requested").and_then(Value::as_u64) == Some(u64::from(all))
        {
            cases.push(Case {
                unit: runner.id().into(),
                name: format!("gen{i:04}"),
                source: CaseSource::Generated { seed },
                inputs: input,
            });
        }
    }
    if !all && cases.len() != 31 {
        return Err(format!("expected 31 Single cases, got {}", cases.len()).into());
    }
    if all {
        for case in goldens::fixtures::load_fixture_dir(root, runner.id())? {
            let input = &case.inputs;
            if input.get("op").and_then(Value::as_str) == Some("gcs")
                && input.get("sop").and_then(Value::as_u64) == Some(15)
                && input.get("requested").and_then(Value::as_u64) == Some(1)
            {
                cases.push(case);
            }
        }
        return Ok(cases);
    }
    let names = [
        "analytic_face_cyl_x_6_15",
        "analytic_face_cyl_x_7_15",
        "analytic_face_halfcyl_x_6_15",
        "analytic_face_halfcyl_x_7_15",
    ];
    let fixtures = goldens::fixtures::load_fixture_dir(root, "intersect.topology_solver_p2")?;
    for name in names {
        cases.push(
            fixtures
                .iter()
                .find(|c| c.name == name)
                .ok_or("missing cylinder case")?
                .clone(),
        );
    }
    Ok(cases)
}

fn solutions(case: &Case, output: &Output) -> Result<Value, &'static str> {
    if output.error.is_some()
        || output.rejected.is_some()
        || output.deferred.is_some()
        || output
            .ints
            .get("status")
            .is_some_and(|s| s.first() != Some(&0))
        || output
            .values
            .get("status")
            .is_some_and(|s| s.first() != Some(&0.0))
    {
        return Ok(Value::Null);
    }
    let topology = case.unit == "intersect.topology_solver_p2";
    let key = if topology { "static" } else { "out" };
    let Some(flat) = output.values.get(key) else {
        return Ok(Value::Null);
    };
    let mut cursor = if topology {
        2 + case
            .inputs
            .get("first")
            .and_then(Value::as_array)
            .ok_or("missing first")?
            .len()
            + case
                .inputs
                .get("second")
                .and_then(Value::as_array)
                .ok_or("missing second")?
                .len()
    } else {
        0
    };
    let count = *flat.get(cursor).ok_or("missing solution count")? as usize;
    cursor += 1;
    let mut result = Vec::new();
    for _ in 0..count {
        let header = flat
            .get(cursor..cursor + 6)
            .ok_or("short solution header")?;
        cursor += 6;
        let width = if topology { header[1] as usize } else { 2 };
        let params = flat.get(cursor..cursor + width).ok_or("short parameters")?;
        cursor += width;
        let mut row = json!({"type":header[0], "num_variables":header[1],
            "num_objects":header[2], "objects":&header[3..5],
            "value":header[5], "parameters":params});
        if header[0] == 2.0 {
            let end = flat
                .get(cursor..cursor + width + 1)
                .ok_or("short range end")?;
            cursor += width + 1;
            row["range_end"] = json!({"value":end[0], "parameters":&end[1..]});
        }
        if topology {
            row["nodes"] = json!(flat.get(cursor..cursor + 5).ok_or("short nodes")?);
            cursor += 5;
        }
        result.push(row);
    }
    if cursor != flat.len() {
        return Err("unparsed output tail");
    }
    Ok(json!(result))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let root = Path::new(
        args.get(1)
            .ok_or("usage: audit ROOT [--inputs-only|--replay FILE]")?,
    );
    let all = args.get(2).map(String::as_str) == Some("--all-inputs-only");
    let inputs_only = all || args.get(2).map(String::as_str) == Some("--inputs-only");
    let cases = if args.get(2).map(String::as_str) == Some("--replay") {
        let text = std::fs::read_to_string(args.get(3).ok_or("missing replay file")?)?;
        text.lines()
            .filter(|line| !line.trim().is_empty())
            .map(|line| {
                let row: Value = serde_json::from_str(line)?;
                Ok(Case {
                    unit: row["unit"].as_str().ok_or("missing unit")?.into(),
                    name: row["case"].as_str().ok_or("missing name")?.into(),
                    source: CaseSource::Inline,
                    inputs: row["input"].clone(),
                })
            })
            .collect::<Result<Vec<_>, Box<dyn std::error::Error>>>()?
    } else {
        selected(root, all)?
    };
    let mut stdout = io::BufWriter::new(io::stdout().lock());
    for case in cases {
        let mut row = json!({"unit":case.unit, "case":case.name, "input":case.inputs});
        if !inputs_only {
            let runner: &dyn UnitDiff = match case.unit.as_str() {
                "curve.curve_p9" => &goldens::units::curve_curve_p9::Runner,
                "intersect.topology_solver_p2" => {
                    &goldens::units::intersect_topology_solver_p2::Runner
                }
                _ => return Err("unsupported replay unit".into()),
            };
            let start = std::time::Instant::now();
            match runner.run_rust(&case) {
                Ok(output) => {
                    row["output"] = json!({"values":output.values, "ints":output.ints,
                        "error":output.error, "rejected":output.rejected,
                        "deferred":output.deferred.as_ref().map(|d| json!({"unit":d.unit,"symbol":d.symbol}))});
                    match solutions(&case, &output) {
                        Ok(decoded) => row["solutions"] = decoded,
                        Err(reason) => row["decode_error"] = json!(reason),
                    }
                }
                Err(error) => row["harness_error"] = json!(error.to_string()),
            }
            row["elapsed_seconds"] = json!(start.elapsed().as_secs_f64());
        }
        writeln!(stdout, "{row}")?;
        stdout.flush()?;
    }
    Ok(())
}
