"""Observe the nodal-round regression without changing kernel or test sources.

From either the base or corrected checkout, build first:
  CARGO_BUILD_JOBS=2 nice -n 10 cargo test --profile goldens -p solidcore \
      --test wr669_exact_nodal_quartic_round --no-run
  CARGO_BUILD_JOBS=2 nice -n 10 cargo build --profile goldens -p diagnostics
Then run this script by its absolute path, from that checkout:
  python3 PATH/observe.py before OUTPUT
  python3 PATH/observe.py after OUTPUT
The second invocation belongs in the corrected checkout. Output is evidence only.
"""

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

stage, directory = sys.argv[1:]
assert stage in ("before", "after")
root = Path.cwd()
output = Path(directory).resolve()
output.mkdir(parents=True, exist_ok=True)
dependencies = root / "bin/goldens/deps"


def compile_source(source, executable, libraries, test=False):
    args = ["nice", "-n", "10", "rustc", "--edition=2021", str(source), "-O",
            "-L", "dependency=" + str(dependencies), "-o", str(executable)]
    if test:
        args.append("--test")
    for name in libraries:
        library = max(dependencies.glob("lib" + name + "-*.rlib"),
                      key=lambda path: path.stat().st_mtime)
        args.extend(["--extern", name + "=" + str(library)])
    subprocess.run(args, check=True)


original = (root / "crates/solidcore/tests/wr669_exact_nodal_quartic_round.rs").read_text()
source = "\n".join(line for line in original.splitlines() if not line.startswith("//!"))
source += r'''
#[test]
fn export_nodal_body() {
    let (mut k, mut body, res, rounds) = pipes(true);
    res.unwrap();
    let attributes = solidcore::geom::attribute::AttributeStore::default();
    let records = solidcore::topology::topology::NotifyRecords::new();
    let key_of = |entity| body.owner_identity(entity);
    let env = solidcore::topology::brep_data_write::FromBrepEnv {
        attributes: &attributes, key_of: &key_of, records: &records,
        curve_key: &|_| solidcore::geom::attribute::ObjectKey(0),
        uv_curve_key: &|_| solidcore::geom::attribute::ObjectKey(0),
        surface_key: &|_| solidcore::geom::attribute::ObjectKey(0),
    };
    let path = std::path::PathBuf::from(std::env::var("SOLIDCORE_2508_EXPORT").unwrap());
    let bytes = solidcore::topology::api_files::write_brep_to_file(&k, Some(&body), &env, true).unwrap();
    std::fs::write(&path, bytes).unwrap();
    let s = face::surface(&body, rounds[0]).unwrap().unwrap();
    let (ball, rails) = off_ball(&k, s);
    let check = solidcore::topology::brep::validation::check_brep(
        &mut k, &mut body, solidcore::topology::brep::validation::CheckLevel::Level2
    ).unwrap();
    let result = serde_json::json!({
        "ball_error_mm": ball, "rail_error_mm": rails,
        "rounds": rounds.len(), "validity": format!("{check:?}")
    });
    std::fs::write(path.with_extension("measurements.json"), serde_json::to_vec_pretty(&result).unwrap()).unwrap();
    eprintln!("{result}");
}
'''
probe = output / ("observer-" + stage + ".rs")
probe.write_text(source)
executable = output / ("observer-" + stage)
compile_source(probe, executable, ["solidcore", "serde_json"], test=True)
body = output / ("nodal-" + stage + ".smb")
environment = dict(os.environ, SOLIDCORE_2508_EXPORT=str(body))
subprocess.run(["nice", "-n", "10", str(executable), "export_nodal_body", "--nocapture"],
               env=environment, check=True)

# The diagnostic renderer is unchanged except for an explicit camera frame. The
# frame is identical before and after, derived from the round's position and extent.
renderer = output / ("renderer-" + stage)
shutil.copytree(root / "crates/diagnostics/src", renderer, dirs_exist_ok=True)
render = renderer / "render.rs"
text = render.read_text()
needle = "fn fit(mesh: &Mesh, edges: &[EdgeInfo], view: View, width: usize, height: usize) -> Frame {"
camera = r'''
    if let Ok(raw) = std::env::var("SOLIDCORE_2508_VIEWPORT") {
        let p: Vec<f64> = raw.split(',').map(|x| x.parse().unwrap()).collect();
        let [x,y,_] = view.project([p[0],p[1],p[2]]);
        let scale = width as f64 / p[3];
        return Frame {scale,offset:[width as f64/2.0-x*scale,height as f64/2.0-y*scale],view};
    }
'''
assert needle in text
render.write_text(text.replace(needle, needle + camera))
executable = output / ("diagnostic-" + stage)
compile_source(renderer / "main.rs", executable,
               ["solidcore", "solidcore_convert", "solidcore_tess", "miniz_oxide", "serde_json"])
for kind, size in [("detail", "1000x800"), ("whole", "1000x800")]:
    destination = output / kind
    environment = dict(os.environ)
    if kind == "detail":
        environment["SOLIDCORE_2508_VIEWPORT"] = "2,2,0,6"
    subprocess.run(["nice", "-n", "10", str(executable), str(body), "--out", str(destination),
                    "--view", "-140,10", "--highlight", "6", "--size", size],
                   env=environment, check=True)
print(json.dumps({"stage": stage, "output": str(output)}))
