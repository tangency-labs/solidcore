"""Replay the native cone/sphere certificate without changing a Cargo target.

Example (from repository root, using the exact already-built kernel):
python plan/findings/diagnostics/projected_pair_cone_sphere_native_audit.py \
  --rlib <built-libsolidcore.rlib> --deps <matching-deps-directory> \
  --revision c7b84de47 --rustc 'rustup run 1.98.0 rustc'

The driver copies the reviewed residual helper mechanically, changing only import
paths so it links to the public kernel. It prints both source hashes and asserts
that the copy exactly matches that specified transformation. No predicate edits.
"""
import argparse
import hashlib
from pathlib import Path
import shlex
import subprocess
import tempfile

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--rlib', type=Path, required=True)
p.add_argument('--deps', type=Path, required=True)
p.add_argument('--revision', default='c7b84de47')
p.add_argument('--rustc', default='rustc')
a = p.parse_args()
repo = Path(__file__).resolve().parents[3]
source = subprocess.check_output([
    'git', '-C', str(repo), 'show',
    a.revision + ':crates/solidcore/src/intersect/topology_solver/projected_residual.rs',
])
text = source.decode()
transformed = text.replace(
    'use crate::blend::fillet_section::certificate::knot_cuts;',
    'use crate::knot_cuts;',
).replace('use crate::', 'use solidcore::').replace(
    'use solidcore::knot_cuts;', 'use crate::knot_cuts;',
)
with tempfile.TemporaryDirectory(prefix='native-cone-sphere-proof-') as folder:
    out = Path(folder)
    copied = out / 'projected_residual_audit_copy.rs'
    copied.write_bytes(transformed.encode('utf-8'))
    assert copied.read_text(encoding='utf-8') == transformed
    print('residual_source_sha256=' + hashlib.sha256(source).hexdigest(), flush=True)
    print('residual_import_adapted_sha256=' + hashlib.sha256(copied.read_bytes()).hexdigest(), flush=True)
    audit = out / 'audit.rs'
    audit.write_bytes(Path(__file__).with_suffix('.rs').read_bytes())
    binary = out / 'audit'
    subprocess.run(shlex.split(a.rustc) + [
        '--edition', '2021', str(audit), '--extern', 'solidcore=' + str(a.rlib.resolve()),
        '-L', 'dependency=' + str(a.deps.resolve()), '-o', str(binary),
    ], check=True)
    subprocess.run([str(binary)], check=True)
