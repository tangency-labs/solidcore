"""Validate and hash the immutable #3110 baseline artifact directory."""
import argparse
import json
import pathlib
import hashlib
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('artifact_directory', type=pathlib.Path)
root = p.parse_args().artifact_directory
inputs=[json.loads(s) for s in (root/'cf110-inputs.jsonl').read_text().splitlines()]
outputs=[json.loads(s) for s in (root/'cf110-before/outputs.jsonl').read_text().splitlines()]
assert len(inputs)==len(outputs)==2575
assert [(r['case'],r['input']) for r in inputs]==[(r['case'],r['input']) for r in outputs]
assert len({r['case'] for r in inputs})==2575
focal=[r for r in outputs if r['focal_native_pair']]
assert len(focal)==14
assert not any('decode_error' in r for r in outputs)
harness_errors = sum('harness_error' in r for r in outputs)
files=['cf110-inputs.jsonl','cf110-before/outputs.jsonl','cf110-before/intersect.topology_solver_p2.debug','affine_face_output_audit_3110.rs','audit-cf110']
manifest={'baseline_revision':'cf11091a9','rustc':'1.98.0','cases':2575,'native_pair_cases':14,'generated_tree_controls':600,'fixture_cases':1975,'inputs_match_outputs':True,'decode_errors':0,'harness_errors':harness_errors,'sha256':{p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in files}}
(root/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
(root/'cf110-native14-inputs.jsonl').write_text(''.join(json.dumps(r,separators=(',',':'))+'\n' for r in inputs if r['focal_native_pair']))
(root/'cf110-before/native14-outputs.jsonl').write_text(''.join(json.dumps(r,separators=(',',':'))+'\n' for r in focal))
print(json.dumps(manifest))
for r in focal: print(r['case'],r['output']['values'].get('status'),r['solutions'])
