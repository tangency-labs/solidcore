import json,math,argparse,sys
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
parser=argparse.ArgumentParser(description='Audit reflected offset metadata from frozen defining inputs; no solver tolerance.')
parser.add_argument('--native',type=Path,default=Path(__file__).with_name('projected_offset_reversal_native_prior.jsonl'))
parser.add_argument('--require-success',action='store_true')
a=parser.parse_args()
base=Path(__file__).parent
inputs=list(map(json.loads,(base/'projected_offset_reversal_inputs.jsonl').read_text(encoding='utf-8-sig').splitlines()))
changes={r['case']:r for r in map(json.loads,(base/'projected_offset_reversal_before_prior.jsonl').read_text(encoding='utf-8-sig').splitlines())}
native={r['case']:r for r in map(json.loads,a.native.read_text(encoding='utf-8-sig').splitlines())}
proof=[]
for r in inputs:
 name=r['case'];i=r['input'];record={'case':name,'before':changes[name]['before'],'prior_parent':changes[name]['after'],'current':{k:native[name]['native'][k] for k in ['status','out']}};assert native[name]['native']['status']==record['current']['status'];assert native[name]['native']['out']==record['current']['out']
 if i['op']=='continuities':
  before=record['before']['out'];out=before[:2]+before[2:][::-1];assert out==record['current']['out'];record['expected_metadata']={'status':[0],'out':out};record['proof']='same indexed continuity labels in reverse parameter order; absolute inherited nonzero-offset labels are not re-certified';record['index_correspondence']=list(range(int(before[1])-1,-1,-1))
 else:
  raw=i['base'].get('bspline',{}).get('u',[0.,0.,0.,1.,1.,1.]);pairs=sorted(Counter(raw).items());layers=[];bad=[]
  for s in [i['spec']]+([i['outer']] if i['outer'] else []):
   if s['nrev']%2:
    S=s['a']+s['b'];mapped=[]
    for k,m in pairs[::-1]:
     t=S-k;mapped.append((t,m))
     if S-t!=k:
      if S-t>k: left=t;right=math.nextafter(t,math.inf)
      else:left=math.nextafter(t,-math.inf);right=t
      assert math.nextafter(left,math.inf)==right
      assert S-left>k>S-right
      bad.append({'defining_knot':k,'rounded_reflection':t,'roundtrip':S-t,'adjacent_parameters':[left,right],'mapped_adjacent':[S-left,S-right],'exact_inverse_fraction':str(F(S)-F(k))})
    layers.append({'sum':S,'source_pairs':pairs,'reflected_pairs':mapped});pairs=mapped
  request=i.get('ivl');mapped_query=None
  if request:
   lo,hi=sorted(request);pairs=[(k,m) for k,m in pairs if lo<=k<=hi]
   s=i['outer'] or i['spec'];S=s['a']+s['b'];mapped_query=sorted([S-lo,S-hi]) if s['nrev']%2 else [lo,hi]
  out=[float(len(pairs))]+[k for k,m in pairs]+[float(len(pairs)) if i['mults'] else 0.]+([float(m) for k,m in pairs] if i['mults'] else [])
  record.update(expected_metadata={'status':[0],'out':out},layers=layers,requested_interval=request,immediate_base_query=mapped_query,nonrepresentable_breaks=bad)
  if record['current']['status']==[0]:assert out==record['current']['out']
  else:assert bad and record['current']['status']==[1007];record['current_is_api_regression']=True
 proof.append(record)
if a.require_success:assert all(r['current']['status']==[0] for r in proof)
print('\n'.join(map(json.dumps,proof)))
print('Verified 19 records; broad metadata refusals remaining: '+str(sum(r['current']['status']!=[0] for r in proof)),file=sys.stderr)
