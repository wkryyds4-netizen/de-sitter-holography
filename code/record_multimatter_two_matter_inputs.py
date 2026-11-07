"""Read-only supplemental verification of the separate November 6 branch."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/susskind-multimatter-moving-frame/prior-two-matter-evidence.json'
if OUT.exists():raise FileExistsError(OUT)
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text(encoding='utf-8-sig'))
parent='results/susskind-two-matter/local-evidence.json';manifest=read(parent)
verified={parent:sha(parent)};historical={}
for p,h in manifest['hashes'].items():
 if p in ('README.md','PLAN.md','NEXT.md','PROJECT.md','AGENTS.md') or p.startswith('journal/'):
  historical[p]=h
 else:
  assert sha(p)==h,p
  verified[p]=h
science=read('results/susskind-two-matter/validation.json')
assert science['status']=='passed' and len(science['symbolic']['checks'])==38
for p,h in science['input_sha256'].items():
 assert sha(p)==h,p
 verified[p]=h
tensor=read('results/susskind-two-matter/tensor-check.json')
assert tensor['status']=='passed' and len(tensor['checks'])==3
assert sha('code/check_two_matter_tensor.py')==tensor['input_sha256']
current='results/susskind-multimatter-moving-frame/local-evidence.json'
assert sha(current)=='8e0bd77324ef54c01934a200862d271a473279d7b16d43582329c0025f0c4b7d'
result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),
 reason='The November 6 two-matter branch is not transitively included in the recent moving-frame manifest; verify its separate sealed science and proof inputs explicitly.',
 current_evidence={current:sha(current)},verified_prior_hashes=verified,
 historical_mutable_entrances=historical,
 recorder_sha256=sha(Path(__file__).relative_to(ROOT)),
 checks=41,scientific_drivers_rerun=False)
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',verified_files=len(verified),sha256=sha(OUT.relative_to(ROOT)))))
