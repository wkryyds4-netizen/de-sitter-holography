"""One-time artifact/evidence seal; never runs a scientific checker."""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,re,subprocess
ROOT=Path(__file__).resolve().parents[1]
STAGE='susskind-multimatter-moving-frame'
OUT=ROOT/f'results/{STAGE}/local-evidence.json'
if OUT.exists():raise FileExistsError(OUT)
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text(encoding='utf-8-sig'))
science=f'results/{STAGE}/checks.json'
assert sha(science)=='4f3621ca842f4fac8db04b276567991204fb3ac65bf825053ee0500fddf084db'
d=read(science);assert d['status']=='passed'
assert sum(d['exact_counts'].values())==1764
assert len(d['rows'])==12 and len(d['precision_60_decimal'])==8
for p,h in d['input_sha256'].items():assert sha(p)==h,p
parent='results/susskind-chord-module-frame/local-evidence.json'
assert sha(parent)=='e53e52972b31177f8fb000197d6db8cde04d7f8b369ab4c0913aa3c24f430d50'
verified={parent:sha(parent)};historical={}
mutable={'README.md','PLAN.md','NEXT.md','PROJECT.md','AGENTS.md'}
def verify(p,h):
 if isinstance(h,dict):
  for sub,digest in h.items():verify(sub,digest)
 elif p in mutable or p.startswith('journal/'):historical[p]=h
 else:
  assert sha(p)==h,p
  verified[p]=h
for key in ('hashes','verified_parent_scientific_hashes','scientific_output_hashes','input_hashes_verified','proof_inputs'):
 for p,h in read(parent).get(key,{}).items():verify(p,h)
pdf=f'output/pdf/de-sitter-{STAGE}.pdf'
build=f'tmp/pdfs/{STAGE}-build';render=f'tmp/pdfs/{STAGE}-render'
assert sha(pdf)==sha(f'{build}/main.pdf')
log=(ROOT/f'results/{STAGE}/main.log').read_text(errors='replace')
assert 'Output written on' in log and not re.search('undefined|Overfull|Rerun to get|Fatal error',log)
info=subprocess.check_output(['pdfinfo',str(ROOT/pdf)]).decode(errors='replace')
assert int(re.search(r'Pages:\s+(\d+)',info).group(1))==6
passes={};renders={}
for i in (1,2):
 p=f'{build}/pass{i}.txt';raw=(ROOT/p).read_bytes()
 txt=raw.decode('utf-16') if raw.startswith((b'\xff\xfe',b'\xfe\xff')) else raw.decode()
 assert 'Output written on' in txt and 'Fatal error' not in txt
 passes[p]=sha(p)
for i in range(1,7):
 p=f'{render}/page-{i}.png';assert (ROOT/p).read_bytes().startswith(b'\x89PNG')
 renders[p]=sha(p)
files=['code/multimatter_moving_frame.py','code/check_multimatter_moving_frame.py',
 'code/build_multimatter_moving_frame.ps1','code/record_multimatter_moving_frame.py',
 f'manuscript/{STAGE}/main.tex',f'notes/{STAGE}.md',f'literature/{STAGE}-sources.md',
 pdf,f'results/{STAGE}/main.log',science]
result=dict(status='Ordered two-matter moving finite-word addition established; full research goal active',
 created_at_utc=datetime.now(timezone.utc).isoformat(),
 scope=d['scope'],limits='Fixed-lambda DS first; epsilon->0, m->infinity, epsilon*m->0; fixed internal levels and finite words. No identical-matter completion or mass selection.',
 scientific_output_hashes={science:sha(science)},input_hashes_verified=d['input_sha256'],
 exact_counts=d['exact_counts'],numerical_summary=dict(
  moving_cases=len(d['rows']),
  max_word_refinement=max(r['word_refinement_absolute'] for r in d['rows']),
  max_input_tail_squared_bound=max(r['input_tail_squared_bound'] for r in d['rows']),
  max_60_decimal_relative=max(float(r['relative']) for r in d['precision_60_decimal']),
  max_cancellation={k:max(r['commutators'][k]['cancellation'] for r in d['rows']) for k in ('BP','TQ')}),
 precision_scope='60 decimal sparse outer-ladder composition at m=128,j=2,k1=1,k2=2; 8 quantities including BP-iE and TQ-iQ defects.',
 limitations=d['limitations'],environment=d['environment'],
 sources='HVX2511.08743v2 factorization and swapping formula added as prior work; exact read regions in literature note. Existing two-matter and one-matter scientific inputs unchanged.',
 previous_evidence={parent:sha(parent)},verified_parent_scientific_hashes=verified,
 historical_mutable_entrances=historical,hashes={p:sha(p) for p in files},
 build_pass_hashes=passes,render_hashes=renders,
 build=dict(command='./code/build_multimatter_moving_frame.ps1',passes=2,pages=6,
  render='pdftoppm -r 90 -png',visual_review='All six final pages viewed: equations, table, references, footnote and margins checked.',no_undefined_or_overfull=True),
 remaining='Crossed physical metric and singular odd-channel normalization in the moving frame, actual identical-singlet four-point dynamics; selected ordering, opposite momentum/mass, four-bilocal and finite-G.',
 git_commit_created=False,git_push_performed=False,external_reply_sent=False)
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',current_files=len(files),parent_scientific_files=len(verified),
 pdf_sha256=sha(pdf),evidence_sha256=sha(OUT.relative_to(ROOT)))))
