"""One-time read-only evidence seal. Does not run scientific drivers."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib,json,re,subprocess
ROOT=Path(__file__).resolve().parents[1]
STAGE='susskind-thermal-chord-window'
OUT=ROOT/f'results/{STAGE}/local-evidence.json'
if OUT.exists():raise FileExistsError(OUT)
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text(encoding='utf-8-sig'))
science={f'results/{STAGE}/window-checks.json':'d36161666759fda1bfe351a0a09be636692931705fef7facdecf27a10a4cd6b0',
 f'results/{STAGE}/precision-checks.json':'a5a06530821b404234e0ebc57bac299b6ae63513fc21c263dc5569e094d2e94c',
 f'results/{STAGE}/tail-bounds.json':'0875b6f3d3475d3e521211ed13d1ad64444b293ca4fbca8feb34ab8d2e6f7939'}
for p,h in science.items():
 assert sha(p)==h,p
 assert read(p)['status']=='passed'
 for src,digest in read(p)['input_sha256'].items():assert sha(src)==digest,src
d=read(f'results/{STAGE}/window-checks.json')
assert sum(d['exact_counts'].values())==4152 and len(d['rows'])==4
prec=read(f'results/{STAGE}/precision-checks.json')
tails=read(f'results/{STAGE}/tail-bounds.json')
assert len(prec['rows'])==len(tails['rows'])==4
verified={};historical={}
mutable={'README.md','PLAN.md','NEXT.md','PROJECT.md','AGENTS.md'}
def verify(p,h,origin):
 if isinstance(h,dict):
  for sub,digest in h.items():verify(sub,digest,origin)
 elif p in mutable or p.startswith('journal/'):
  historical.setdefault(p,[]).append(dict(origin=origin,recorded_sha256=h))
 else:
  assert sha(p)==h,(origin,p)
  verified[p]=h
parent='results/susskind-preparation-modular-clock/local-evidence.json'
assert sha(parent)=='6b7a75c24e923934b3f63a5e6b8ab24eefa14cb96204ee57fdc73e9e89423e6d'
verified[parent]=sha(parent)
for key in ['hashes','proof_inputs','verified_parent_scientific_hashes',
            'scientific_output_hashes','input_hashes_verified']:
 for p,h in read(parent).get(key,{}).items():verify(p,h,parent)
pdf=f'output/pdf/de-sitter-{STAGE}.pdf'
build=f'tmp/pdfs/{STAGE}-build';render=f'tmp/pdfs/{STAGE}-render'
assert sha(pdf)==sha(f'{build}/main.pdf')
log=(ROOT/f'results/{STAGE}/main.log').read_text(errors='replace')
assert 'Output written on' in log and not re.search('undefined|Overfull|Rerun to get|Fatal error',log)
info=subprocess.check_output(['pdfinfo',str(ROOT/pdf)]).decode(errors='replace')
assert int(re.search(r'Pages:\s+(\d+)',info).group(1))==6
passes={}
for i in (1,2):
 p=f'{build}/pass{i}.txt';raw=(ROOT/p).read_bytes()
 txt=raw.decode('utf-16') if raw.startswith((b'\xff\xfe',b'\xfe\xff')) else raw.decode()
 assert 'Output written on' in txt and 'Fatal error' not in txt
 passes[p]=sha(p)
renders={}
for i in range(1,7):
 p=f'{render}/page-{i}.png';assert (ROOT/p).read_bytes().startswith(b'\x89PNG')
 renders[p]=sha(p)
files=['code/thermal_chord_window.py','code/check_thermal_chord_window.py',
 'code/check_thermal_chord_precision.py','code/record_thermal_chord_tail_bounds.py',
 'code/build_thermal_chord_window.ps1','code/record_thermal_chord_window.py',
 f'manuscript/{STAGE}/main.tex',f'notes/{STAGE}.md',f'literature/{STAGE}-sources.md',
 pdf,f'results/{STAGE}/main.log',*science]
diagnostics=[str(p.relative_to(ROOT)).replace('\\','/') for p in sorted((ROOT/f'results/{STAGE}/diagnostics').glob('*'))]
files+=diagnostics
result=dict(status='Thermal endpoint joint window and single-action boost matching established; overall research goal active',
 created_at_utc=datetime.now(timezone.utc).isoformat(),
 scope='Fixed-lambda tracial one-insertion DS module followed by epsilon->0, m->infinity, epsilon*m^2->0. Actual charged curvature is d=2; full selected rho not inserted.',
 scientific_output_hashes=science,exact_counts=d['exact_counts'],
 corrections=dict(longdouble=prec['correction'],tail_underflow=tails['correction'],
   scipy_expm='Long and split expm diagnostics did not provide the requested reference precision; use positive Taylor and supplemental 50-decimal values.'),
 numerical_summary=dict(max_decimal_relative_difference=max(v for r in prec['rows'] for v in r['relative_difference'].values()),
   max_split_expm_norm_difference=max(r['split_expm_norm_difference'] for r in d['rows']),
   max_d2_boost_cancellation=max(r['observables']['2']['boost_sum_cancellation'] for r in d['rows']),
   jacobi_cutoff_bound=max(float(r['normalized_cutoff_vector_bound']) for r in tails['rows']),
   scaled_taylor_bound=max(float(r['scaled_vector_taylor_tail_bound']) for r in tails['rows'])),
 environment=d['environment'],previous_evidence={parent:sha(parent)},
 verified_parent_scientific_hashes=verified,historical_mutable_entrances=historical,
 hashes={p:sha(p) for p in files},build_pass_hashes=passes,render_hashes=renders,
 build=dict(command='./code/build_thermal_chord_window.ps1',passes=2,pages=6,
   render='pdftoppm -r 90 -png',visual_review='All six final pages viewed; equations, table, references, footnote and margins checked.',
   no_undefined_or_overfull=True),
 remaining='Moving-space map on the full matter module; selected rho mixed generator matrix elements and descendant commutator norms; null translations, mass, four-bilocal and finite-G dynamics.',
 git_commit_created=False,git_push_performed=False,external_reply_sent=False)
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',current_files=len(files),parent_scientific_files=len(verified),
 pdf_sha256=sha(pdf),evidence_sha256=sha(OUT.relative_to(ROOT)))))
