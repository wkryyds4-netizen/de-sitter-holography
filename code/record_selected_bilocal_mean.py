"""One-time read-only evidence verification; never reruns scientific drivers."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,re,subprocess
ROOT=Path(__file__).resolve().parents[1]
STAGE='susskind-selected-bilocal-mean'
OUT=ROOT/f'results/{STAGE}/local-evidence.json'
if OUT.exists():raise FileExistsError(OUT)
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text(encoding='utf-8-sig'))
science={f'results/{STAGE}/projection-checks.json':'25b14d70ef08df4a58cf11f9638fafe408f61b5b862b304235263e16bda7099a',
         f'results/{STAGE}/time-checks.json':'4a5c7ddb00b17ec821d46894afe4282dbaa6a1b715a876df81b2451d74263da5'}
inputs={};data=[]
for p,h in science.items():
 assert sha(p)==h,p
 d=read(p);assert d['status']=='passed'
 for src,digest in d['input_sha256'].items():
  assert sha(src)==digest,src
  inputs[src]=digest
 data.append(d)
assert sum(data[0]['counts'].values())==320
assert all(x['passed'] for x in data[0]['rows'])
assert data[1]['symbol_checks']==72 and len(data[1]['numerical'])==12
assert sum(len(r['rows']) for r in data[1]['exact_time'])==10
derived=[]
for row in data[1]['numerical']:
 fine=complex(*row['value']);coarse=complex(*row['coarse']['value'])
 diff=abs(fine-coarse);assert diff==row['refinement_difference']
 derived.append(dict(q=row['q'],t=row['t'],T=row['T'],absolute_difference=diff,
                     relative_difference=diff/max(abs(fine),1e-300),
                     relative_imaginary_part=abs(fine.imag)/max(abs(fine.real),1e-300),
                     exact_arithmetic_degree_tail=row['degree_tail'],cancellation=row['cancellation']))
verified={};parents={};historical={}
mutable={'README.md','NEXT.md','PLAN.md','PROJECT.md','AGENTS.md'}
def verify(p,h,origin):
 if isinstance(h,dict):
  for sub,digest in h.items():verify(sub,digest,origin)
 elif p in mutable or p.startswith('journal/'):
  historical.setdefault(p,[]).append(dict(origin=origin,recorded_sha256=h))
 else:
  assert sha(p)==h,(origin,p)
  verified[p]=h
for stage in ['susskind-neutral-selected-response','susskind-sector-bilocal',
              'susskind-memory-wick-descendants','susskind-slater-higher-wick']:
 p=f'results/{stage}/local-evidence.json';d=read(p)
 parents[p]=sha(p);verified[p]=sha(p)
 for key in ['hashes','proof_inputs','verified_parent_scientific_hashes','scientific_output_hashes','input_hashes_verified']:
  for src,h in d.get(key,{}).items():verify(src,h,p)
proof=['notes/susskind-bilocal-multisource-proof.md','notes/susskind-memory-neutral-marked-series.md',
       'manuscript/susskind-memory-neutral-marked-series/main.tex','notes/susskind-memory-joint-regression-proof.md',
       'notes/susskind-sector-bilocal.md','notes/susskind-memory-wick-descendants.md',
       'notes/susskind-slater-higher-wick.md']
pdf=f'output/pdf/de-sitter-{STAGE}.pdf';build=f'tmp/pdfs/{STAGE}-build';render=f'tmp/pdfs/{STAGE}-render'
assert sha(pdf)==sha(f'{build}/main.pdf')
log=(ROOT/f'results/{STAGE}/main.log').read_text(errors='replace')
assert 'Output written on' in log and not re.search('Overfull|undefined|Rerun to get|Fatal error',log)
passes={}
for i in (1,2):
 p=f'{build}/pass{i}.txt';raw=(ROOT/p).read_bytes()
 txt=raw.decode('utf-16') if raw.startswith((b'\xff\xfe',b'\xfe\xff')) else raw.decode()
 assert 'Output written on' in txt and 'Fatal error' not in txt
 passes[p]=sha(p)
info=subprocess.check_output(['pdfinfo',str(ROOT/pdf)]).decode(errors='replace')
assert int(re.search(r'Pages:\s+(\d+)',info).group(1))==7
renders={}
for i in range(1,8):
 p=f'{render}/page-{i}.png';assert (ROOT/p).read_bytes().startswith(b'\x89PNG');renders[p]=sha(p)
files=['code/neutral_tower_projection.py','code/check_neutral_tower_projection.py',
       'code/check_selected_bilocal_mean.py','code/build_selected_bilocal_mean.ps1',
       'code/record_selected_bilocal_mean.py',f'manuscript/{STAGE}/main.tex',f'notes/{STAGE}.md',
       f'literature/{STAGE}-sources.md',*science,pdf,f'results/{STAGE}/main.log']
result=dict(status='All-degree single-insertion projection and fixed-lambda selected bilocal mean proved in the manuscript; full research goal active',
 created_at_utc=datetime.now(timezone.utc).isoformat(),
 exact_scope='P T_vw=kappa_v H_(v+w)+sqrt(lambda) sum_(i<v) Wick(h^i U h^(v+w-1-i)); P J_d=0. Lower-chaos recursion and actual bilocal time filter give a convergent two-energy formula for the first moment.',
 short_time='h-time, center T=0: sigma B(0,t)=(2/3)*sqrt(q)*lambda^2*q*(1-q^2)*t^4+O(t^6).',
 limits='Finite-N exact pinching, half-filled fixed-lambda DS, then removal of degree cutoff; no new small-lambda, growing-time, spatial/mass or finite-G result.',
 scientific_output_hashes=science,input_hashes_verified=inputs,
 checks=dict(projection=data[0]['counts'],CKY_polynomials=72,force_commutator_Taylor=10,bilocal_t4=2,physical_time_points=12),
 environment=dict(projection_python=data[0]['python'],time=data[1]['environment']),
 saved_data_relative_analysis=derived,
 max_relative_refinement=max(x['relative_difference'] for x in derived),
 max_degree_tail=max(x['exact_arithmetic_degree_tail'] for x in derived),
 proof_inputs={p:sha(p) for p in proof},previous_evidence=parents,
 verified_parent_scientific_hashes=verified,historical_mutable_entrances=historical,
 hashes={p:sha(p) for p in files},build_pass_hashes=passes,render_hashes=renders,
 build=dict(command='./code/build_selected_bilocal_mean.ps1',passes=2,pages=7,render='pdftoppm -r 90 -png',
            visual_review='All seven final pages viewed after analytic convergence sentence and bibliography location correction; formulas, references, margins and pagination checked.',
            no_undefined_or_overfull=True),
 caveats='Exact-arithmetic degree bound excludes roundoff and energy/time quadrature. Refinement checks are diagnostics. Rational source checks use independent formal q and lambda. Projection is not multiplicative.',
 remaining='Two-insertion selected bilocal fluctuations with connected marks and orthogonal components, then spatial/mass generators and finite-G. Arbitrary selected interlacings remain open.',
 git_commit_created=False,git_push_performed=False,external_reply_sent=False)
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',current_files=len(files),parent_scientific_files=len(verified),
 max_relative_refinement=result['max_relative_refinement'],pdf_sha256=sha(pdf),sha256=sha(OUT.relative_to(ROOT)))))
