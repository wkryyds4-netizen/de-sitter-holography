"""Verify saved response evidence and seal once; no scientific reruns."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,re,subprocess,math

ROOT=Path(__file__).resolve().parents[1]
STAGE='susskind-selected-rindler-response'
OUT=ROOT/f'results/{STAGE}/local-evidence.json'
if OUT.exists():raise FileExistsError(OUT)
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text(encoding='utf-8-sig'))
science=f'results/{STAGE}/response-checks.json'
assert sha(science)=='814039c5cc5e3366db8d6957bfef66e80efee02adf409c7e1df0ac59e587fee0'
data=read(science)
assert data['status']=='passed' and data['laurent_cancellations']==36
assert len(data['exact_taylor'])==8 and len(data['rows'])==4
for p,h in data['input_sha256'].items():assert sha(p)==h,p
for row in data['rows']:
    h=row['h'];L=row['dimensions'][0]
    for point in row['values']:
        z=2*abs(point['u'])/math.sqrt(h*(-math.expm1(-h)))
        assert z*z<=L*(L-1)
        if point['u']:assert point['tail_error_log10']<-25
parents={};verified={};historical={}
mutable={'README.md','NEXT.md','PLAN.md','PROJECT.md','AGENTS.md'}
def verify(p,h,origin):
    if isinstance(h,dict):
        for sub,digest in h.items():verify(sub,digest,origin)
    elif p in mutable or p.startswith('journal/'):
        historical.setdefault(p,[]).append(dict(origin=origin,recorded_sha256=h))
    else:
        assert sha(p)==h,(origin,p)
        verified[p]=h
for stage in ['susskind-sector-bilocal','susskind-slater-selected-state']:
    p=f'results/{stage}/local-evidence.json';d=read(p)
    parents[p]=sha(p);verified[p]=sha(p)
    for key in ['hashes','proof_inputs','verified_parent_scientific_hashes','scientific_output_hashes','input_hashes_verified']:
        for path,h in d.get(key,{}).items():verify(path,h,p)
proofs=['notes/susskind-translations.md','notes/susskind-bulk-target.md',
        'notes/susskind-rindler-packet.md','notes/susskind-bilocal-kms-polarization.md',
        'code/compute_slater_selected_curvature.py']
pdf=f'output/pdf/de-sitter-{STAGE}.pdf'
build=f'tmp/pdfs/{STAGE}-build';render=f'tmp/pdfs/{STAGE}-render'
assert sha(pdf)==sha(f'{build}/main.pdf')
log=(ROOT/f'results/{STAGE}/main.log').read_text(errors='replace')
assert 'Output written on' in log and not re.search('Overfull|undefined|Rerun to get|Fatal error',log)
passes={}
for i in (1,2):
    p=f'{build}/pass{i}.txt';raw=(ROOT/p).read_bytes()
    text=raw.decode('utf-16') if raw.startswith((b'\xff\xfe',b'\xfe\xff')) else raw.decode()
    assert 'Output written on' in text and 'Fatal error' not in text
    passes[p]=sha(p)
info=subprocess.check_output(['pdfinfo',str(ROOT/pdf)]).decode(errors='replace')
assert int(re.search(r'Pages:\s+(\d+)',info).group(1))==5
renders={}
for i in range(1,6):
    p=f'{render}/page-{i}.png';assert (ROOT/p).read_bytes().startswith(b'\x89PNG')
    renders[p]=sha(p)
files=['code/selected_rindler_response.py','code/check_selected_rindler_response.py',
       'code/record_selected_rindler_response.py','code/build_selected_rindler_response.ps1',
       f'manuscript/{STAGE}/main.tex',f'notes/{STAGE}.md',f'literature/{STAGE}-sources.md',
       science,pdf,f'results/{STAGE}/main.log']
result=dict(status='Exact source identity and tested semiclassical prediction recorded; full goal active',
    created_at_utc=datetime.now(timezone.utc).isoformat(),
    exact_scope='Selected density equals normalized-source second response; charged restriction uses the physical Hamiltonian, full tower uses e-W. Exact finite Taylor coefficients checked through degree 14.',
    asymptotic_scope='R(u)=sech(u)^4 [u tanh(u)-5 tanh(u)^2] is a semiclassical prediction conditional on convergence of the twice-differentiated thermal saddle. Uniform differentiated remainder not proved.',
    scientific_output_hashes={science:sha(science)},input_hashes_verified=data['input_sha256'],
    environment=data['environment'],parameters=data['parameters'],
    checks=dict(taylor_coefficients=8,laurent_cancellations=36,finite_deformation_times=28),
    previous_evidence=parents,verified_parent_scientific_hashes=verified,
    historical_mutable_entrances=historical,proof_inputs={p:sha(p) for p in proofs},
    hashes={p:sha(p) for p in files},render_hashes=renders,build_pass_hashes=passes,
    build=dict(command='./code/build_selected_rindler_response.ps1',passes=2,pages=5,
               render='pdftoppm -r 90 -png',no_undefined_or_overfull=True,
               visual_review='All five final pages viewed; formulas, table, citations and page boundaries checked.'),
    remaining='Differentiated q-Mehler saddle estimate, mixed e/W sources, microscopic spatial/mass generators and finite-G dynamics.',
    git_commit_created=False,git_push_performed=False,external_reply_sent=False,
    recorder_scope='Saved-byte and build verification. Floating calculations are empirical diagnostics, not interval certificates; no old driver or recorder executed.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',current_files=len(files),parent_scientific_files=len(verified),sha256=sha(OUT.relative_to(ROOT)))))
