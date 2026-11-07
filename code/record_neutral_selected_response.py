"""Read-only verification and one-time sealing of the neutral response stage."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,re,subprocess
ROOT=Path(__file__).resolve().parents[1]
STAGE='susskind-neutral-selected-response'
OUT=ROOT/f'results/{STAGE}/local-evidence.json'
if OUT.exists():raise FileExistsError(OUT)
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text(encoding='utf-8-sig'))
science={
 f'results/{STAGE}/response-checks.json':'475535bab4df2947c1de0d886ae42672e478becd612ffb66338afc14b09af843',
 f'results/{STAGE}/spectral-bounds.json':'f5cd67a57050f1564c25e1ceb89a25f26a47acac561945eb5d656bc1a42a9991',
 f'results/{STAGE}/mixed-source-checks.json':'29754be4bb5b8673cc1f76c96c07684ecfc543987f5ade90039f42cb64283b51'}
inputs={};data=[]
for p,h in science.items():
    assert sha(p)==h,p
    d=read(p);assert d['status']=='passed'
    for src,digest in d['input_sha256'].items():
        assert sha(src)==digest,src
        inputs[src]=digest
    data.append(d)
assert sum(x['three_block_identities'] for x in data[0]['exact'])==330
assert sum(x['identities'] for x in data[0]['selected_transfer'])==168
assert len(data[0]['independent'])==14
assert len(data[1]['rows'])==56
assert sum(x['identities'] for x in data[2]['exact'])==90
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
for stage in ['susskind-weighted-mehler-limit','susskind-sector-bilocal','susskind-slater-selected-state','susskind-slater-orbit-density']:
    p=f'results/{stage}/local-evidence.json';d=read(p)
    parents[p]=sha(p);verified[p]=sha(p)
    for key in ['hashes','proof_inputs','verified_parent_scientific_hashes','scientific_output_hashes','input_hashes_verified']:
        for src,h in d.get(key,{}).items():verify(src,h,p)
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
assert int(re.search(r'Pages:\s+(\d+)',info).group(1))==8
renders={}
for i in range(1,9):
    p=f'{render}/page-{i}.png';assert (ROOT/p).read_bytes().startswith(b'\x89PNG')
    renders[p]=sha(p)
files=['code/neutral_selected_response.py','code/check_neutral_selected_response.py',
       'code/check_neutral_mixed_source.py','code/record_neutral_spectral_bounds.py',
       'code/record_neutral_selected_response.py','code/build_neutral_selected_response.ps1',
       f'manuscript/{STAGE}/main.tex',f'notes/{STAGE}.md',f'literature/{STAGE}-sources.md',
       *science,pdf,f'results/{STAGE}/main.log']
result=dict(status='Exact neutral selected response and compact-time asymptotics proved in the manuscript; full research goal active',
 created_at_utc=datetime.now(timezone.utc).isoformat(),
 exact_scope='S(u,v)=tilted_d1(u-v)+sqrt(q)*lambda*(1+q)*B(u,v); B is a two-step Mehler Gram kernel. Odd source gives sigma({h,U(u)})=-sqrt(q*lambda)*H(u). Actual first neutral slope interface inherited.',
 asymptotic_scope='Fixed compact string times and finite derivatives: S=g(u-v)+epsilon/2*(R1(u-v)+g(u)g(v))+O(epsilon^2); odd mean=-2sqrt(epsilon)g(u)+O(epsilon^(3/2)). Equal trace/selected mixture retains epsilon*g(u)g(v)/4.',
 limits='Finite-N exact pinching, half-filled fixed-lambda DS, then epsilon=4lambda to zero. No N-uniform, growing-time, complete-bilocal, mass or finite-G claim.',
 scientific_output_hashes=science,input_hashes_verified=inputs,
 environment=[d['environment'] for d in data],parameters=data[0]['parameters'],
 checks=dict(three_block=330,selected_even=168,selected_odd=90,independent_angular=14,two_time=28,odd_time=20,analytic_bounds=56),
 previous_evidence=parents,verified_parent_scientific_hashes=verified,historical_mutable_entrances=historical,
 hashes={p:sha(p) for p in files},build_pass_hashes=passes,render_hashes=renders,
 build=dict(command='./code/build_neutral_selected_response.ps1',passes=2,pages=8,render='pdftoppm -r 90 -png',
            visual_review='All eight final pages viewed after odd-source addition; equations, references, page boundaries and layout checked.',
            no_undefined_or_overfull=True),
 recorder_scope='Saved-byte and PDF verification only. No successful scientific driver or prior recorder rerun. Exact-arithmetic analytic spectral bounds do not include roundoff; angular quadrature not interval certified.',
 remaining='Higher neutral residuals and charged-neutral cross terms through the complete bilocal filter, spatial/mass generators and finite-G; arbitrary selected interlacings remain open.',
 git_commit_created=False,git_push_performed=False,external_reply_sent=False)
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',current_files=len(files),parent_scientific_files=len(verified),sha256=sha(OUT.relative_to(ROOT)))))
