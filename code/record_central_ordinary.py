"""Saved provenance verification only; no successful scientific driver imports."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,re,subprocess
ROOT=Path(__file__).resolve().parent.parent
STAGE='susskind-central-ordinary';DIR=ROOT/'results'/STAGE;OUT=DIR/'local-evidence.json'
if OUT.exists():raise FileExistsError(OUT)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def rel(p):return p.relative_to(ROOT).as_posix()
expected={'exact-checks.json':'54286ec6f07967d3c60770c44ffd93799dd4eac26cb346dbdcab49165280f45c',
 'local-coverage.json':'e7321618b31f4bffe489a4f2eb61f379ee56e443dcfc84dbb912c372195b7add',
 'fourier-checks.json':'eceb91cd3d59939dad487a968f5a8370ade15aaff8d7690bbe317662a27961f3'}
science={};inputs={}
for name,h in expected.items():
    p=DIR/name;assert sha(p)==h;data=read(p);assert data['status']=='passed';science[name]=data
    for path,digest in data['input_sha256'].items():assert sha(ROOT/path)==digest,path;inputs[path]=digest
assert all(science['local-coverage.json']['counts'][f'{r}_source_cases']>0 for r in (2,3,4))
assert science['exact-checks.json']['counts'].get('site_vs_geometric_central_sources',0)==0
parentpath=ROOT/'results/susskind-slater-higher-wick/local-evidence.json'
assert sha(parentpath)=='14c5583cca8d9e6b25ff3cccacfc7157f6ef98a33424519a891ba93519417918'
parent=read(parentpath);verified={}
for key in ['hashes','verified_parent_scientific_hashes']:
    for path,h in parent[key].items():assert sha(ROOT/path)==h,path;verified[path]=h
v=parent['parent_manifest'];assert sha(ROOT/v['path'])==v['sha256'];verified[v['path']]=v['sha256']
pdf=ROOT/f'output/pdf/de-sitter-{STAGE}.pdf';log=DIR/'main.log'
build=ROOT/f'tmp/pdfs/{STAGE}-build';render=ROOT/f'tmp/pdfs/{STAGE}-render'
assert sha(pdf)==sha(build/'main.pdf')
txt=log.read_text(errors='replace')
assert 'Output written on' in txt and not re.search(r'undefined|Overfull|Fatal error|Rerun to get',txt)
passes={}
for i in (1,2):
    p=build/f'pass{i}.txt';raw=p.read_bytes()
    txt=raw.decode('utf-16') if raw.startswith((b'\xff\xfe',b'\xfe\xff')) else raw.decode(errors='replace')
    assert 'Output written on' in txt and 'Fatal error' not in txt;passes[rel(p)]=sha(p)
info=subprocess.check_output(['pdfinfo',str(pdf)]).decode(errors='replace')
assert int(re.search(r'Pages:\s+(\d+)',info).group(1))==6
before=read(ROOT/'tmp/pdfs/susskind-central-ordinary-before-clarification.json');renders={}
changed=[]
for i in range(1,7):
    p=render/f'page-{i}.png';assert p.read_bytes().startswith(b'\x89PNG');renders[rel(p)]=sha(p)
    if before[p.name]!=sha(p):changed.append(i)
assert changed==[2,3]
proofs=['notes/susskind-memory-charge-tilt-proof.md','notes/susskind-memory-joint-regression-proof.md',
 'notes/susskind-bilocal-fourth-tail-charged-proof.md','notes/susskind-bilocal-fourth-tail-operator-draft.md',
 'notes/susskind-bilocal-fourth-tail-diagram-draft.md','literature/susskind-bilocal-clifford-riesz-sources.md']
files=[Path(__file__),ROOT/'code/central_ordinary_sources.py',ROOT/'code/check_central_ordinary_sources.py',
 ROOT/'code/check_central_ordinary_local.py',ROOT/'code/compute_central_ordinary_fourier.py',
 ROOT/'code/build_central_ordinary.ps1',ROOT/f'notes/{STAGE}.md',ROOT/f'literature/{STAGE}-sources.md',
 ROOT/f'manuscript/{STAGE}/main.tex',pdf,log,DIR/'diagnostics/initial-local-coverage-gap.txt']+[DIR/n for n in expected]
result=dict(status='Local research milestone; complete de Sitter goal remains active',
 created_at_utc=datetime.now(timezone.utc).isoformat(),
 scientific_scope='Exact half-charge Fourier conditioning of ordinary D source blocks, all fixed central moments, same-copy ordinary-to-Wick strong bridge in fixed Schatten norms, and transfer of fixed exactly pinched coefficients to selected preparations. Uniform actual degree tail remains open.',
 scientific_output_hashes={rel(DIR/n):h for n,h in expected.items()},input_hashes_verified=inputs,
 parameters={n:d['parameters'] for n,d in science.items()},
 checks={'exact':science['exact-checks.json']['counts'],'local':science['local-coverage.json']['counts'],'Fourier_rows':len(science['fourier-checks.json']['rows'])},
 parent_manifest={'path':rel(parentpath),'sha256':sha(parentpath)},verified_parent_scientific_hashes=verified,
 historical_mutable_entrances=parent['historical_mutable_entrances'],
 proof_inputs={p:sha(ROOT/p) for p in proofs},hashes={rel(p):sha(p) for p in files},
 render_hashes=renders,build_pass_hashes=passes,
 build={'command':'./code/build_central_ordinary.ps1','passes':2,'pages':6,'no_undefined_or_overfull':True,
        'render':'pdftoppm -r 90 -png','visual_review':'All six pages viewed; normalization/support clarification followed by final two-pass build and reinspection of changed pages 2/3. Other four render hashes unchanged.'},
 environment={n:d['environment'] for n,d in science.items()},
 diagnostics='First checker returned valid Fourier/Gram/scalar checks but zero local layouts; its output is retained unchanged and not counted as local coverage. Separate supplemental checker passed 5220 local cases with explicit nonzero coverage assertions. All three output-producing drivers ran once; no sealed science or parent manifests rewritten. Initial and final PDF builds passed.',
 git_commit_created=False,git_push_performed=False,external_reply_sent=False,
 recorder_scope='Saved-byte integrity, nonvacuous coverage accounting and PDF build checks; not an independent formal proof of the Fourier bounds or infinite-time convergence.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({'status':'passed','current_files':len(files),'parent_scientific_files':len(verified),'sha256':sha(OUT)}))
