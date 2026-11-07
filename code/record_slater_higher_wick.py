"""Read-only saved evidence verifier; exclusive creation of this stage manifest."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,re,subprocess
ROOT=Path(__file__).resolve().parent.parent
STAGE='susskind-slater-higher-wick';DIR=ROOT/'results'/STAGE
OUT=DIR/'local-evidence.json'
if OUT.exists():raise FileExistsError(OUT)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def rel(p):return p.relative_to(ROOT).as_posix()
science=DIR/'exact-checks.json'
assert sha(science)=='27b3557c6ce6c68605eb1a89ecbb7e23ebe7f11d00fc40e87f84f8425f845143'
data=read(science);assert data['status']=='passed'
for p,h in data['input_sha256'].items():assert sha(ROOT/p)==h,p
parentpath=ROOT/'results/susskind-slater-selected-state/local-evidence.json'
# No parent recorder is imported.
assert sha(parentpath)=='d387e521347c020daa48b4dde8c1291be5d308a71caed4c89a7c637566524f01'
parent=read(parentpath);verified={}
for key in ['hashes','verified_parent_scientific_hashes']:
    for p,h in parent[key].items():assert sha(ROOT/p)==h,p;verified[p]=h
p=parent['parent_manifest'];assert sha(ROOT/p['path'])==p['sha256'];verified[p['path']]=p['sha256']
pdf=ROOT/f'output/pdf/de-sitter-{STAGE}.pdf'
build=ROOT/f'tmp/pdfs/{STAGE}-build';render=ROOT/f'tmp/pdfs/{STAGE}-render'
assert sha(pdf)==sha(build/'main.pdf')
log=DIR/'main.log';logtext=log.read_text(errors='replace')
assert 'Output written on' in logtext and not re.search(r'undefined|Overfull|Fatal error|Rerun to get',logtext)
passes={}
for i in (1,2):
    p=build/f'pass{i}.txt';raw=p.read_bytes()
    txt=raw.decode('utf-16') if raw.startswith((b'\xff\xfe',b'\xfe\xff')) else raw.decode(errors='replace')
    assert 'Output written on' in txt and 'Fatal error' not in txt;passes[rel(p)]=sha(p)
info=subprocess.check_output(['pdfinfo',str(pdf)]).decode(errors='replace')
assert int(re.search(r'Pages:\s+(\d+)',info).group(1))==6
renders={}
for i in range(1,7):
    p=render/f'page-{i}.png';assert p.read_bytes().startswith(b'\x89PNG');renders[rel(p)]=sha(p)
proofs=['notes/susskind-memory-joint-regression-proof.md',
        'notes/susskind-bilocal-fourth-tail-diagram-draft.md',
        'notes/susskind-bilocal-fourth-tail-operator-draft.md',
        'notes/susskind-bilocal-fourth-tail-charged-proof.md',
        'literature/susskind-bilocal-clifford-riesz-sources.md']
files=[Path(__file__),ROOT/'code/slater_higher_wick_sources.py',
       ROOT/'code/check_slater_higher_wick_sources.py',ROOT/'code/build_slater_higher_wick.ps1',
       ROOT/f'notes/{STAGE}.md',ROOT/f'literature/{STAGE}-sources.md',
       ROOT/f'manuscript/{STAGE}/main.tex',pdf,log,science,
       DIR/'diagnostics/initial-build-failure.txt',DIR/'diagnostics/initial-recorder-failure.txt']
result=dict(status='Local research milestone; complete de Sitter goal remains active',
    created_at_utc=datetime.now(timezone.utc).isoformat(),
    scientific_scope='All fixed disjoint Wick marks/subblocks and chaos constraints: boundary density identity via finite local reduction; charged tower f(h), fixed-degree central spectral/pinching transfer. Ordinary internal-pair central bridge and uniform actual tails remain open.',
    scientific_output_hashes={rel(science):sha(science)},input_hashes_verified=data['input_sha256'],
    parameters=data['parameters'],checks=data['counts'],density_cases=len(data['density_rows']),
    legacy_comparisons=data['legacy_comparisons'],
    parent_manifest={'path':rel(parentpath),'sha256':sha(parentpath)},
    verified_parent_scientific_hashes=verified,
    historical_mutable_entrances=parent['historical_mutable_entrances'],
    proof_inputs={p:sha(ROOT/p) for p in proofs},hashes={rel(p):sha(p) for p in files},
    render_hashes=renders,build_pass_hashes=passes,
    build={'command':'./code/build_slater_higher_wick.ps1','passes':2,'pages':6,
           'no_undefined_or_overfull':True,'render':'pdftoppm -r 90 -png',
           'visual_review':'All six final pages viewed at 90 dpi; no clipping, overlapping formulas or unresolved references observed.'},
    environment=data['environment'],
    diagnostics='Fresh scientific checker passed its first full run. Initial LaTeX pass failed on a missing backslash before left; failure log preserved, corrected and two passes completed. Initial recorder stopped before output on UTF-16 PowerShell pass logs; BOM-aware reading corrected without changing logs. No successful scientific output or parent manifest rewritten.',
    git_commit_created=False,git_push_performed=False,external_reply_sent=False,
    recorder_scope='Saved-byte integrity and build verification; general conclusions also require the manuscript reduction and cited prior proofs.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({'status':'passed','current_files':len(files),'parent_scientific_files':len(verified),'sha256':sha(OUT)}))
