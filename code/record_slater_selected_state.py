"""Read-only provenance and saved rational-interval consistency checks."""
from pathlib import Path
from datetime import datetime,timezone
from fractions import Fraction
import hashlib,json,platform,re,subprocess
ROOT=Path(__file__).resolve().parent.parent
STAGE='susskind-slater-selected-state';DIR=ROOT/'results'/STAGE;OUT=DIR/'local-evidence.json'
if OUT.exists():raise FileExistsError(OUT)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def rel(p):return p.relative_to(ROOT).as_posix()
expected={'exact-checks.json':'65d3cac3fe494ab0770ff97df005452c6e3ecdd2cb6e444f713a0d7ae1a5dea6',
          'curvature-dynamics.json':'b45dd55f29af635552289b500c79f4a3301a7932c73414fb5a3845907bbef382',
          'sign-certificate.json':'635e988df224a84157a9a382640901d96c1d96167db3248566f55ab0fb633f2a'}
science={};inputs={}
for name,digest in expected.items():
    p=DIR/name;assert sha(p)==digest;data=read(p);assert data['status']=='passed';science[name]=data
    for path,h in data['input_sha256'].items():assert sha(ROOT/path)==h,path;inputs[path]=h
    for path,h in data.get('figure_sha256',{}).items():assert sha(ROOT/path)==h,path
for row in science['sign-certificate.json']['rows']:
    lower,upper=Fraction(row['lower']),Fraction(row['upper'])
    partial,error=Fraction(row['partial_sum']),Fraction(row['error_bound'])
    assert lower==partial-error and upper==partial+error and error>0
    assert upper<0 if row['state']=='one-selected' else lower>0
parentpath=ROOT/'results/susskind-slater-composite-transfer/local-evidence.json'
assert sha(parentpath)=='d3203a3fa10a17fa72c3eaed96968730fce3f2de28b53a64a4b15529934953c3'
parent=read(parentpath);verified={}
for key in ['hashes','verified_parent_scientific_hashes']:
    for path,h in parent[key].items():assert sha(ROOT/path)==h,path;verified[path]=h
p=parent['parent_manifest'];assert sha(ROOT/p['path'])==p['sha256'];verified[p['path']]=p['sha256']
pdf=ROOT/f'output/pdf/de-sitter-{STAGE}.pdf';log=DIR/'main.log'
build=ROOT/f'tmp/pdfs/{STAGE}-build';render=ROOT/f'tmp/pdfs/{STAGE}-render'
assert sha(pdf)==sha(build/'main.pdf')
text=log.read_text(errors='replace')
assert 'Output written on' in text and not re.search(r'undefined|Overfull|Fatal error|Rerun to get',text)
passes={}
for i in (1,2):
    p=build/f'pass{i}.txt';text=p.read_text(errors='replace')
    assert 'Output written on' in text and 'did not succeed' not in text;passes[rel(p)]=sha(p)
info=subprocess.check_output(['pdfinfo',str(pdf)]).decode(errors='replace')
assert int(re.search(r'Pages:\s+(\d+)',info).group(1))==7
renders={}
for i in range(1,8):
    p=render/f'page-{i}.png';assert p.read_bytes().startswith(b'\x89PNG\r\n\x1a\n');renders[rel(p)]=sha(p)
files=[Path(__file__),ROOT/'code/slater_selected_covariant_state.py',
       ROOT/'code/check_slater_selected_covariant_state.py',ROOT/'code/compute_slater_selected_curvature.py',
       ROOT/'code/certify_slater_selected_curvature.py',ROOT/'code/build_slater_selected_state.ps1',
       ROOT/f'notes/{STAGE}.md',ROOT/f'literature/{STAGE}-sources.md',ROOT/f'manuscript/{STAGE}/main.tex',
       pdf,log,DIR/'selected-curvature.png']+[DIR/name for name in expected]
parameters={name:data['parameters'] for name,data in science.items()}
result=dict(status='Local research milestone; the complete de Sitter research goal remains active',
    created_at_utc=datetime.now(timezone.utc).isoformat(),
    scientific_scope='One-selected-state bounded density on the full h/u/x covariant algebra, stationary h/x restriction, tracial curvature marginal, actual fixed-time pinched-curvature kernel and rational strict-sign distinction at physical q=1/4,t=2. Higher bilocal descendants and the bulk dictionary remain open.',
    scientific_output_hashes={rel(DIR/name):digest for name,digest in expected.items()},
    input_hashes_verified=inputs,parameters=parameters,
    parameter_hash=hashlib.sha256(json.dumps(parameters,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
    parent_manifest=dict(path=rel(parentpath),sha256=sha(parentpath)),verified_parent_scientific_hashes=verified,
    historical_mutable_entrances=parent['historical_mutable_entrances'],
    hashes={rel(p):sha(p) for p in files},render_hashes=renders,build_pass_hashes=passes,
    build=dict(command='./code/build_slater_selected_state.ps1',passes=2,pages=7,
        no_undefined_or_overfull=True,render='pdftoppm -r 90 -png',
        visual_review='All seven pages viewed; finite-N normalization and rotation wording clarified, then final pages 1 and 3 re-viewed. Other five rendered page hashes unchanged. No clipping or formula/citation defects observed.'),
    checks=dict(density=2186,boundary_source=2953,reduced_hj=1022,pure_j=26,
                physical_q=3,time_rows=18,rational_sign_intervals=2),
    environment={name:data['environment'] for name,data in science.items()},
    diagnostics='All three fresh scientific drivers passed on their first full run and were not rerun. Both initial and final two-pass PDF builds passed. Final manuscript wording changes affected pages 1 and 3 only. Successful parent science and manifests were not rewritten.',
    git_commit_created=False,git_push_performed=False,external_reply_sent=False,
    recorder_scope='Saved file/build integrity and rational interval consistency; not an independent proof of general source identities or finite-N convergence.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',current_files=len(files),parent_scientific_files=len(verified),pages=7,sha256=sha(OUT))))
