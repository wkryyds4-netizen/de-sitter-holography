"""Read-only stage manifest; no scientific or plotting generator execution."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,re,subprocess,platform
ROOT=Path(__file__).resolve().parent.parent
STAGE='susskind-slater-stationary'
DIR=ROOT/f'results/{STAGE}'
OUT=DIR/'local-evidence.json'
if OUT.exists():raise FileExistsError(OUT)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rel(p):return p.relative_to(ROOT).as_posix()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
science=DIR/'stationary-checks.json'
assert sha(science)=='ec0c798d12dbb09f6f8020d939636f105637f98007b402284f9e5f6c4317f9f9'
data=read(science);assert data['status']=='passed'
for path,digest in data['input_sha256'].items():assert sha(ROOT/path)==digest,path
for path,digest in data['figure_sha256'].items():assert sha(ROOT/path)==digest,path
layout=read(DIR/'layout-evidence.json')
assert layout['source_sha256']==sha(science)
assert layout['script_sha256']==sha(ROOT/'code/plot_slater_stationary.py')
assert layout['figure_sha256']==sha(DIR/'finite-time-states-layout.png')
parentpath='results/susskind-slater-transfer/local-evidence.json'
assert sha(ROOT/parentpath)=='136226a7a720632fcd5bb92eac67a92724fbfd510a6dde70b701bd7faafac58f'
parent=read(ROOT/parentpath)
for path,digest in parent['hashes'].items():assert sha(ROOT/path)==digest,path
build=ROOT/f'tmp/pdfs/{STAGE}-build'
pdf=ROOT/f'output/pdf/de-sitter-{STAGE}.pdf'
log=DIR/'main.log'
assert sha(pdf)==sha(build/'main.pdf')
txt=log.read_text(errors='replace')
assert 'Output written on' in txt and not re.search(r'undefined|Overfull|Fatal error|Rerun to get',txt)
passes={}
for i in (1,2):
    p=build/f'pass{i}.txt';t=p.read_text(errors='replace')
    assert 'Output written on' in t and 'did not succeed' not in t
    passes[rel(p)]=sha(p)
info=subprocess.check_output(['pdfinfo',str(pdf)]).decode(errors='replace')
assert int(re.search(r'Pages:\s+(\d+)',info).group(1))==8
renders={}
for i in range(1,9):
    p=ROOT/f'tmp/pdfs/{STAGE}-render/page-{i}.png'
    assert p.read_bytes().startswith(b'\x89PNG\r\n\x1a\n')
    renders[rel(p)]=sha(p)
changed=[Path(x['Path']).name for x in read(build/'pre-final-render-hashes.json')
         if sha(Path(x['Path']))!=x['Hash'].lower()]
assert changed==['page-5.png','page-6.png','page-8.png']
files=[Path(__file__),ROOT/'code/check_slater_stationary.py',ROOT/'code/plot_slater_stationary.py',
       ROOT/'code/build_slater_stationary.ps1',ROOT/f'manuscript/{STAGE}/main.tex',
       ROOT/f'notes/{STAGE}.md',ROOT/f'literature/{STAGE}-sources.md',pdf,log,science,
       DIR/'layout-evidence.json',DIR/'finite-time-states.png',DIR/'finite-time-states-layout.png',
       DIR/'diagnostics/initial-build-failure.txt']
result=dict(status='Local research milestone; full de Sitter goal remains active',
            created_at_utc=datetime.now(timezone.utc).isoformat(),
            scientific_scope='Norm algebra and all stationary states in the existing h,M representation; no normal stationary density; no finite-beta KMS on full B(H); finite-time averaging; exact covariant fixed-sector Schur interface. No completed joint bilocal extension or finite-G dictionary.',
            scientific_output_hashes={rel(science):sha(science)},
            input_hashes_verified=data['input_sha256'],
            parameter_hash=hashlib.sha256(json.dumps(data['parameters'],sort_keys=True,separators=(',',':')).encode()).hexdigest(),
            previous_evidence={parentpath:dict(sha256=sha(ROOT/parentpath),verified_scientific_hashes=parent['hashes'])},
            hashes={rel(p):sha(p) for p in files},render_hashes=renders,build_pass_hashes=passes,
            checks=dict(q_cases=3,quadrature_rows=36,independent_time_integrals=3,exact_occupancy_systems=5,
                        max_quadrature_difference=max(x['max_L_difference'] for x in data['rows'])),
            build=dict(command='./code/build_slater_stationary.ps1',passes=2,final_exit_code_observed=0,pages=8,
                       review='All eight pages viewed; final pages 5,6,8 re-inspected after spacing/tick fixes; other five final render hashes unchanged.',
                       no_undefined_references_or_overfull=True,render='pdftoppm -r 90 -png'),
            environment=dict(recorder_python=platform.python_version(),science=data['environment']),
            diagnostics='Initial AA macro conflict retained; corrected macro name. Final layout rendered from sealed science JSON, no scientific rerun. Finite diagnostics are not interval certification or proof of infinite-dimensional classification.',
            git_commit_created=False,git_push_performed=False,external_reply_sent=False,
            recorder_scope='Read-only integrity and build evidence; does not certify analytic proofs.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',current_files=len(files),parent_files=len(parent['hashes']),pages=8,sha256=sha(OUT))))
