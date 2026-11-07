"""Read-only proof/check/build integrity record; no scientific rerun."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,platform,re,subprocess
ROOT=Path(__file__).resolve().parent.parent
STAGE='susskind-slater-orbit-density';DIR=ROOT/'results'/STAGE
OUT=DIR/'local-evidence.json'
if OUT.exists():raise FileExistsError(OUT)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rel(p):return p.relative_to(ROOT).as_posix()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
science=DIR/'exact-checks.json'
assert sha(science)=='95313475713ec714a8d525d8a0d32315afb92ed7ee41983bcd96719a4c4120ac'
data=read(science);assert data['status']=='passed'
for path,value in data['input_sha256'].items():assert sha(ROOT/path)==value,path
parentpath='results/susskind-slater-neutral-slope/local-evidence.json'
assert sha(ROOT/parentpath)=='cf278d649e5c6d76bebe925ab79ff894b6555a806538e066096adbfadefafc82'
parent=read(ROOT/parentpath)
for path,value in parent['hashes'].items():assert sha(ROOT/path)==value,path
pdf=ROOT/f'output/pdf/de-sitter-{STAGE}.pdf';log=DIR/'main.log';build=ROOT/f'tmp/pdfs/{STAGE}-build'
assert sha(pdf)==sha(build/'main.pdf')
txt=log.read_text(errors='replace')
assert 'Output written on' in txt and not re.search(r'undefined|Overfull|Fatal error|Rerun to get',txt)
passes={}
for i in (1,2):
    p=build/f'pass{i}.txt';txt=p.read_text(errors='replace')
    assert 'Output written on' in txt and 'did not succeed' not in txt
    passes[rel(p)]=sha(p)
info=subprocess.check_output(['pdfinfo',str(pdf)]).decode(errors='replace')
assert int(re.search(r'Pages:\s+(\d+)',info).group(1))==6
renders={}
for i in range(1,7):
    p=ROOT/f'tmp/pdfs/{STAGE}-render/page-{i}.png'
    assert p.read_bytes().startswith(b'\x89PNG\r\n\x1a\n');renders[rel(p)]=sha(p)
changed=[Path(x['Path']).name for x in read(build/'pre-final-render-hashes.json')
         if sha(Path(x['Path']))!=x['Hash'].lower()]
assert changed==['page-1.png'],changed
files=[Path(__file__),ROOT/'code/check_slater_orbit_density.py',ROOT/'code/build_slater_orbit_density.ps1',
       ROOT/f'notes/{STAGE}.md',ROOT/f'literature/{STAGE}-sources.md',ROOT/f'manuscript/{STAGE}/main.tex',
       pdf,log,science,DIR/'diagnostics/initial-check-failure.txt',DIR/'diagnostics/initial-build-failure.txt']
result=dict(status='Local research milestone; full de Sitter goal remains active',
    created_at_utc=datetime.now(timezone.utc).isoformat(),
    scientific_scope='Finite-N uniform orbit-density moments for homogeneous linear stabilizer-covariant preparations, including M-dagger. Exact neutral-slope pinching transferred in the one-selected prepared state and covariant post-dressings. Not yet arbitrary nonlinear/interlaced selected words or the full bilocal dictionary.',
    scientific_output_hashes={rel(science):sha(science)},input_hashes_verified=data['input_sha256'],
    parameters=data['parameters'],parameter_hash=hashlib.sha256(json.dumps(data['parameters'],sort_keys=True,separators=(',',':')).encode()).hexdigest(),
    previous_evidence={parentpath:dict(sha256=sha(ROOT/parentpath),verified_scientific_hashes=parent['hashes'])},
    hashes={rel(p):sha(p) for p in files},render_hashes=renders,build_pass_hashes=passes,
    checks=data['counts'],
    build=dict(command='./code/build_slater_orbit_density.ps1',passes=2,final_exit_code_observed=0,pages=6,
        review='All six pages viewed at 90 dpi; final page 1 re-inspected after explicit L2 integrability wording. Other five render hashes unchanged.',
        no_undefined_references_or_overfull=True,render='pdftoppm -r 90 -png'),
    environment=dict(recorder_python=platform.python_version(),science=data['environment']),
    diagnostics='Initial exact comparison required canonicalization of complex expressions; no mathematical identity changed. Initial LaTeX mathscr macro unavailable, replaced by mathcal. Failed diagnostics retained. Successful science output not regenerated.',
    git_commit_created=False,git_push_performed=False,external_reply_sent=False,
    recorder_scope='Integrity/build evidence; finite diagnostics do not independently certify the general analytic theorem.')
with OUT.open('x',encoding='utf-8') as stream:json.dump(result,stream,indent=2);stream.write('\n')
print(json.dumps(dict(status='passed',current_files=len(files),parent_files=len(parent['hashes']),pages=6,sha256=sha(OUT))))
