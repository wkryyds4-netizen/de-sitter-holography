"""Read-only provenance for the polynomial-density stage; no scientific reruns."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,re,subprocess,platform
ROOT=Path(__file__).resolve().parent.parent
STAGE='susskind-slater-polynomial-density';DIR=ROOT/'results'/STAGE
OUT=DIR/'local-evidence.json'
if OUT.exists():raise FileExistsError(OUT)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rel(p):return p.relative_to(ROOT).as_posix()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
science=DIR/'exact-checks.json'
assert sha(science)=='644911ef2ffb91c3b6b5d4419fe757653d0ba6f1decc09053888983c6b2e4951'
data=read(science);assert data['status']=='passed'
for p,h in data['input_sha256'].items():assert sha(ROOT/p)==h,p
parentpath='results/susskind-slater-orbit-density/local-evidence.json'
assert sha(ROOT/parentpath)=='0864d7bca7abbc5bdfb86a23e8e21b08bb97c7e852dc3b3a28a752b6a5bea028'
parent=read(ROOT/parentpath)
for p,h in parent['hashes'].items():assert sha(ROOT/p)==h,p
pdf=ROOT/f'output/pdf/de-sitter-{STAGE}.pdf';log=DIR/'main.log'
build=ROOT/f'tmp/pdfs/{STAGE}-build';render=ROOT/f'tmp/pdfs/{STAGE}-render'
assert sha(pdf)==sha(build/'main.pdf')
t=log.read_text(errors='replace')
assert 'Output written on' in t and not re.search(r'undefined|Overfull|Fatal error|Rerun to get',t)
passes={}
for i in (1,2):
    p=build/f'pass{i}.txt';t=p.read_text(errors='replace')
    assert 'Output written on' in t and 'did not succeed' not in t
    passes[rel(p)]=sha(p)
info=subprocess.check_output(['pdfinfo',str(pdf)]).decode(errors='replace')
assert int(re.search(r'Pages:\s+(\d+)',info).group(1))==8
renders={}
for i in range(1,9):
    p=render/f'page-{i}.png';assert p.read_bytes().startswith(b'\x89PNG\r\n\x1a\n')
    renders[rel(p)]=sha(p)
old=read(render/'pre-final-hashes.json')
assert [name for name,h in old.items() if sha(render/name)!=h]==['page-6.png','page-7.png']
files=[Path(__file__),ROOT/'code/check_slater_polynomial_density.py',ROOT/'code/build_slater_polynomial_density.ps1',
       ROOT/f'notes/{STAGE}.md',ROOT/f'literature/{STAGE}-sources.md',ROOT/f'manuscript/{STAGE}/main.tex',pdf,log,science]
result=dict(status='Local research milestone; full de Sitter goal remains active',
    created_at_utc=datetime.now(timezone.utc).isoformat(),
    scientific_scope='Uniform orbit-density moments for fixed Gaussian-degree stabilizer-covariant preparations in half-filled double scaling. Exact neutral-slope pinching on the selected polynomial core and within a covariant block between polynomial preparations. Two normalized quadratic response limits. Not all interlaced nonpolynomial pinched moments or the full bulk dictionary.',
    scientific_output_hashes={rel(science):sha(science)},input_hashes_verified=data['input_sha256'],
    parameters=data['parameters'],parameter_hash=hashlib.sha256(json.dumps(data['parameters'],sort_keys=True,separators=(',',':')).encode()).hexdigest(),
    previous_evidence={parentpath:dict(sha256=sha(ROOT/parentpath),verified_scientific_hashes=parent['hashes'])},
    hashes={rel(p):sha(p) for p in files},render_hashes=renders,build_pass_hashes=passes,
    checks=dict(seed_weights=data['seed_weights'],partial_weights=data['partial_weights'],
                factorial_rational=data['factorial_rational_checks'],near_weights=len(data['near_weight_rows']),
                conditional_Wick=data['conditional_Wick_checks'],rectangular_moments=len(data['rectangular_tensors'])),
    build=dict(command='./code/build_slater_polynomial_density.ps1',passes=2,final_exit_code_observed=0,pages=8,
        review='All eight pages viewed at 90 dpi; final pages 6 and 7 re-inspected after clarifying fixed polynomial coefficients and occupied index notation. Other six page hashes unchanged.',
        no_undefined_references_or_overfull=True,render='pdftoppm -r 90 -png'),
    environment=dict(recorder_python=platform.python_version(),science=data['environment']),
    diagnostics='Scientific driver passed on first run and was not rerun. Both builds passed; final changes were notation only. pdftotext emitted a MiKTeX logger permission warning during ancillary text extraction; PDF rendering and visual review succeeded.',
    git_commit_created=False,git_push_performed=False,external_reply_sent=False,
    recorder_scope='File/build integrity. Finite algebra checks do not independently certify general analytic proofs.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',current_files=len(files),parent_files=len(parent['hashes']),pages=8,sha256=sha(OUT))))
