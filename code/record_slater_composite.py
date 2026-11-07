"""Read-only provenance for the same-J selected composite stage."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,re,subprocess,platform
ROOT=Path(__file__).resolve().parent.parent
STAGE='susskind-slater-composite';DIR=ROOT/'results'/STAGE;OUT=DIR/'local-evidence.json'
if OUT.exists():raise FileExistsError(OUT)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rel(p):return p.relative_to(ROOT).as_posix()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
science=DIR/'exact-checks.json'
assert sha(science)=='dcb7de18e4231bd129b23d82494909b47107103863c4af83012d3097b9eb603b'
data=read(science);assert data['status']=='passed'
for p,h in data['input_sha256'].items():assert sha(ROOT/p)==h,p
parents={};parentfiles={};historical_entrances={}
for stage in ['susskind-slater-polynomial-density','susskind-bilocal-four-point',
              'susskind-singlet-multimark','susskind-slater-joint-words']:
    p=f'results/{stage}/local-evidence.json';parent=read(ROOT/p)
    verified={}
    for path,digest in parent['hashes'].items():
        if path in {'README.md','PLAN.md','NEXT.md'}:
            historical_entrances.setdefault(p,{})[path]=dict(historical_sha256=digest,
                reason='Mutable project entrance; historical snapshot retained, not compared to current bytes.')
            continue
        assert sha(ROOT/path)==digest,path;parentfiles[path]=digest;verified[path]=digest
    parents[p]=dict(sha256=sha(ROOT/p),verified_scientific_hashes=verified)
assert parents['results/susskind-slater-polynomial-density/local-evidence.json']['sha256']=='484bf963aa1b457a187e3e05a1cb32d23a5006f81c7eb9636e5fb56bb89da075'
pdf=ROOT/f'output/pdf/de-sitter-{STAGE}.pdf';log=DIR/'main.log'
build=ROOT/f'tmp/pdfs/{STAGE}-build';render=ROOT/f'tmp/pdfs/{STAGE}-render'
assert sha(pdf)==sha(build/'main.pdf')
t=log.read_text(errors='replace')
assert 'Output written on' in t and not re.search(r'undefined|Overfull|Fatal error|Rerun to get',t)
passes={}
for i in (1,2):
    p=build/f'pass{i}.txt';t=p.read_text(errors='replace')
    assert 'Output written on' in t and 'did not succeed' not in t;passes[rel(p)]=sha(p)
info=subprocess.check_output(['pdfinfo',str(pdf)]).decode(errors='replace')
assert int(re.search(r'Pages:\s+(\d+)',info).group(1))==8
renders={}
for i in range(1,9):
    p=render/f'page-{i}.png';assert p.read_bytes().startswith(b'\x89PNG\r\n\x1a\n');renders[rel(p)]=sha(p)
files=[Path(__file__),ROOT/'code/check_slater_composite_sources.py',ROOT/'code/slater_composite_sources.py',
       ROOT/'code/build_slater_composite.ps1',ROOT/f'notes/{STAGE}.md',ROOT/f'literature/{STAGE}-sources.md',
       ROOT/f'manuscript/{STAGE}/main.tex',pdf,log,science,DIR/'diagnostics/initial-recorder-failure.txt']
result=dict(status='Local research milestone; full de Sitter goal remains active',
    created_at_utc=datetime.now(timezone.utc).isoformat(),
    scientific_scope='Constructive all-fixed-word same-J H/W/J/M moments; direct central-sector mixed-q identification; exact curvature and neutral-slope pinching within covariant blocks between arbitrary fixed polynomial preparations; bounded normal covariant representation on the enlarged GNS. Closed selected-composite transfer, interlaced nonlinear pinching and the bulk dictionary remain open.',
    scientific_output_hashes={rel(science):sha(science)},input_hashes_verified=data['input_sha256'],
    parameters=data['parameters'],parameter_hash=hashlib.sha256(json.dumps(data['parameters'],sort_keys=True,separators=(',',':')).encode()).hexdigest(),
    previous_evidence=parents,historical_mutable_entrances=historical_entrances,
    hashes={rel(p):sha(p) for p in files},render_hashes=renders,build_pass_hashes=passes,
    checks=dict(finite_CAR=sum(x['comparisons'] for x in data['finite_CAR']),negative=sum(x['negative'] for x in data['finite_CAR']),
                all_G_sources=data['central_source_comparisons'],central_words=data['central_words'],
                old_selected_words=data['old_selected_words'],symbolic_moments=len(data['symbolic_selected_moments']),
                rational_Gram_cases=len(data['Gram_checks']),Gram_size=len(data['Gram_preparations'])),
    build=dict(command='./code/build_slater_composite.ps1',passes=2,final_exit_code_observed=0,pages=8,
        review='All eight final pages viewed at 90 dpi; no clipping, missing references or formula overflow observed.',
        no_undefined_references_or_overfull=True,render='pdftoppm -r 90 -png'),
    environment=dict(recorder_python=platform.python_version(),science=data['environment']),
    diagnostics='Fresh scientific driver passed on first run and was not rerun. Exploratory pure-helper evaluations preceded it without scientific output writes. First two-pass LaTeX build passed. Initial recorder stopped on historical mutable NEXT.md; audit found only README/PLAN/NEXT changed in two older manifests. These historical entrances are explicitly excluded from current-byte verification; all their scientific files are verified.',
    git_commit_created=False,git_push_performed=False,external_reply_sent=False,
    recorder_scope='File/build integrity, not independent certification of arbitrary-order analytic proofs.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',current_files=len(files),unique_parent_files=len(parentfiles),pages=8,sha256=sha(OUT))))
