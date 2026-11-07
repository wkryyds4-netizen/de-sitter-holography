"""Read-only integrity recorder for the neutral-slope milestone."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,platform,re,subprocess

ROOT=Path(__file__).resolve().parent.parent
STAGE='susskind-slater-neutral-slope';DIR=ROOT/'results'/STAGE
OUT=DIR/'local-evidence.json'
if OUT.exists():raise FileExistsError(OUT)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rel(p):return p.relative_to(ROOT).as_posix()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
expected={
    'rank-checks.json':'7c8fa223ca99157f1f196f22187b025f028ef16e2c4263742badb4186d7c0123',
    'amplification-checks.json':'a168acdb46e530f879bc62d2251d935f49b773ce3d04eec7790eff5ecc139cb5',
    'dephasing-checks.json':'47109f31b526f3797f5a83d12758065acb76b6d281fab064bcdfc0745672013b'}
science={};inputs={}
for name,digest in expected.items():
    p=DIR/name;assert sha(p)==digest,name
    data=read(p);assert data['status']=='passed';science[name]=data
    for path,value in data['input_sha256'].items():
        assert sha(ROOT/path)==value,path
        inputs[path]=value
    for path,value in data.get('figure_sha256',{}).items():assert sha(ROOT/path)==value,path
parentpath='results/susskind-slater-covariant/local-evidence.json'
assert sha(ROOT/parentpath)=='0c26a1bffe0cd006a252d7e229851cb486c0022892564aed773e20da0a980aea'
parent=read(ROOT/parentpath)
for path,value in parent['hashes'].items():assert sha(ROOT/path)==value,path
pdf=ROOT/f'output/pdf/de-sitter-{STAGE}.pdf';log=DIR/'main.log'
build=ROOT/f'tmp/pdfs/{STAGE}-build'
assert sha(pdf)==sha(build/'main.pdf')
txt=log.read_text(errors='replace')
assert 'Output written on' in txt and not re.search(r'undefined|Overfull|Fatal error|Rerun to get',txt)
passes={}
for i in (1,2):
    p=build/f'pass{i}.txt';txt=p.read_text(errors='replace')
    assert 'Output written on' in txt and 'did not succeed' not in txt
    passes[rel(p)]=sha(p)
info=subprocess.check_output(['pdfinfo',str(pdf)]).decode(errors='replace')
assert int(re.search(r'Pages:\s+(\d+)',info).group(1))==8
renders={}
for i in range(1,9):
    p=ROOT/f'tmp/pdfs/{STAGE}-render/page-{i}.png'
    assert p.read_bytes().startswith(b'\x89PNG\r\n\x1a\n');renders[rel(p)]=sha(p)
files=[Path(__file__),ROOT/'code/check_slater_neutral_slope.py',
       ROOT/'code/check_slater_covariant_amplification.py',ROOT/'code/compute_slater_neutral_dephasing.py',
       ROOT/'code/build_slater_neutral_slope.ps1',ROOT/f'notes/{STAGE}.md',
       ROOT/f'literature/{STAGE}-sources.md',ROOT/f'manuscript/{STAGE}/main.tex',
       pdf,log,DIR/'neutral-slope-dephasing.png']+[DIR/name for name in expected]
result=dict(status='Local research milestone; full de Sitter goal remains active',
    created_at_utc=datetime.now(timezone.utc).isoformat(),
    scientific_scope='Microscopic unpinched neutral slope in all selected words/time; full ordinary q-Fock amplification and compact quotient; bounded densities for finite-grade selected vectors; finite-window ordered microscopic limit. Exact finite-N pinching is transferred only in purely covariant words, not yet general M-excited words.',
    scientific_output_hashes={rel(DIR/name):digest for name,digest in expected.items()},
    input_hashes_verified=inputs,
    parameters={name:data['parameters'] for name,data in science.items()},
    parameter_hashes={name:hashlib.sha256(json.dumps(data['parameters'],sort_keys=True,separators=(',',':')).encode()).hexdigest() for name,data in science.items()},
    previous_evidence={parentpath:dict(sha256=sha(ROOT/parentpath),verified_scientific_hashes=parent['hashes'])},
    hashes={rel(p):sha(p) for p in files},render_hashes=renders,build_pass_hashes=passes,
    checks=dict(exact_rank_systems=20,all_rank_DS_rows=12,
        amplification=science['amplification-checks.json']['counts'],dephasing_q_cases=3,
        dephasing_cutoff_time_rows=36,time_integral_difference=science['dephasing-checks.json']['independent_time_integral']['absolute_difference']),
    build=dict(command='./code/build_slater_neutral_slope.ps1',passes=2,final_exit_code_observed=0,pages=8,
        review='All eight final pages viewed at 90 dpi, including plotted figure; no layout revisions required.',
        no_undefined_references_or_overfull=True,render='pdftoppm -r 90 -png'),
    environment=dict(recorder_python=platform.python_version(),science={name:data['environment'] for name,data in science.items()}),
    diagnostics='Full rank tails and finite coefficient deviations retained. Rational algebra fixtures do not impose the physical curve. Jacobi cutoffs carry exact-arithmetic Taylor bounds; roundoff/infinite-product quadrature tails are not interval enclosed. No successful scientific driver rerun.',
    git_commit_created=False,git_push_performed=False,external_reply_sent=False,
    recorder_scope='Integrity/build evidence, not independent proof certification.')
with OUT.open('x',encoding='utf-8') as stream:json.dump(result,stream,indent=2);stream.write('\n')
print(json.dumps(dict(status='passed',current_files=len(files),parent_files=len(parent['hashes']),pages=8,sha256=sha(OUT))))
