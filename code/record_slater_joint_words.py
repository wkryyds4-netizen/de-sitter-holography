"""Read-only provenance and final PDF record; no scientific main is executed."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,platform,re,subprocess
ROOT=Path(__file__).resolve().parent.parent
STAGE='susskind-slater-joint-words'
OUT=ROOT/f'results/{STAGE}/local-evidence.json'
EXPECTED={
 'site-checks.json':'864f4342efb51562d8a877ed2a1f9be14d7d1fa9d9d4e19d01a3d84a200c6f80',
 'word-checks.json':'f4413f05fe251d1ff1cffc927f83af5b013efba32a7bc23c2d4320e93ce0e3bd',
 'time-bound-checks.json':'1e244b4d232e55008a0243cc4a4ecb2afc134b26364d12b6d034d285c4469f34'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rel(p):return p.relative_to(ROOT).as_posix()
def read(p):return json.loads(p.read_text(encoding='utf-8'))
if OUT.exists():raise FileExistsError(OUT)
science={};inputs={};versions={};parameter_hashes={};data={}
for name,digest in EXPECTED.items():
    p=ROOT/f'results/{STAGE}/{name}';assert sha(p)==digest
    j=read(p);data[name]=j;assert j['status']=='passed'
    for name2,digest2 in j['input_sha256'].items():
        assert sha(ROOT/name2)==digest2,(name2,'input hash');inputs[name2]=digest2
    science[rel(p)]=digest;versions[rel(p)]=j['environment']
    parameter_hashes[rel(p)]=hashlib.sha256(json.dumps(j['parameters'],sort_keys=True,separators=(',',':')).encode()).hexdigest()
site=data['site-checks.json']
assert not site['failures'] and site['quadratic_checks']==1814
assert sum(x['checks']['CAR_vs_site_assignments'] for x in site['fixtures'])==1461120
assert sum(x['checks']['polynomial_coefficients'] for x in site['fixtures'])==1215
assert sum(x['signs'].get('-1',0) for x in site['fixtures'])==80
assert len(data['word-checks.json']['sixth_words'])==729
assert data['time-bound-checks.json']['covariance_scalar_checks']==348
parent_name='results/susskind-slater-joint-fourth/local-evidence.json'
parent=read(ROOT/parent_name);parent_hashes={}
for path,digest in parent['hashes'].items():
    assert sha(ROOT/path)==digest,(path,'parent hash');parent_hashes[path]=digest
build=ROOT/f'tmp/pdfs/{STAGE}-build';pdf=ROOT/f'output/pdf/de-sitter-{STAGE}.pdf'
log=ROOT/f'results/{STAGE}/main.log';txt=log.read_text(errors='replace')
assert 'Output written on' in txt and not re.search(r'undefined|Overfull|Fatal error|Rerun to get',txt)
assert sha(pdf)==sha(build/'main.pdf')
passes={}
for i in (1,2):
    p=build/f'pass{i}.txt';txt=p.read_text(errors='replace')
    assert 'Output written on' in txt and 'did not succeed' not in txt
    passes[rel(p)]=sha(p)
info=subprocess.run(['pdfinfo',str(pdf)],capture_output=True,check=True).stdout.decode(errors='replace')
pages=int(re.search(r'Pages:\s+(\d+)',info).group(1));assert pages==9
renders={}
for i in range(1,pages+1):
    p=ROOT/f'tmp/pdfs/{STAGE}-render/page-{i}.png'
    assert p.read_bytes().startswith(b'\x89PNG\r\n\x1a\n');renders[rel(p)]=sha(p)
files=[Path(__file__),ROOT/'code/slater_joint_site_polynomial.py',ROOT/'code/check_slater_joint_sites.py',
 ROOT/'code/check_slater_joint_word_law.py',ROOT/'code/check_slater_joint_time_bound.py',
 ROOT/'code/build_slater_joint_words.ps1',ROOT/f'manuscript/{STAGE}/main.tex',ROOT/f'notes/{STAGE}.md',
 ROOT/'literature/susskind-slater-joint-word-sources.md',pdf,log]+[ROOT/p for p in science]
result=dict(status='Local milestone; full de Sitter objective remains active',created_at_utc=datetime.now(timezone.utc).isoformat(),
 scientific_scope='All fixed same-J Slater word limits, their positive GNS representation, essential self-adjointness of h, and actual annealed finite-time polynomial-insertion correlators on compact complex time sets. Not quenched concentration, a stationary observer dictionary, or finite-G scattering.',
 parameters='Half filling, standard CAR, overlapping Hermitian complex Gaussian p-body model, fixed Slater state, fixed lambda>0, h normalized by exact sector variance.',
 scientific_output_hashes=science,input_hashes_verified=inputs,parameter_hashes=parameter_hashes,
 previous_evidence={parent_name:dict(sha256=sha(ROOT/parent_name),verified_scientific_hashes=parent_hashes)},
 checks=dict(CAR_label_assignments=1461120,coefficient_extractions=1215,negative_CAR_assignments=80,
 quadratic_graphs=1814,sixth_words=729,rational_Gram_cases=4,covariance_counts=348,tail_diagnostics=9),
 hashes={rel(p):sha(p) for p in files},render_hashes=renders,build_pass_hashes=passes,
 build=dict(command='./code/build_slater_joint_words.ps1',passes=2,final_script_exit_code_observed=0,pages=9,
 no_undefined_references_or_overfull=True,render='pdftoppm -r 90 -png',
 review='All nine final pages viewed by root; no clipping, overlap or broken references.'),
 environment=dict(recorder_python=platform.python_version(),scientific_environments=versions),
 git_commit_created=False,git_push_performed=False,external_reply_sent=False,
 recorder_scope='Read-only hash/status/coverage/build verification. Does not prove analytic statements or rerun scientific generators.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',scientific_outputs=len(science),current_files=len(files),
 parent_files=len(parent_hashes),pages=pages,sha256=sha(OUT))))
