"""Read-only stage evidence recorder. Does not execute scientific generators."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib,json,platform,re,subprocess

ROOT=Path(__file__).resolve().parent.parent
STAGE='susskind-slater-joint-fourth'
OUT=ROOT/f'results/{STAGE}/local-evidence.json'
EXPECTED={
 'exact-checks.json':'5021e7ce3cdf836a7d0b75025de546b0ce253beb02cf546287ccb6be215920f0',
 'gram-checks.json':'55d3aef7e7a051b28471360922604c4edf2adf5e2acb4039e4b08c28d12c0b6c',
 'covariance-checks.json':'8e4ab41f5a2e9db3286ea8c13a355d80a8b82600d142ac89240a7a16915d18dc'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rel(p):return p.relative_to(ROOT).as_posix()
def read(p):return json.loads(p.read_text(encoding='utf-8'))
if OUT.exists():raise FileExistsError(OUT)
science={};inputs={};versions={}
for name,digest in EXPECTED.items():
    p=ROOT/f'results/{STAGE}/{name}'
    assert sha(p)==digest,(p,'sealed hash')
    j=read(p);assert j['status']=='passed'
    for name2,digest2 in j['input_sha256'].items():
        assert sha(ROOT/name2)==digest2,(name2,'input hash')
        inputs[name2]=digest2
    science[rel(p)]=digest;versions[rel(p)]=j['environment']
parents={}
for name in [f'results/{x}/local-evidence.json' for x in ('susskind-slater-ladder','susskind-bilocal-four-point')]:
    j=read(ROOT/name);verified={};excluded=[]
    for path,digest in j['hashes'].items():
        if path in ('README.md','NEXT.md','PLAN.md','PROJECT.md','AGENTS.md') or path.startswith('journal/'):
            excluded.append(path);continue
        assert sha(ROOT/path)==digest,(path,'parent science')
        verified[path]=digest
    parents[name]=dict(sha256=sha(ROOT/name),verified_scientific_hashes=verified,excluded_mutable_or_private=excluded)
build=ROOT/f'tmp/pdfs/{STAGE}-build'
pdf=ROOT/f'output/pdf/de-sitter-{STAGE}.pdf'
log=ROOT/f'results/{STAGE}/main.log'
logtext=log.read_text(errors='replace')
assert 'Output written on' in logtext
assert not re.search(r'undefined|Overfull|Fatal error|Rerun to get',logtext)
assert sha(pdf)==sha(build/'main.pdf')
passes={}
for i in (1,2):
    p=build/f'pass{i}.txt';txt=p.read_text(errors='replace')
    assert 'Output written on' in txt and 'did not succeed' not in txt
    passes[rel(p)]=sha(p)
info=subprocess.run(['pdfinfo',str(pdf)],capture_output=True,check=True).stdout.decode(errors='replace')
pages=int(re.search(r'Pages:\s+(\d+)',info).group(1));assert pages==7
renders={}
for i in range(1,pages+1):
    p=ROOT/f'tmp/pdfs/{STAGE}-render/page-{i}.png'
    assert p.read_bytes().startswith(b'\x89PNG\r\n\x1a\n')
    renders[rel(p)]=sha(p)
files=[Path(__file__),ROOT/'code/check_slater_joint_fourth.py',ROOT/'code/check_slater_joint_gram.py',
 ROOT/'code/check_slater_joint_covariance.py',ROOT/'code/build_slater_joint_fourth.ps1',
 ROOT/f'manuscript/{STAGE}/main.tex',ROOT/f'notes/{STAGE}.md',
 ROOT/'literature/susskind-slater-joint-sources.md',pdf,log]
files.extend(ROOT/p for p in science)
hashes={rel(p):sha(p) for p in files}
diagnostics={rel(build/p):sha(build/p) for p in ('initial-sandbox-pass1.txt','initial-main.log')}
result=dict(status='Local research milestone; long-term de Sitter goal remains active',
 created_at_utc=datetime.now(timezone.utc).isoformat(),
 scope='Exact finite two-pair table, its fixed-degree limits, degree<=2 positive quotient and first Hamiltonian actions. Not an all-word joint representation or finite-time/finite-G solution.',
 parameters='Same original unit Hermitian complex Gaussian J, standard CAR, annealed half-filled Slater state, p^2/N -> fixed lambda>0; h normalized by exact sector variance; no joint self-averaging asserted.',
 scientific_output_hashes=science,input_hashes_verified=inputs,previous_evidence=parents,
 checks=dict(finite_CAR_label_assignments=2000499,fourth_words=324,pairing_formulas=972,
 second_words=36,covariance_columns=116,symbolic_identities=62,literature_scalar_comparisons=12),
 build=dict(command='./code/build_slater_joint_fourth.ps1',passes=2,final_script_exit_code_observed=0,
 pages=pages,no_undefined_references_or_overfull=True,render_command='pdftoppm -r 95 -png',
 visual_review='All seven final pages inspected by root using view_image; no clipping, overlap or broken references.',
 initial_attempt='Sandbox attempt produced a PDF but returned nonzero. Preserved diagnostics; final two-pass build succeeded with approved MiKTeX access.'),
 hashes=hashes,render_hashes=renders,build_pass_hashes=passes,diagnostic_hashes=diagnostics,
 environment=dict(recorder_python=platform.python_version(),scientific_environments=versions),
 publication=dict(this_stage_committed=False,this_stage_pushed=False,
 prior_authorized_commit='6e0f40af31f2ebf0ff0cf8f842244aaa04d82462',
 prior_authorized_issue_reply='https://github.com/wkryyds4-netizen/de-sitter-holography/issues/1#issuecomment-6708558378'),
 recorder_scope='File provenance, saved-output status, counts and build checks only; does not prove analytic claims or rerun science.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',scientific_outputs=len(science),current_files=len(hashes),
 parent_files=sum(len(x['verified_scientific_hashes']) for x in parents.values()),pages=pages,sha256=sha(OUT))))
