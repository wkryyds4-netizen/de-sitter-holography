"""Read-only verification followed by one new provenance manifest."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,platform,re,subprocess
ROOT=Path(__file__).resolve().parent.parent
STAGE='susskind-slater-composite-transfer'
DIR=ROOT/'results'/STAGE;OUT=DIR/'local-evidence.json'
if OUT.exists():raise FileExistsError(OUT)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def rel(p):return p.relative_to(ROOT).as_posix()
science=DIR/'exact-checks.json'
assert sha(science)=='f13bd76194cb358ae1f471db1a4b567ab77854040b8b78619db7e8a1249487a2'
data=read(science);assert data['status']=='passed'
for p,h in data['input_sha256'].items():assert sha(ROOT/p)==h,p
parentpath=ROOT/'results/susskind-slater-composite/local-evidence.json'
assert sha(parentpath)=='c5f9fbdeb165a7e7373e8abedd3310873f2a055765643e884229ce20a1a9bad8'
parent=read(parentpath);verified={}
for p,h in parent['hashes'].items():assert sha(ROOT/p)==h,p;verified[p]=h
for reference,entry in parent['previous_evidence'].items():
    assert sha(ROOT/reference)==entry['sha256'],reference
    verified[reference]=entry['sha256']
    for p,h in entry['verified_scientific_hashes'].items():assert sha(ROOT/p)==h,p;verified[p]=h
pdf=ROOT/f'output/pdf/de-sitter-{STAGE}.pdf'
build=ROOT/f'tmp/pdfs/{STAGE}-build';render=ROOT/f'tmp/pdfs/{STAGE}-render'
assert sha(pdf)==sha(build/'main.pdf')
log=DIR/'main.log';text=log.read_text(errors='replace')
assert 'Output written on' in text and not re.search(r'undefined|Overfull|Fatal error|Rerun to get',text)
passes={}
for i in (1,2):
    p=build/f'pass{i}.txt';text=p.read_text(errors='replace')
    assert 'Output written on' in text and 'did not succeed' not in text
    passes[rel(p)]=sha(p)
info=subprocess.check_output(['pdfinfo',str(pdf)]).decode(errors='replace')
assert int(re.search(r'Pages:\s+(\d+)',info).group(1))==6
renders={}
for i in range(1,7):
    p=render/f'page-{i}.png';assert p.read_bytes().startswith(b'\x89PNG\r\n\x1a\n');renders[rel(p)]=sha(p)
files=[Path(__file__),ROOT/'code/slater_composite_transfer.py',ROOT/'code/check_slater_composite_transfer.py',
       ROOT/'code/build_slater_composite_transfer.ps1',ROOT/f'notes/{STAGE}.md',
       ROOT/f'literature/{STAGE}-sources.md',ROOT/f'manuscript/{STAGE}/main.tex',
       pdf,log,science,DIR/'diagnostics/initial-reference-stop.txt']
result=dict(status='Local research milestone; long-term de Sitter goal remains active',
    created_at_utc=datetime.now(timezone.utc).isoformat(),
    scientific_scope='All-fixed-word geometric source classification and explicit open-endpoint/link recurrence; positive physical form on the reachable module, with GNS adjoints. No minimal orthogonal basis, general polynomial complexity, enlarged selected algebra classification, or bulk spatial/mass dictionary is asserted.',
    scientific_output_hashes={rel(science):sha(science)},input_hashes_verified=data['input_sha256'],
    parameters=data['parameters'],parameter_hash=hashlib.sha256(json.dumps(data['parameters'],sort_keys=True,separators=(',',':')).encode()).hexdigest(),
    parent_manifest=dict(path=rel(parentpath),sha256=sha(parentpath)),verified_parent_scientific_hashes=verified,
    historical_mutable_entrances=parent['historical_mutable_entrances'],
    hashes={rel(p):sha(p) for p in files},render_hashes=renders,build_pass_hashes=passes,
    checks={k:data[k] for k in ['source_exponents','maximum_nonzero_source_order','source_moment_comparisons','old_six_letter_words','long_word_checks','Gram_checks','unpruned_checks']},
    build=dict(command='./code/build_slater_composite_transfer.ps1',passes=2,pages=6,
        no_undefined_references_or_overfull=True,render='pdftoppm -r 90 -png',
        visual_review='All six final pages viewed at 90 dpi; no clipping, formula overflow, missing citations or layout defects observed.'),
    environment=dict(recorder_python=platform.python_version(),science=data['environment']),
    diagnostics='First checker interrupted during an impractical H^24 reference enumeration, before any science JSON existed; independent reference changed to the verified Touchard--Riordan sum. Next run completed successfully and was not rerun. First two-pass PDF build succeeded. Exploratory single-preparation j^n comparisons are recorded only as a next-task clue, not certified general results.',
    git_commit_created=False,git_push_performed=False,external_reply_sent=False,
    recorder_scope='Read-only byte/build verification; finite tests are not independent certification of arbitrary-order analytic theorems.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',current_files=len(files),parent_scientific_files=len(verified),pages=6,sha256=sha(OUT))))
