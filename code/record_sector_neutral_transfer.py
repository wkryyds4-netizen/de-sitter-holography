"""Read-only verification of saved science and proof inputs, then fresh manifest."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import re
import subprocess

ROOT=Path(__file__).resolve().parents[1]
STAGE='susskind-sector-neutral-transfer'
DIR=ROOT/'results'/STAGE
OUT=DIR/'local-evidence.json'
if OUT.exists():
    raise FileExistsError(OUT)


def sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def read(path):
    return json.loads((ROOT/path).read_text(encoding='utf-8-sig'))


science=f'results/{STAGE}/interface-checks.json'
assert sha(science)=='6f10a97e7b014b0e267dbdcb2c5686bd758849d9d535bf34ca840050edc69635'
data=read(science)
assert data['status']=='passed'
assert data['counts']==dict(finite_rows=15,nonzero_filter_rows=6,
                            exact_comparisons=56,selected_rows=18)
for p,h in data['input_sha256'].items():
    assert sha(p)==h,p
parent_stages=[
    'susskind-sector-covariance-tangent','susskind-sector-gradient-tail',
    'susskind-memory-wick-rank','susskind-memory-central-window',
    'susskind-memory-neutral-marked-series','susskind-memory-wick-descendants',
    'susskind-slater-higher-wick','susskind-slater-polynomial-density',
]
parents={}
verified={}
historical={}
mutable={'README.md','PLAN.md','NEXT.md','PROJECT.md'}
for stage in parent_stages:
    path=f'results/{stage}/local-evidence.json'
    parent=read(path)
    parents[path]=dict(sha256=sha(path),checked_scientific_files=[])
    for key in ('hashes','proof_inputs','verified_parent_scientific_hashes'):
        for p,h in parent.get(key,{}).items():
            if p in mutable or p.startswith('journal/'):
                historical.setdefault(p,[]).append(dict(manifest=path,recorded_sha256=h))
                continue
            assert sha(p)==h,(path,p)
            verified[p]=h
            parents[path]['checked_scientific_files'].append(p)
            if p.endswith('.json'):
                obj=read(p)
                if isinstance(obj,dict):
                    for name,digest in obj.get('input_sha256',{}).items():
                        assert sha(name)==digest,(p,name)
                        verified[name]=digest
    verified[path]=sha(path)
proofs=[
    'manuscript/susskind-memory-wick-rank/main.tex',
    'manuscript/susskind-memory-central-window/main.tex',
    'manuscript/susskind-memory-neutral-marked-series/main.tex',
    'manuscript/susskind-memory-wick-descendants/main.tex',
    'manuscript/susskind-slater-higher-wick/main.tex',
    'manuscript/susskind-slater-polynomial-density/main.tex',
    'manuscript/susskind-sector-covariance-tangent/main.tex',
    'manuscript/susskind-sector-gradient-tail/main.tex',
    'notes/susskind-bilocal-fourth-tail-operator-draft.md',
]
pdf=f'output/pdf/de-sitter-{STAGE}.pdf'
build=f'tmp/pdfs/{STAGE}-build'
render=f'tmp/pdfs/{STAGE}-render'
assert sha(pdf)==sha(f'{build}/main.pdf')
log=(DIR/'main.log').read_text(errors='replace')
assert 'Output written on' in log
assert not re.search(r'undefined|Overfull|Fatal error|Rerun to get',log)
passes={}
for i in (1,2):
    p=f'{build}/pass{i}.txt'
    raw=(ROOT/p).read_bytes()
    text=raw.decode('utf-16') if raw.startswith((b'\xff\xfe',b'\xfe\xff')) else raw.decode()
    assert 'Output written on' in text and 'Fatal error' not in text
    passes[p]=sha(p)
info=subprocess.check_output(['pdfinfo',str(ROOT/pdf)]).decode(errors='replace')
assert int(re.search(r'Pages:\s+(\d+)',info).group(1))==8
renders={}
for i in range(1,9):
    p=f'{render}/page-{i}.png'
    assert (ROOT/p).read_bytes().startswith(b'\x89PNG')
    renders[p]=sha(p)
files=[
    'code/sector_neutral_transfer.py','code/check_sector_neutral_transfer.py',
    'code/record_sector_neutral_transfer.py','code/build_sector_neutral_transfer.ps1',
    f'manuscript/{STAGE}/main.tex',f'notes/{STAGE}.md',
    f'literature/{STAGE}-sources.md',pdf,f'results/{STAGE}/main.log',science,
]
result=dict(
    status='Neutral component completed in stated sector/preparation scope; full research goal active',
    created_at_utc=datetime.now(timezone.utc).isoformat(),
    scientific_scope='Fixed-degree composite covariance filter via graph interpolation; same-copy compact-central neutral finite-time strong transfer in every fixed finite Schatten order; half-filled central exact pinching and controlled selected moments.',
    scientific_output_hashes={science:sha(science)},input_hashes_verified=data['input_sha256'],
    checks=data['counts'],parameters=data['parameters'],environment=data['environment'],
    previous_evidence=parents,verified_parent_scientific_hashes=verified,
    historical_mutable_entrances=historical,
    proof_inputs={p:sha(p) for p in proofs},hashes={p:sha(p) for p in files},
    render_hashes=renders,build_pass_hashes=passes,
    build=dict(command='./code/build_sector_neutral_transfer.ps1',passes=2,pages=8,
               render='pdftoppm -r 90 -png',no_undefined_or_overfull=True,
               visual_review='All eight pages viewed. Top-chaos degree hypothesis clarified; after final rebuild page 3 rechecked and the other seven render hashes equal the previously viewed pages.'),
    diagnostics='Scientific checker passed its first run; no previous scientific drivers or recorders rerun. Initial clean PDF retained locally before clarifying the top-chaos degree hypothesis.',
    remaining='Sectorwise total-memory/collective-return reduction and bilocal reconstruction; arbitrary interlaced nonpolynomial selected blocks; spatial/mass dictionary and finite-G dynamics.',
    git_commit_created=False,git_push_performed=False,external_reply_sent=False,
    recorder_scope='Saved-byte and build validation; not proof-assistant certification of analytic estimates.')
with OUT.open('x',encoding='utf-8') as f:
    json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({'status':'passed','current_files':len(files),
                  'parent_scientific_files':len(verified),'sha256':sha(OUT.relative_to(ROOT))}))
