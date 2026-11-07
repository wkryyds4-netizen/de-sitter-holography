"""Read saved results/builds and seal a fresh manifest without rerunning science."""
from datetime import datetime,timezone
from pathlib import Path
import hashlib,json,re,subprocess

ROOT=Path(__file__).resolve().parents[1]
STAGE='susskind-sector-memory-return'
OUT=ROOT/'results'/STAGE/'local-evidence.json'
if OUT.exists():raise FileExistsError(OUT)
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text(encoding='utf-8-sig'))

science=f'results/{STAGE}/interface-checks.json'
assert sha(science)=='05599623446e6af7dc9a9cc00f1ab175703db3ae53ce8bf6705c45c1b6adf3b6'
data=read(science)
assert data['status']=='passed'
assert data['counts']==dict(formal_degrees=21,finite_rows=9,nonzero_lower_chaos_rows=9,
                           dilation_times=15,schatten_bounds=90)
for p,h in data['input_sha256'].items():assert sha(p)==h,p

parents={};verified={};historical={}
mutable={'README.md','PLAN.md','NEXT.md','PROJECT.md','AGENTS.md'}
def verify(p,h,manifest):
    if p in mutable or p.startswith('journal/'):
        historical.setdefault(p,[]).append(dict(manifest=manifest,recorded_sha256=h))
        return
    assert sha(p)==h,(manifest,p)
    verified[p]=h
for stage in ['susskind-sector-neutral-transfer','susskind-sector-gradient-tail',
              'susskind-memory-return-expansion','susskind-memory-return-background',
              'susskind-slater-polynomial-density']:
    path=f'results/{stage}/local-evidence.json';parent=read(path)
    parents[path]=dict(sha256=sha(path))
    verified[path]=sha(path)
    for key in ['hashes','proof_inputs','input_hashes_verified','scientific_output_hashes',
                'verified_parent_scientific_hashes']:
        for p,h in parent.get(key,{}).items():
            verify(p,h,path)
            if p.endswith('.json') and p not in mutable and not p.startswith('journal/'):
                obj=read(p)
                if isinstance(obj,dict):
                    for name,digest in obj.get('input_sha256',{}).items():verify(name,digest,p)

proofs=[
    'notes/susskind-bilocal-fourth-tail-operator-draft.md',
    'manuscript/susskind-sector-gradient-tail/main.tex',
    'manuscript/susskind-sector-neutral-transfer/main.tex',
    'manuscript/susskind-memory-return-expansion/main.tex',
    'manuscript/susskind-memory-return-background/main.tex',
    'manuscript/susskind-slater-polynomial-density/main.tex',
]
pdf=f'output/pdf/de-sitter-{STAGE}.pdf'
build=f'tmp/pdfs/{STAGE}-build';render=f'tmp/pdfs/{STAGE}-render'
assert sha(pdf)==sha(f'{build}/main.pdf')
log=(ROOT/f'results/{STAGE}/main.log').read_text(errors='replace')
assert 'Output written on' in log
assert not re.search(r'undefined|Overfull|Fatal error|Rerun to get',log)
passes={}
for i in (1,2):
    p=f'{build}/pass{i}.txt';raw=(ROOT/p).read_bytes()
    text=raw.decode('utf-16') if raw.startswith((b'\xff\xfe',b'\xfe\xff')) else raw.decode()
    assert 'Output written on' in text and 'Fatal error' not in text
    passes[p]=sha(p)
info=subprocess.check_output(['pdfinfo',str(ROOT/pdf)]).decode(errors='replace')
assert int(re.search(r'Pages:\s+(\d+)',info).group(1))==7
renders={}
initial=read(f'tmp/pdfs/{STAGE}-initial/render-hashes.json')
changed=[]
for i in range(1,8):
    p=f'{render}/page-{i}.png'
    assert (ROOT/p).read_bytes().startswith(b'\x89PNG')
    renders[p]=sha(p)
    if sha(p)!=initial[f'page-{i}.png']:changed.append(i)
assert changed==[3]
files=[
    'code/sector_memory_return.py','code/check_sector_memory_return.py',
    'code/record_sector_memory_return.py','code/build_sector_memory_return.ps1',
    f'manuscript/{STAGE}/main.tex',f'notes/{STAGE}.md',
    f'literature/{STAGE}-sources.md',pdf,f'results/{STAGE}/main.log',science,
]
result=dict(status='Sectorwise memory-return suppression completed; full research goal active',
            created_at_utc=datetime.now(timezone.utc).isoformat(),
            scientific_scope='Exact force-defect inverse; compact-central all-finite-Schatten time bound; finite-parameter conserved reference; actual total-memory minus ordinary field O(p^-1/2); fixed polynomial preparation error.',
            scientific_output_hashes={science:sha(science)},input_hashes_verified=data['input_sha256'],
            checks=data['counts'],parameters=data['parameters'],environment=data['environment'],
            previous_evidence=parents,verified_parent_scientific_hashes=verified,
            historical_mutable_entrances=historical,proof_inputs={p:sha(p) for p in proofs},
            hashes={p:sha(p) for p in files},render_hashes=renders,build_pass_hashes=passes,
            build=dict(command='./code/build_sector_memory_return.ps1',passes=2,pages=7,
                       render='pdftoppm -r 90 -png',no_undefined_or_overfull=True,
                       visual_review='All seven pages viewed; after clarifying the common covariance lower bound, page 3 re-viewed and the other six pixel hashes unchanged.'),
            diagnostics='First scientific run passed and was sealed. Initial PDF/build retained before notation clarification; no successful scientific drivers or old recorders rerun.',
            remaining='Central four-u regression and charged finite-alpha interface; actual selected bilocal reconstruction, arbitrary nonpolynomial interlacings, spatial/mass dictionary and finite-G dynamics.',
            git_commit_created=False,git_push_performed=False,external_reply_sent=False,
            recorder_scope='Saved-byte integrity and build verification, not proof-assistant certification.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',current_files=len(files),parent_scientific_files=len(verified),
                     sha256=sha(OUT.relative_to(ROOT)))))
