"""Seal saved tangent evidence; do not execute any scientific driver."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
STAGE = 'susskind-sector-covariance-tangent'
OUT = ROOT/'results'/STAGE/'local-evidence.json'
if OUT.exists():
    raise FileExistsError(OUT)


def sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


science_path = 'results/susskind-sector-neutral-transfer/tangent-checks.json'
science_hash = '2bcc489f4dd8a5d81d0c0a77f6de4a16af81997be506247dc282113e30989e1c'
assert sha(science_path) == science_hash
data = json.loads((ROOT/science_path).read_text())
assert data['status'] == 'passed'
assert data['counts']['path_equalities'] == 1824
assert data['counts']['norm_bounds'] == 1824
assert data['counts']['zero_reference_cases'] == 608
assert len(data['scalar_rows']) == 60
for p, h in data['input_sha256'].items():
    assert sha(p) == h, p
parent_path = 'results/susskind-sector-gradient-tail/local-evidence.json'
parent_hash = '6a4543be7dfb804201ceeee2e2e27b96ea908c1bd6c3f5f26aee6c1a17d7c15e'
assert sha(parent_path) == parent_hash
parent = json.loads((ROOT/parent_path).read_text())
inherited = {}
for key in ('hashes', 'proof_inputs', 'verified_parent_scientific_hashes'):
    for p, h in parent[key].items():
        assert sha(p) == h, p
        inherited[p] = h
proofs = [
    'notes/susskind-bilocal-fourth-tail-diagram-draft.md',
    'notes/susskind-bilocal-fourth-tail-operator-draft.md',
]
pdf = f'output/pdf/de-sitter-{STAGE}.pdf'
build = f'tmp/pdfs/{STAGE}-build'
render = f'tmp/pdfs/{STAGE}-render'
assert sha(pdf) == sha(f'{build}/main.pdf')
log = (ROOT/f'results/{STAGE}/main.log').read_text(errors='replace')
assert 'Output written on' in log
assert not re.search(r'undefined|Overfull|Fatal error|Rerun to get', log)
passes = {}
for i in (1, 2):
    p = f'{build}/pass{i}.txt'
    raw = (ROOT/p).read_bytes()
    txt = raw.decode('utf-16') if raw.startswith((b'\xff\xfe', b'\xfe\xff')) else raw.decode()
    assert 'Output written on' in txt and 'Fatal error' not in txt
    passes[p] = sha(p)
info = subprocess.check_output(['pdfinfo', str(ROOT/pdf)]).decode(errors='replace')
assert int(re.search(r'Pages:\s+(\d+)', info).group(1)) == 4
renders = {}
for i in range(1, 5):
    p = f'{render}/page-{i}.png'
    assert (ROOT/p).read_bytes().startswith(b'\x89PNG')
    renders[p] = sha(p)
files = [
    'code/record_sector_covariance_tangent.py',
    'code/build_sector_covariance_tangent.ps1',
    f'manuscript/{STAGE}/main.tex', f'notes/{STAGE}.md',
    f'literature/{STAGE}-sources.md', f'results/{STAGE}/main.log', pdf,
    science_path, *data['input_sha256'],
]
files += [
    p.relative_to(ROOT).as_posix()
    for p in sorted((ROOT/'results/susskind-sector-neutral-transfer/diagnostics').iterdir())
    if p.is_file()
]
result = {
    'created_at_utc': datetime.now(timezone.utc).isoformat(),
    'status': 'Covariance tangent substage; full neutral transfer remains open',
    'scientific_output_hashes': {science_path: science_hash},
    'input_hashes_verified': data['input_sha256'],
    'checks': data['counts'], 'scalar_rows': len(data['scalar_rows']),
    'parameters': data['parameters'], 'environment': data['environment'],
    'parent_manifest': {'path': parent_path, 'sha256': parent_hash},
    'verified_parent_scientific_hashes': inherited,
    'proof_inputs': {p: sha(p) for p in proofs},
    'hashes': {p: sha(p) for p in files},
    'render_hashes': renders, 'build_pass_hashes': passes,
    'build': {'command': './code/build_sector_covariance_tangent.ps1',
              'passes': 2, 'pages': 4, 'no_undefined_or_overfull': True,
              'render': 'pdftoppm -r 90 -png',
              'visual_review': 'All four pages viewed; duplicate equation numbering corrected and affected page rechecked.'},
    'scope': 'Finite exact path expansion and fixed-degree local-size tangent conditional on inherited pure-deficit bounds. Saved-byte audit, not formal proof verification.',
    'remaining': 'Composite covariance/chaos filter, central exact pinching, full neutral series and memory/bilocal return; bulk dictionary and finite G.',
}
with OUT.open('x', encoding='utf-8') as f:
    json.dump(result, f, indent=2)
    f.write('\n')
print(json.dumps({'status': 'passed', 'current_files': len(files),
                  'parent_files': len(inherited), 'sha256': sha(OUT.relative_to(ROOT))}))
