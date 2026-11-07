"""Read saved scientific and build evidence; never rerun scientific drivers."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
STAGE = 'susskind-sector-gradient-tail'
DIR = ROOT/'results'/STAGE
OUT = DIR/'local-evidence.json'
if OUT.exists():
    raise FileExistsError(OUT)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def rel(path):
    return path.relative_to(ROOT).as_posix()


expected = {
    'results/sector-cp-gradient/finite-car-checks.json':
        '255d2ef0993c3b073cd9eab3d7900f9ba7ebc88626b1051dfb1afc3f1b0a087f',
    'results/susskind-sector-gradient-tail/graded-checks.json':
        '719fd3b024e20edb89212a220c58ff15a3af9888c6d25fb9ec7715652f8bb65c',
}
science = {}
inputs = {}
for path, digest in expected.items():
    assert sha(ROOT/path) == digest, path
    data = read(ROOT/path)
    assert data['status'] == 'passed'
    science[path] = data
    for name, h in data.get('input_sha256', {}).items():
        assert sha(ROOT/name) == h, name
        inputs[name] = h
old = science['results/sector-cp-gradient/finite-car-checks.json']
assert old['code_sha256'] == sha(ROOT/'code/check_sector_cp_gradient.py')
inputs['code/check_sector_cp_gradient.py'] = old['code_sha256']
assert old['exact_integer_equalities'] == 1248
assert len(old['psd_checks']) == 168 and len(old['schatten_checks']) == 324
new = science['results/susskind-sector-gradient-tail/graded-checks.json']
assert new['counts'] == {
    'matrix_units': 340, 'graded_derivative_commutators': 5008,
    'repeated_shift_identities': 840, 'normal_monomial_euler': 340,
    'normal_monomial_deletions': 2504, 'nonzero_deletions': 1252}

parentpath = ROOT/'results/susskind-central-ordinary/local-evidence.json'
assert sha(parentpath) == '51eedf4105ae14f9630c08d0e3c0002afe4ab9336bdde5e8e0ee028788108606'
parent = read(parentpath)
verified = {}
for key in ('hashes', 'verified_parent_scientific_hashes'):
    for path, h in parent[key].items():
        assert sha(ROOT/path) == h, path
        verified[path] = h
v = parent['parent_manifest']
assert sha(ROOT/v['path']) == v['sha256']
verified[v['path']] = v['sha256']

pdf = ROOT/f'output/pdf/de-sitter-{STAGE}.pdf'
log = DIR/'main.log'
build = ROOT/f'tmp/pdfs/{STAGE}-build'
render = ROOT/f'tmp/pdfs/{STAGE}-render'
assert sha(pdf) == sha(build/'main.pdf')
txt = log.read_text(errors='replace')
assert 'Output written on' in txt
assert not re.search(r'undefined|Overfull|Fatal error|Rerun to get', txt)
passes = {}
for i in (1, 2):
    p = build/f'pass{i}.txt'
    raw = p.read_bytes()
    text = raw.decode('utf-16') if raw.startswith((b'\xff\xfe', b'\xfe\xff')) else raw.decode(errors='replace')
    assert 'Output written on' in text and 'Fatal error' not in text
    passes[rel(p)] = sha(p)
info = subprocess.check_output(['pdfinfo', str(pdf)]).decode(errors='replace')
assert int(re.search(r'Pages:\s+(\d+)', info).group(1)) == 8
renders = {}
for i in range(1, 9):
    p = render/f'page-{i}.png'
    assert p.read_bytes().startswith(b'\x89PNG')
    renders[rel(p)] = sha(p)
proofs = [
    'notes/susskind-bilocal-fourth-tail-diagram-draft.md',
    'notes/susskind-bilocal-fourth-tail-operator-draft.md',
    'notes/susskind-bilocal-fourth-tail-charged-proof.md',
    'notes/susskind-central-ordinary.md',
    'notes/susskind-slater-polynomial-density.md',
    'literature/susskind-bilocal-clifford-riesz-sources.md',
]
files = [Path(__file__), ROOT/'code/check_sector_graded_gradient.py',
         ROOT/'code/check_sector_cp_gradient.py', ROOT/'code/build_sector_gradient_tail.ps1',
         ROOT/f'notes/{STAGE}.md', ROOT/f'literature/{STAGE}-sources.md',
         ROOT/f'manuscript/{STAGE}/main.tex', pdf, log,
         DIR/'diagnostics/initial-build-failure.txt', DIR/'diagnostics/initial-overfull.log']
files += [ROOT/p for p in expected]
result = {
    'status': 'Local research milestone; full de Sitter goal remains active',
    'created_at_utc': datetime.now(timezone.utc).isoformat(),
    'scientific_scope': 'Sector CP gradient transport; graded half-unit shift; pure marked-Wick gradient bound; rectangular CP product estimates; degree-uniform actual ordinary charged tail on |m-N/2|<=L0 p and finite-time transfer to previously controlled fixed polynomial preparations.',
    'scientific_output_hashes': expected, 'input_hashes_verified': inputs,
    'parameters': {p: d['parameters'] for p, d in science.items()},
    'checks': {'graded': new['counts'], 'channel_exact': 1248,
               'channel_numerical_psd': 168, 'channel_numerical_schatten': 324},
    'parent_manifest': {'path': rel(parentpath), 'sha256': sha(parentpath)},
    'verified_parent_scientific_hashes': verified,
    'historical_mutable_entrances': parent['historical_mutable_entrances'],
    'proof_inputs': {p: sha(ROOT/p) for p in proofs},
    'hashes': {rel(p): sha(p) for p in files},
    'render_hashes': renders, 'build_pass_hashes': passes,
    'build': {'command': './code/build_sector_gradient_tail.ps1',
              'passes': 2, 'pages': 8, 'no_undefined_or_overfull': True,
              'render': 'pdftoppm -r 90 -png',
              'visual_review': 'All eight final rendered pages viewed; equations, references and page boundaries checked.'},
    'environment': {p: d['environment'] for p, d in science.items()},
    'diagnostics': 'Initial build failed on undefined mathscr, then a long Jacobi display overflowed. Both logs retained; command replaced with mathcal and display split. Final two-pass build clean. Scientific drivers each succeeded once and were not rerun.',
    'remaining': 'Actual selected neutral/rank/filter strong bridge and memory-return sector audit; complete selected bilocal; spatial/mass bulk dictionary and finite-G dynamics.',
    'git_commit_created': False, 'git_push_performed': False, 'external_reply_sent': False,
    'recorder_scope': 'Saved-byte integrity, finite-check accounting and PDF build evidence; not a proof-assistant verification of the analytic estimates.'
}
with OUT.open('x', encoding='utf-8') as f:
    json.dump(result, f, indent=2)
    f.write('\n')
print(json.dumps({'status': 'passed', 'current_files': len(files),
                  'parent_scientific_files': len(verified), 'sha256': sha(OUT)}))
