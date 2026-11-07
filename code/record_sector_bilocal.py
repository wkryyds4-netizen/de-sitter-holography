"""Seal saved bilocal evidence; never execute scientific drivers."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
STAGE = 'susskind-sector-bilocal'
OUT = ROOT / 'results' / STAGE / 'local-evidence.json'
if OUT.exists():
    raise FileExistsError(OUT)

def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()

def read(path):
    return json.loads((ROOT / path).read_text(encoding='utf-8-sig'))

science = f'results/{STAGE}/interface-checks.json'
assert sha(science) == '0442b1715b4d80fb55703aaf6f474773e98f7708765e2bb4cd6e2ef671ea80f3'
data = read(science)
assert data['status'] == 'passed'
assert data['counts'] == dict(source_pairings=124, colored_words=1640,
                             replica_words=400, number_remainders=105,
                             charged_rows=10, bilocal_times=6, bilocal_bounds=18)
for path, digest in data['input_sha256'].items():
    assert sha(path) == digest, path

parents, verified, historical = {}, {}, {}
mutable = {'README.md', 'PLAN.md', 'NEXT.md', 'PROJECT.md', 'AGENTS.md'}

def verify(path, digest, origin):
    if path in mutable or path.startswith('journal/'):
        historical.setdefault(path, []).append(dict(manifest=origin, recorded_sha256=digest))
        return
    assert sha(path) == digest, (origin, path)
    verified[path] = digest

for stage in ['susskind-sector-memory-return', 'susskind-sector-gradient-tail',
              'susskind-sector-neutral-transfer', 'susskind-central-ordinary',
              'susskind-slater-higher-wick', 'susskind-slater-polynomial-density',
              'susskind-memory-two-component', 'susskind-bilocal-reconstruction']:
    path = f'results/{stage}/local-evidence.json'
    parent = read(path)
    parents[path] = dict(sha256=sha(path))
    verified[path] = sha(path)
    for key in ['hashes', 'proof_inputs', 'input_hashes_verified',
                'scientific_output_hashes', 'verified_parent_scientific_hashes']:
        for name, digest in parent.get(key, {}).items():
            if isinstance(digest, dict):
                for nested_path, nested_digest in digest.items():
                    verify(nested_path, nested_digest, path + ':' + name)
                continue
            verify(name, digest, path)
            if name.endswith('.json') and name not in mutable and not name.startswith('journal/'):
                obj = read(name)
                if isinstance(obj, dict):
                    for p, h in obj.get('input_sha256', {}).items():
                        verify(p, h, name)

proofs = ['notes/susskind-bilocal-fourth-tail-operator-draft.md',
          'notes/susskind-bilocal-reconstruction-proof.md',
          'notes/susskind-bilocal-full-four-point-transfer.md',
          'manuscript/susskind-memory-uniform-bound/main.tex']
pdf = f'output/pdf/de-sitter-{STAGE}.pdf'
build, render = f'tmp/pdfs/{STAGE}-build', f'tmp/pdfs/{STAGE}-render'
assert sha(pdf) == sha(f'{build}/main.pdf')
log = (ROOT / f'results/{STAGE}/main.log').read_text(errors='replace')
assert 'Output written on' in log
assert not re.search(r'undefined|Overfull|Fatal error|Rerun to get', log)
passes = {}
for i in (1, 2):
    path = f'{build}/pass{i}.txt'
    raw = (ROOT / path).read_bytes()
    text = raw.decode('utf-16') if raw.startswith((b'\xff\xfe', b'\xfe\xff')) else raw.decode()
    assert 'Output written on' in text and 'Fatal error' not in text
    passes[path] = sha(path)
info = subprocess.check_output(['pdfinfo', str(ROOT / pdf)]).decode(errors='replace')
assert int(re.search(r'Pages:\s+(\d+)', info).group(1)) == 10
renders = {}
for i in range(1, 11):
    path = f'{render}/page-{i:02d}.png'
    assert (ROOT / path).read_bytes().startswith(b'\x89PNG')
    renders[path] = sha(path)
files = ['code/sector_bilocal.py', 'code/check_sector_bilocal.py',
         'code/record_sector_bilocal.py', 'code/build_sector_bilocal.ps1',
         f'manuscript/{STAGE}/main.tex', f'notes/{STAGE}.md',
         f'literature/{STAGE}-sources.md', pdf, f'results/{STAGE}/main.log', science]
result = dict(
    status='Fixed-charge bilocal stage completed; full research goal active',
    created_at_utc=datetime.now(timezone.utc).isoformat(),
    scientific_scope='Compact-central actual memory bounds and O_r(p^-1) bilocal reconstruction; half-filled colored regression, common local memory, and finite covariant bilocal blocks between controlled fixed polynomial preparations.',
    scientific_output_hashes={science: sha(science)}, input_hashes_verified=data['input_sha256'],
    checks=data['counts'], parameters=data['parameters'], environment=data['environment'],
    parameter_note='C_m=1 is the CAR-fixture normalization; the independent degenerate bath fixture is not a SYK ensemble.',
    previous_evidence=parents, verified_parent_scientific_hashes=verified,
    historical_mutable_entrances=historical, proof_inputs={p: sha(p) for p in proofs},
    hashes={p: sha(p) for p in files}, render_hashes=renders, build_pass_hashes=passes,
    build=dict(command='./code/build_sector_bilocal.ps1', passes=2, pages=10,
               render='pdftoppm -r 90 -png', no_undefined_or_overfull=True,
               visual_review='All ten final rendered pages viewed; equations, references and page boundaries checked.'),
    diagnostics='Saved successful science verified by hash, not rerun. First recorder attempt stopped before writing output on a nested historical input-hash schema; nested maps are now read explicitly. Rational algebra and finite float64 diagnostics are not interval certificates or double-scaling fits.',
    remaining='Arbitrary selected factors interlaced inside nonpolynomial blocks; static-observer spatial/mass dictionary and finite-G dynamics. No rate for colored conditional regression.',
    git_commit_created=False, git_push_performed=False, external_reply_sent=False,
    recorder_scope='Saved-byte integrity and build verification, not proof-assistant certification.')
with OUT.open('x', encoding='utf-8') as stream:
    json.dump(result, stream, indent=2)
    stream.write('\n')
print(json.dumps(dict(status='passed', current_files=len(files),
                     parent_scientific_files=len(verified), sha256=sha(OUT.relative_to(ROOT)))))
