"""Read-only provenance recorder; never executes scientific generators."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, platform, re, subprocess

ROOT = Path(__file__).resolve().parent.parent
STAGE = 'susskind-slater-transfer'
OUT = ROOT / f'results/{STAGE}/local-evidence.json'
EXPECTED = {
    'exact-checks.json': '9169dc569325988ee624bc5cae80b0935737cda4b058c7a3025167e338441c6a',
    'canonical-checks.json': 'de1fda36e686f83b7f06798e7c30a2a470e69c52f80ac965361a98e32a953f0e',
    'spectral-checks.json': 'a7bcf2201f3b0bf4dd935a1a459ee6ffbf0b6a60fa103d091ddca0601172e263',
}
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def rel(p): return p.relative_to(ROOT).as_posix()
def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
if OUT.exists(): raise FileExistsError(OUT)
science, inputs, versions, parameters, counts = {}, {}, {}, {}, {}
for name, digest in EXPECTED.items():
    p = ROOT / f'results/{STAGE}/{name}'
    assert sha(p) == digest, name
    data = read(p)
    assert data['status'] == 'passed'
    for path, expected in data['input_sha256'].items():
        assert sha(ROOT / path) == expected, path
        inputs[path] = expected
    science[rel(p)] = digest
    versions[name] = data['environment']
    parameters[name] = hashlib.sha256(json.dumps(data['parameters'], sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    if 'counts' in data: counts[name] = data['counts']
parent_path = 'results/susskind-slater-joint-words/local-evidence.json'
assert sha(ROOT / parent_path) == '869aa11a2827556bb5bd7718558d35232ec46c71971d0296826c47f8b287878b'
parent = read(ROOT / parent_path)
for path, digest in parent['hashes'].items():
    assert sha(ROOT / path) == digest, path
build = ROOT / f'tmp/pdfs/{STAGE}-build'
pdf = ROOT / f'output/pdf/de-sitter-{STAGE}.pdf'
log = ROOT / f'results/{STAGE}/main.log'
assert sha(pdf) == sha(build / 'main.pdf')
logtext = log.read_text(errors='replace')
assert 'Output written on' in logtext
assert not re.search(r'undefined|Overfull|Fatal error|Rerun to get', logtext)
passes = {}
for i in (1, 2):
    p = build / f'pass{i}.txt'
    txt = p.read_text(errors='replace')
    assert 'Output written on' in txt and 'did not succeed' not in txt
    passes[rel(p)] = sha(p)
info = subprocess.run(['pdfinfo', str(pdf)], capture_output=True, check=True).stdout.decode(errors='replace')
assert int(re.search(r'Pages:\s+(\d+)', info).group(1)) == 9
renders = {}
for i in range(1, 10):
    p = ROOT / f'tmp/pdfs/{STAGE}-render/page-{i}.png'
    assert p.read_bytes().startswith(b'\x89PNG\r\n\x1a\n')
    renders[rel(p)] = sha(p)
changed = [Path(row['Path']).name for row in read(build / 'pre-final-render-hashes.json')
           if sha(Path(row['Path'])) != row['Hash'].lower()]
assert changed == ['page-5.png'], changed
figure = ROOT / f'results/{STAGE}/selected-spectrum.png'
assert sha(figure) == '2d0aed06b6286ccf55be25dda914f2e5287f046476cc99b2a4593598832f8dd2'
files = [Path(__file__), pdf, log, figure] + [ROOT / p for p in science]
files += [ROOT / ('code/' + name) for name in (
    'slater_open_transfer.py', 'check_slater_open_transfer.py',
    'check_slater_transfer_canonical.py', 'compute_slater_transfer_spectrum.py',
    'build_slater_transfer.ps1')]
files += [ROOT / f'manuscript/{STAGE}/main.tex', ROOT / f'notes/{STAGE}.md',
          ROOT / f'literature/{STAGE}-sources.md']
files += sorted((ROOT / f'results/{STAGE}/diagnostics').glob('*'))
result = dict(
    status='Local milestone sealed before the user-authorized consolidated publication; long-term goal active',
    created_at_utc=datetime.now(timezone.utc).isoformat(),
    scientific_scope='Same-J annealed Slater open transfer, positive graded metric quotient, q-Jacobi towers, compact selected operator, Type I closure, finite-core undressing, selected spectral weight and single-insertion compression. No stationary observer or finite-G dictionary claim.',
    scientific_output_hashes=science, input_hashes_verified=inputs,
    parameter_hashes=parameters, checks=counts,
    previous_evidence={parent_path: dict(sha256=sha(ROOT / parent_path), verified_scientific_hashes=parent['hashes'])},
    hashes={rel(p): sha(p) for p in files}, render_hashes=renders, build_pass_hashes=passes,
    build=dict(command='./code/build_slater_transfer.ps1', passes=2, pages=9,
               final_script_exit_code_observed=0, no_undefined_references_or_overfull=True,
               review='All nine pages inspected in prior continuation; final unit-vector clarification changed only page 5, re-inspected in this continuation. Other eight render hashes identical.',
               changed_final_pages=changed, render='pdftoppm -r 90 -png'),
    environment=dict(recorder_python=platform.python_version(), scientific_environments=versions),
    diagnostics='Initial float contamination in exact checker preserved; corrected Fraction(1) run passed before sealing. Spectral JSON retains Python -Infinity at t=0 logarithmic zero-tail bound. Floating diagnostics are not interval certification.',
    recorder_scope='Read-only hash, status and build verification; no scientific rerun or independent proof certification.')
with OUT.open('x', encoding='utf-8') as f:
    json.dump(result, f, indent=2)
    f.write('\n')
print(json.dumps(dict(status='passed', science=len(science), parent_files=len(parent['hashes']), files=len(files), pages=9, sha256=sha(OUT))))
