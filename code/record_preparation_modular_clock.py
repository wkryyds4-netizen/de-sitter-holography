"""One-time read-only seal for modular-clock results; never reruns science."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, re, subprocess

ROOT = Path(__file__).resolve().parents[1]
STAGE = 'susskind-preparation-modular-clock'
OUT = ROOT/f'results/{STAGE}/local-evidence.json'
if OUT.exists():
    raise FileExistsError(OUT)
def sha(p): return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p): return json.loads((ROOT/p).read_text(encoding='utf-8-sig'))
science = f'results/{STAGE}/clock-checks.json'
assert sha(science) == 'eef7dc91d4280ee0d3955aea75de20a536a243931f2916f3fb941d103da65cc8'
d = read(science)
assert d['status'] == 'passed' and len(d['rows']) == 22 and len(d['energy_windows']) == 4
for p, h in d['input_sha256'].items(): assert sha(p) == h, p
assert all(r['fine']['ratio'] >= r['analytic_ratio_lower'] > 1 for r in d['energy_windows'])
verified, parents, historical = {}, {}, {}
mutable = {'README.md', 'PLAN.md', 'NEXT.md', 'PROJECT.md', 'AGENTS.md'}
def verify(p, h, origin):
    if isinstance(h, dict):
        for sub, digest in h.items(): verify(sub, digest, origin)
    elif p in mutable or p.startswith('journal/'):
        historical.setdefault(p, []).append(dict(origin=origin, recorded_sha256=h))
    else:
        assert sha(p) == h, (origin, p)
        verified[p] = h
for stage in ['susskind-selected-bilocal-covariance', 'susskind-slater-selected-state',
              'susskind-slater-higher-wick']:
    p = f'results/{stage}/local-evidence.json'
    parents[p] = sha(p); verified[p] = sha(p)
    for key in ['hashes', 'proof_inputs', 'verified_parent_scientific_hashes',
                'scientific_output_hashes', 'input_hashes_verified']:
        for src, h in read(p).get(key, {}).items(): verify(src, h, p)
assert parents['results/susskind-selected-bilocal-covariance/local-evidence.json'] == '735aa31bd3013df3ced7f6761ee83d9c8f9b54a3a53a7bbc23f7961e975a3bd9'
proof = ['notes/susskind-slater-selected-state.md', 'notes/susskind-slater-higher-wick.md',
         'notes/susskind-translations.md', 'notes/susskind-rindler-packet.md',
         'notes/susskind-selected-bilocal-covariance.md']
pdf = f'output/pdf/de-sitter-{STAGE}.pdf'
build, render = f'tmp/pdfs/{STAGE}-build', f'tmp/pdfs/{STAGE}-render'
assert sha(pdf) == sha(f'{build}/main.pdf')
log = (ROOT/f'results/{STAGE}/main.log').read_text(errors='replace')
assert 'Output written on' in log and not re.search('undefined|Overfull|Rerun to get|Fatal error', log)
passes = {}
for i in (1, 2):
    p = f'{build}/pass{i}.txt'; raw = (ROOT/p).read_bytes()
    txt = raw.decode('utf-16') if raw.startswith((b'\xff\xfe', b'\xfe\xff')) else raw.decode()
    assert 'Output written on' in txt and 'Fatal error' not in txt
    passes[p] = sha(p)
info = subprocess.check_output(['pdfinfo', str(ROOT/pdf)]).decode(errors='replace')
assert int(re.search(r'Pages:\s+(\d+)', info).group(1)) == 6
renders = {}
for i in range(1, 7):
    p = f'{render}/page-{i}.png'
    assert (ROOT/p).read_bytes().startswith(b'\x89PNG')
    renders[p] = sha(p)
files = ['code/preparation_modular_clock.py', 'code/check_preparation_modular_clock.py',
         'code/build_preparation_modular_clock.ps1', 'code/record_preparation_modular_clock.py',
         f'manuscript/{STAGE}/main.tex', f'notes/{STAGE}.md',
         f'literature/{STAGE}-sources.md', science, pdf, f'results/{STAGE}/main.log']
result = dict(
    status='Explicit selected modular spectral kernel and finite-window Hamiltonian non-KMS witness established in the manuscript; overall research goal active',
    created_at_utc=datetime.now(timezone.utc).isoformat(),
    scope='Fixed-lambda half-filled DS covariant algebra with established bounded invertible selected density; standard tracial modular theory. No new finite-N modular convergence, spatial generator, small-lambda modular limit or finite-G result.',
    inherited_not_new='Selected density, q-Gaussian rotation, positivity and nonstationarity commutator were already established in slater-selected-state; no old scientific drivers rerun.',
    scientific_output_hashes={science:sha(science)}, input_hashes_verified=d['input_sha256'],
    checks=dict(complex_time_points=22, energy_windows=4, independent_representations=['Jacobi','angular products'],
                orders=[128,192], window_orders=[32,64], product_cutoff=192),
    numerical_summary=dict(max_cross_absolute=max(r['cross_absolute'] for r in d['rows']),
        max_cross_relative=max(r['cross_relative'] for r in d['rows']),
        max_refinement_absolute=max(r['refinement_absolute'] for r in d['rows']),
        max_cancellation=max(r['cancellation'] for r in d['rows'])),
    limitations=d['scope'], environment=d['environment'],
    proof_inputs={p:sha(p) for p in proof}, previous_evidence=parents,
    verified_parent_scientific_hashes=verified, historical_mutable_entrances=historical,
    hashes={p:sha(p) for p in files}, build_pass_hashes=passes, render_hashes=renders,
    build=dict(command='./code/build_preparation_modular_clock.ps1', passes=2, pages=6,
               render='pdftoppm -r 90 -png', visual_review='All six final pages viewed; equations, references, tables and margins checked.',
               no_undefined_or_overfull=True),
    remaining='Thermally regulated chord preparations and the smooth large-chord limit for spatial generators, with original Hamiltonian boost and actual selected preparation kept explicit; four-bilocal and finite-G dynamics.',
    git_commit_created=False, git_push_performed=False, external_reply_sent=False)
with OUT.open('x', encoding='utf-8') as f:
    json.dump(result, f, indent=2); f.write('\n')
print(json.dumps(dict(status='passed', current_files=len(files), parent_scientific_files=len(verified),
                     pdf_sha256=sha(pdf), evidence_sha256=sha(OUT.relative_to(ROOT)))))
