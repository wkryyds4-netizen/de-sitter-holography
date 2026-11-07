"""One-time evidence recorder. Reads saved results; never runs scientific drivers."""
from pathlib import Path
from datetime import datetime, timezone
from fractions import Fraction
import hashlib, json, re, subprocess

ROOT = Path(__file__).resolve().parents[1]
STAGE = 'susskind-selected-bilocal-covariance'
OUT = ROOT / f'results/{STAGE}/local-evidence.json'
if OUT.exists():
    raise FileExistsError(OUT)

def sha(p):
    return hashlib.sha256((ROOT / p).read_bytes()).hexdigest()

def read(p):
    return json.loads((ROOT / p).read_text(encoding='utf-8-sig'))

science = {f'results/{STAGE}/{n}': h for n, h in {
    'exact-checks.json': '2e2249e93ddbb2637676191598e6e2578343de14fdf05b96431d63feef3b8641',
    'gram-checks.json': 'd65c19537be469110424325d7eeb6ce824bad2942dc6ea7f67bf9476df3ba3d2',
    'sign-witness.json': '0f3c3485c4f98197e9e060a7df8630c43e893c5b927314411302effc539a83a6',
}.items()}
data, inputs = [], {}
for p, h in science.items():
    assert sha(p) == h, p
    d = read(p)
    assert d['status'] == 'passed', p
    for src, digest in d['input_sha256'].items():
        assert sha(src) == digest, src
        inputs[src] = digest
    data.append(d)
exact, gram, witness = data
assert exact['counts'] == dict(full_selected=386, odd_boundary=16, even_boundary=16,
                              curvature_marginal=2, positive_residual_Gram=2)
assert all(r['passed'] for r in exact['rows'])
minors = [v for r in exact['rows'] if r['kind'] == 'positive_residual_Gram'
          for v in r['principal_minors']]
assert len(minors) == 14 and all(Fraction(v['value']) >= 0 for v in minors)
assert [len(gram[k]) for k in ['diagonal', 'matrices', 'short_time', 'previous_means']] == [30, 2, 6, 10]
assert Fraction(witness['centered_excess']) < 0 < Fraction(witness['selected_variance'])
assert all(min(m['connected_excess_eigenvalues']) < 0 for m in gram['matrices'])

verified, parents, historical = {}, {}, {}
mutable = {'README.md', 'NEXT.md', 'PLAN.md', 'PROJECT.md', 'AGENTS.md'}
def verify(p, h, origin):
    if isinstance(h, dict):
        for sub, digest in h.items():
            verify(sub, digest, origin)
    elif p in mutable or p.startswith('journal/'):
        historical.setdefault(p, []).append(dict(origin=origin, recorded_sha256=h))
    else:
        assert sha(p) == h, (origin, p)
        verified[p] = h

for stage in ['susskind-selected-bilocal-mean', 'susskind-neutral-selected-response',
              'susskind-sector-bilocal', 'susskind-slater-higher-wick']:
    p = f'results/{stage}/local-evidence.json'
    parents[p] = sha(p)
    verified[p] = sha(p)
    for key in ['hashes', 'proof_inputs', 'verified_parent_scientific_hashes',
                'scientific_output_hashes', 'input_hashes_verified']:
        for src, h in read(p).get(key, {}).items():
            verify(src, h, p)
assert parents['results/susskind-selected-bilocal-mean/local-evidence.json'] == '33349572ca4c4d8b8e2f7964a9b8a3d160a6d4aff3bc597f2d739812fffffc29'
proof = ['notes/susskind-bilocal-multisource-proof.md',
         'notes/susskind-memory-joint-regression-proof.md',
         'notes/susskind-memory-wick-descendants.md',
         'notes/susskind-slater-higher-wick.md',
         'notes/susskind-neutral-selected-response.md',
         'notes/susskind-selected-bilocal-mean.md', 'notes/susskind-sector-bilocal.md']
pdf = f'output/pdf/de-sitter-{STAGE}.pdf'
build, render = f'tmp/pdfs/{STAGE}-build', f'tmp/pdfs/{STAGE}-render'
assert sha(pdf) == sha(f'{build}/main.pdf')
log = (ROOT / f'results/{STAGE}/main.log').read_text(errors='replace')
assert 'Output written on' in log and not re.search('Overfull|undefined|Rerun to get|Fatal error', log)
passes = {}
for i in (1, 2):
    p = f'{build}/pass{i}.txt'
    raw = (ROOT / p).read_bytes()
    txt = raw.decode('utf-16') if raw.startswith((b'\xff\xfe', b'\xfe\xff')) else raw.decode()
    assert 'Output written on' in txt and 'Fatal error' not in txt
    passes[p] = sha(p)
info = subprocess.check_output(['pdfinfo', str(ROOT / pdf)]).decode(errors='replace')
assert int(re.search(r'Pages:\s+(\d+)', info).group(1)) == 8
renders = {}
for i in range(1, 9):
    p = f'{render}/page-{i}.png'
    assert (ROOT / p).read_bytes().startswith(b'\x89PNG')
    renders[p] = sha(p)
files = ['code/selected_bilocal_covariance.py', 'code/check_selected_bilocal_covariance.py',
         'code/bilocal_preparation_gram.py', 'code/check_bilocal_preparation_gram.py',
         'code/check_covariance_sign_witness.py', 'code/build_selected_bilocal_covariance.ps1',
         'code/record_selected_bilocal_covariance.py', f'manuscript/{STAGE}/main.tex',
         f'notes/{STAGE}.md', f'literature/{STAGE}-sources.md', *science,
         pdf, f'results/{STAGE}/main.log']
result = dict(
    status='Two-insertion selected preparation identity and full-bilocal covariance reduction established in the manuscript; full research goal remains active',
    created_at_utc=datetime.now(timezone.utc).isoformat(),
    scope='Energy-centered single-mark fields; actual single-selected half-filled fixed-lambda DS state. Full stationary same-disorder spectral baseline retained. Positive uncentered Gram correction minus mean product; connected excess need not be positive.',
    limits='Exact finite-N pinching, then fixed-lambda double scaling. Common center-time decay is taken after DS, with no rate or finite-N recurrence claim.',
    numerical_scope=gram['scope'], scientific_output_hashes=science, input_hashes_verified=inputs,
    checks=dict(exact=exact['counts'], principal_minors=14, physical_diagonal_times=30,
                unequal_time_matrices=2, short_time_cases=6, old_means_read_only=10),
    numerical_summary=dict(
        max_absolute_refinement=max(r['refinement_absolute'] for r in gram['diagonal']),
        max_relative_refinement=max(r['refinement_relative'] for r in gram['diagonal']),
        max_analytic_degree_bound=max(r['propagated_degree_bound'] for r in gram['diagonal']),
        max_subtraction_condition=max(r['subtraction_condition'] for r in gram['diagonal']),
        max_matrix_refinement=max(m['max_refinement_absolute'] for m in gram['matrices']),
        connected_excess_min_eigenvalues=[min(m['connected_excess_eigenvalues']) for m in gram['matrices']]),
    environment=dict(exact=exact['environment'], gram=gram['environment'], witness_python=witness['python']),
    proof_inputs={p: sha(p) for p in proof}, previous_evidence=parents,
    verified_parent_scientific_hashes=verified, historical_mutable_entrances=historical,
    hashes={p: sha(p) for p in files}, build_pass_hashes=passes, render_hashes=renders,
    build=dict(command='./code/build_selected_bilocal_covariance.ps1', passes=2, pages=8,
               render='pdftoppm -r 90 -png', no_undefined_or_overfull=True,
               visual_review='All eight final pages viewed after absolute-time tail and joint-regression section-reference correction; formulas, references and layout checked. Initial two pass logs retained locally.'),
    remaining='Static-patch state, energy/time and spatial/mass generator dictionary; selected four-bilocal source partitions and finite-G dynamics. Absolute full baseline not numerically integrated.',
    publication='User authorized a consolidated Git push and issue reply; this recorder only seals local evidence before publication.')
with OUT.open('x', encoding='utf-8') as f:
    json.dump(result, f, indent=2)
    f.write('\n')
print(json.dumps(dict(status='passed', current_files=len(files), parent_scientific_files=len(verified),
                      pdf_sha256=sha(pdf), evidence_sha256=sha(OUT.relative_to(ROOT)),
                      numerical_summary=result['numerical_summary'])))
