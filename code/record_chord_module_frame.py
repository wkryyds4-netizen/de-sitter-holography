"""One-time read-only seal of existing science, proof and rendered artifact."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, re, subprocess

ROOT = Path(__file__).resolve().parents[1]
STAGE = 'susskind-chord-module-frame'
OUT = ROOT / f'results/{STAGE}/local-evidence.json'
if OUT.exists():
    raise FileExistsError(OUT)

def sha(p):
    return hashlib.sha256((ROOT / p).read_bytes()).hexdigest()

def read(p):
    return json.loads((ROOT / p).read_text(encoding='utf-8-sig'))

science = {
    f'results/{STAGE}/frame-checks.json': 'f2227eca8b2c57f0d9292773bbf04046fb0c4687b8971eb98ab9625bc5615e27',
    f'results/{STAGE}/null-checks.json': '4207be310055a890a839c60c63219920536c230e937d21fc184a9ef441160fb1',
}
inputs = {}
for p, digest in science.items():
    assert sha(p) == digest, p
    result = read(p)
    assert result['status'] == 'passed'
    for src, expected in result['input_sha256'].items():
        assert sha(src) == expected, src
        inputs[src] = expected
d = read(f'results/{STAGE}/frame-checks.json')
n = read(f'results/{STAGE}/null-checks.json')
assert sum(d['exact_counts'].values()) == 2109
assert len(d['rows']) == 48 and len(d['precision_50_decimal']) == 8
assert len(n['rows']) == 8 and len(n['exact_positivity_forms']) == 3

verified, historical = {}, {}
mutable = {'README.md', 'PLAN.md', 'NEXT.md', 'PROJECT.md', 'AGENTS.md'}
def verify(p, h, origin):
    if isinstance(h, dict):
        for sub, digest in h.items():
            verify(sub, digest, origin)
    elif p in mutable or p.startswith('journal/'):
        historical.setdefault(p, []).append(dict(origin=origin, recorded_sha256=h))
    else:
        assert sha(p) == h, (origin, p)
        verified[p] = h

parent = 'results/susskind-thermal-chord-window/local-evidence.json'
assert sha(parent) == 'f5909c457062c3d9dade9871528a422ba7a44dadc1313a02671785cf8bb98234'
verified[parent] = sha(parent)
for key in ('hashes', 'proof_inputs', 'verified_parent_scientific_hashes',
            'scientific_output_hashes', 'input_hashes_verified'):
    for p, h in read(parent).get(key, {}).items():
        verify(p, h, parent)

pdf = f'output/pdf/de-sitter-{STAGE}.pdf'
build = f'tmp/pdfs/{STAGE}-build'
render = f'tmp/pdfs/{STAGE}-render'
assert sha(pdf) == sha(f'{build}/main.pdf')
log = (ROOT / f'results/{STAGE}/main.log').read_text(errors='replace')
assert 'Output written on' in log
assert not re.search('undefined|Overfull|Rerun to get|Fatal error', log)
info = subprocess.check_output(['pdfinfo', str(ROOT / pdf)]).decode(errors='replace')
assert int(re.search(r'Pages:\s+(\d+)', info).group(1)) == 7
passes = {}
for i in (1, 2):
    p = f'{build}/pass{i}.txt'
    raw = (ROOT / p).read_bytes()
    txt = raw.decode('utf-16') if raw.startswith((b'\xff\xfe', b'\xfe\xff')) else raw.decode()
    assert 'Output written on' in txt and 'Fatal error' not in txt
    passes[p] = sha(p)
renders = {}
for i in range(1, 8):
    p = f'{render}/page-{i}.png'
    assert (ROOT / p).read_bytes().startswith(b'\x89PNG')
    renders[p] = sha(p)

files = [
    'code/chord_module_frame.py', 'code/check_chord_module_frame.py',
    'code/check_chord_module_null.py', 'code/build_chord_module_frame.ps1',
    'code/record_chord_module_frame.py', f'manuscript/{STAGE}/main.tex',
    f'notes/{STAGE}.md', f'literature/{STAGE}-sources.md', pdf,
    f'results/{STAGE}/main.log', *science,
]
result = dict(
    status='Moving-core finite-word limit and positive null candidate established; overall goal active',
    created_at_utc=datetime.now(timezone.utc).isoformat(),
    scope='Fixed-lambda half-filled tracial one-matter DS module, then epsilon->0 and m->infinity with epsilon*m->0; actual charged curvature d=2. True thermal endpoint uses stronger epsilon*m^2->0.',
    scientific_output_hashes=science, input_hashes_verified=inputs,
    exact_counts=d['exact_counts'], positivity_form_count=len(n['exact_positivity_forms']),
    numerical_summary=dict(
        moving_cases=len(d['rows']), null_cases=len(n['rows']),
        max_input_tail_squared_bound=max(r['input_tail_squared_bound'] for r in d['rows']),
        max_word_refinement_absolute=max(r['word_refinement_absolute'] for r in d['rows']),
        max_50_decimal_relative=max(float(r['relative']) for r in d['precision_50_decimal']),
        max_commutator_cancellation={key: max(r['residuals'][key]['cancellation'] for r in d['rows'])
                                     for key in ('comm_BP', 'comm_EP', 'comm_BE')},
        max_null_cancellation=max(r['cancellation'] for r in n['rows']),
    ),
    limitations=[
        'Input tail is not a complete unbounded-word numerical error certificate; analytic finite-path tail proof is separate.',
        'EP and BE tiny residuals have cancellation up to 4.03e8 and 1.52e9, respectively; the eight-word high-precision check does not independently certify these residuals.',
        'Finite-q Q is not asserted positive. Limit positivity is a standard Laguerre/discrete-series form.',
        'Full rho bound concerns specified left vector insertion; original selected correlator ordering is still open.',
        'No opposite null momentum, mass determination, multiparticle law or finite-G dynamics established.',
    ],
    environment=d['environment'], previous_evidence={parent: sha(parent)},
    verified_parent_scientific_hashes=verified, historical_mutable_entrances=historical,
    hashes={p: sha(p) for p in files}, build_pass_hashes=passes, render_hashes=renders,
    build=dict(command='./code/build_chord_module_frame.ps1', passes=2, pages=7,
               render='pdftoppm -r 90 -png',
               visual_review='All seven final pages viewed after final edits; formulas, table, footnote, references and margins checked.',
               no_undefined_or_overfull=True),
    remaining='Read established multiparticle coproduct; determine actual Q composition and opposite null action; retain selected ordering and subleading U; mass, four-bilocal and finite-G.',
    git_commit_created=False, git_push_performed=False, external_reply_sent=False,
)
with OUT.open('x', encoding='utf-8') as f:
    json.dump(result, f, indent=2)
    f.write('\n')
print(json.dumps(dict(status='passed', current_files=len(files), parent_scientific_files=len(verified),
                      pdf_sha256=sha(pdf), evidence_sha256=sha(OUT.relative_to(ROOT)))))
