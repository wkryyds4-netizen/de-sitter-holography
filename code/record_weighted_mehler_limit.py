"""Seal saved scientific bytes and the reviewed PDF once; no calculation reruns."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, re, subprocess

ROOT = Path(__file__).resolve().parents[1]
STAGE = 'susskind-weighted-mehler-limit'
OUT = ROOT / f'results/{STAGE}/local-evidence.json'
if OUT.exists():
    raise FileExistsError(OUT)
def sha(p):
    return hashlib.sha256((ROOT / p).read_bytes()).hexdigest()
def read(p):
    return json.loads((ROOT / p).read_text(encoding='utf-8-sig'))

science = {
    f'results/{STAGE}/modular-checks.json': '28364b9e42da4844c38fd03c9062ccafcf820c08f4578ca1b81246d56c11c07f',
    f'results/{STAGE}/trace-mixture-checks.json': 'd86cbb8f250434f72d1cc7d39086fafc8a106bfc0f286f31463908e46d802463',
}
inputs = {}
data = []
for p, h in science.items():
    assert sha(p) == h, p
    d = read(p)
    assert d['status'] == 'passed'
    for src, digest in d['input_sha256'].items():
        assert sha(src) == digest, src
        inputs[src] = digest
    data.append(d)
assert [len(data[0][k]) for k in ['product_modular','boundary_layers','local_expansions','fourier']] == [27,6,18,12]
assert len(data[1]['exact_taylor']) == 8
assert sum(len(r['values']) for r in data[1]['saved_packet_analysis']) == 28
parents, verified, historical = {}, {}, {}
mutable = {'README.md','NEXT.md','PLAN.md','PROJECT.md','AGENTS.md'}
def verify(p, h, origin):
    if isinstance(h, dict):
        for sub, digest in h.items():
            verify(sub, digest, origin)
    elif p in mutable or p.startswith('journal/'):
        historical.setdefault(p, []).append(dict(origin=origin, recorded_sha256=h))
    else:
        assert sha(p) == h, (origin, p)
        verified[p] = h
for stage in ['susskind-selected-rindler-response','susskind-sector-bilocal','susskind-slater-selected-state']:
    p = f'results/{stage}/local-evidence.json'
    d = read(p)
    parents[p] = sha(p)
    verified[p] = sha(p)
    for key in ['hashes','proof_inputs','verified_parent_scientific_hashes','scientific_output_hashes','input_hashes_verified']:
        for src, h in d.get(key, {}).items():
            verify(src, h, p)

pdf = f'output/pdf/de-sitter-{STAGE}.pdf'
build = f'tmp/pdfs/{STAGE}-build'
render = f'tmp/pdfs/{STAGE}-render'
assert sha(pdf) == sha(f'{build}/main.pdf')
log = (ROOT / f'results/{STAGE}/main.log').read_text(errors='replace')
assert 'Output written on' in log and not re.search('Overfull|undefined|Rerun to get|Fatal error', log)
passes = {}
for i in (1,2):
    p = f'{build}/pass{i}.txt'
    raw = (ROOT / p).read_bytes()
    text = raw.decode('utf-16') if raw.startswith((b'\xff\xfe',b'\xfe\xff')) else raw.decode()
    assert 'Output written on' in text and 'Fatal error' not in text
    passes[p] = sha(p)
info = subprocess.check_output(['pdfinfo', str(ROOT / pdf)]).decode(errors='replace')
assert int(re.search(r'Pages:\s+(\d+)', info).group(1)) == 7
renders = {}
for i in range(1,8):
    p = f'{render}/page-{i}.png'
    assert (ROOT / p).read_bytes().startswith(b'\x89PNG')
    renders[p] = sha(p)
tex = (ROOT / f'manuscript/{STAGE}/main.tex').read_text()
assert not any(ord(c)<32 and c not in '\n\r\t' for c in tex)
files = ['code/weighted_mehler_limit.py','code/check_weighted_mehler_limit.py',
         'code/check_mehler_trace_mixture.py','code/record_weighted_mehler_limit.py',
         'code/build_weighted_mehler_limit.ps1',f'manuscript/{STAGE}/main.tex',
         f'notes/{STAGE}.md',f'literature/{STAGE}-sources.md',*science,pdf,f'results/{STAGE}/main.log']
diagnostics = [f'{build}/initial-failed-pass1.txt', f'{build}/initial-failed-main.log']
result = dict(
    status='Weighted compact-time Mehler expansion proved in the manuscript; full research goal active',
    created_at_utc=datetime.now(timezone.utc).isoformat(),
    theorem_scope='Fixed positive integer d, fixed compact real time and fixed finite derivative order: trace g-epsilon R/2 and energy-tilted g+epsilon R/2, both with O(epsilon^2) remainder. Physical selected interpretation here is d=2 curvature.',
    limits='Finite-N exact projection, fixed-lambda double scaling, then epsilon=4lambda to zero at fixed string time; no joint limit or N-uniform rate.',
    scientific_output_hashes=science, input_hashes_verified=inputs,
    environment=[d['environment'] for d in data],
    parameters=dict(d=[1,2,3], product_epsilon=['0.4','0.1','0.04'],fourier_u=[0,0.5,1,2],fourier_cutoff=32,precision_digits=[50,70],paired_images=8),
    checks=dict(product_modular=27,boundary_layers=6,local_expansions=18,fourier=12,exact_trace_coefficients=16,saved_packet_times=28),
    previous_evidence=parents, verified_parent_scientific_hashes=verified,
    historical_mutable_entrances=historical, hashes={p:sha(p) for p in files},
    render_hashes=renders,build_pass_hashes=passes,initial_build_diagnostics={p:sha(p) for p in diagnostics},
    build=dict(command='./code/build_weighted_mehler_limit.ps1',passes=2,pages=7,render='pdftoppm -r 90 -png',no_undefined_or_overfull=True,visual_review='All seven final pages viewed; equations, table, references and page boundaries checked.'),
    remaining='Mixed e,W source for full bilocal, spatial/mass generators, arbitrary selected interlacings, finite-G dynamics; noninteger d and growing times are outside this theorem.',
    recorder_scope='Read-only saved-byte verification; no successful scientific drivers or previous recorders executed. Numerical quadrature is not interval certified. Publication is tracked separately.')
with OUT.open('x', encoding='utf-8') as f:
    json.dump(result, f, indent=2)
    f.write('\n')
print(json.dumps(dict(status='passed',current_files=len(files),parent_scientific_files=len(verified),sha256=sha(OUT.relative_to(ROOT)))))
