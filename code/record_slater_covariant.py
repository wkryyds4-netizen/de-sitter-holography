"""Read-only integrity/build recorder. Does not execute scientific drivers."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, platform, re, subprocess

ROOT = Path(__file__).resolve().parent.parent
STAGE = 'susskind-slater-covariant'
DIR = ROOT / 'results' / STAGE
OUT = DIR / 'local-evidence.json'
if OUT.exists():
    raise FileExistsError(OUT)

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def rel(p):
    return p.relative_to(ROOT).as_posix()

def read(p):
    return json.loads(p.read_text(encoding='utf-8-sig'))

expected = {
    'site-checks.json': '57bb3e63e0dd601e7ca073aae9f80965fbff9bfacc5c1b26327aab8de532c325',
    'word-checks.json': 'dd3f495859d4f7370ca7bf3f04a8b364809f7453ed6ad78454944df418bc65c8',
    'spectral-checks.json': 'b0fc06fd0c1316e24c81ad4785d1de9f5ad466f3f10ee86567d594f9fd994985',
}
science, inputs = {}, {}
for name, digest in expected.items():
    p = DIR / name
    assert sha(p) == digest, name
    data = read(p)
    assert data['status'] == 'passed'
    science[name] = data
    for path, value in data['input_sha256'].items():
        assert sha(ROOT / path) == value, path
        inputs[path] = value
    for path, value in data.get('figure_sha256', {}).items():
        assert sha(ROOT / path) == value, path

parents = {}
for name, digest in (
    ('stationary', '3f48570f38afaf6e384117b440eba63f559d3ddd3774439f6dee9fa3ecac7ff7'),
    ('transfer', '136226a7a720632fcd5bb92eac67a92724fbfd510a6dde70b701bd7faafac58f'),
):
    path = f'results/susskind-slater-{name}/local-evidence.json'
    assert sha(ROOT / path) == digest, path
    parent = read(ROOT / path)
    for file, value in parent['hashes'].items():
        assert sha(ROOT / file) == value, file
    parents[path] = dict(sha256=digest, verified_scientific_hashes=parent['hashes'])

build = ROOT / f'tmp/pdfs/{STAGE}-build'
pdf = ROOT / f'output/pdf/de-sitter-{STAGE}.pdf'
log = DIR / 'main.log'
assert sha(pdf) == sha(build / 'main.pdf')
txt = log.read_text(errors='replace')
assert 'Output written on' in txt
assert not re.search(r'undefined|Overfull|Fatal error|Rerun to get', txt)
passes = {}
for i in (1, 2):
    p = build / f'pass{i}.txt'
    text = p.read_text(errors='replace')
    assert 'Output written on' in text and 'did not succeed' not in text
    passes[rel(p)] = sha(p)
info = subprocess.check_output(['pdfinfo', str(pdf)]).decode(errors='replace')
assert int(re.search(r'Pages:\s+(\d+)', info).group(1)) == 9
renders = {}
for i in range(1, 10):
    p = ROOT / f'tmp/pdfs/{STAGE}-render/page-{i}.png'
    assert p.read_bytes().startswith(b'\x89PNG\r\n\x1a\n')
    renders[rel(p)] = sha(p)
changed = [Path(x['Path']).name for x in read(build / 'pre-final-render-hashes.json')
           if sha(Path(x['Path'])) != x['Hash'].lower()]
assert changed == ['page-5.png'], changed

files = [Path(__file__), ROOT / 'code/build_slater_covariant.ps1',
         ROOT / 'code/slater_covariant_transfer.py',
         ROOT / 'code/check_slater_covariant_sites.py',
         ROOT / 'code/check_slater_covariant_words.py',
         ROOT / 'code/compute_slater_covariant_spectrum.py',
         ROOT / f'notes/{STAGE}.md', ROOT / f'literature/{STAGE}-sources.md',
         ROOT / f'manuscript/{STAGE}/main.tex', pdf, log,
         DIR / 'covariant-channel.png'] + [DIR / name for name in expected]
result = dict(
    status='Local research milestone; full de Sitter goal remains active',
    created_at_utc=datetime.now(timezone.utc).isoformat(),
    scientific_scope='Actual same-J covariant size descendant: all fixed selected mixed words and compact-time convergence, rank replacement, three-letter transfer, positive new direction, distinguished stationary tracial state and its two-color q-Gaussian GNS. No full quotient identification, selected pinching theorem, complete bilocal joint theorem or finite-G dictionary.',
    scientific_output_hashes={rel(DIR / name): digest for name, digest in expected.items()},
    input_hashes_verified=inputs,
    parameter_hashes={name: hashlib.sha256(json.dumps(data['parameters'], sort_keys=True,
                    separators=(',', ':')).encode()).hexdigest() for name, data in science.items()},
    parameters={name: data['parameters'] for name, data in science.items()},
    previous_evidence=parents,
    hashes={rel(p): sha(p) for p in files}, render_hashes=renders, build_pass_hashes=passes,
    checks=dict(CAR_comparisons=sum(x['comparisons'] for x in science['site-checks.json']['fixtures']),
                double_cluster_checks=science['site-checks.json']['cluster_checks'],
                words=science['word-checks.json']['counts'], physical_q_cases=3, time_rows=18),
    build=dict(command='./code/build_slater_covariant.ps1', passes=2, final_exit_code_observed=0,
               pages=9, review='All nine pages viewed at 90 dpi, including plotted figure. Final page 5 re-inspected after equation cross-reference correction; other eight render hashes unchanged.',
               no_undefined_references_or_overfull=True, render='pdftoppm -r 90 -png'),
    environment=dict(recorder_python=platform.python_version(),
                     science={name: data['environment'] for name, data in science.items()}),
    diagnostics='Exact rational transfer fixtures vary gamma and sqrt(lambda) independently; physical spectral fixtures impose lambda=-log(q)/4. Cutoff bound is exact-arithmetic only; floating roundoff is not enclosed. Successful science outputs were not regenerated.',
    git_commit_created=False, git_push_performed=False, external_reply_sent=False,
    recorder_scope='Read-only integrity and build evidence, not independent certification of analytic proofs.')
with OUT.open('x', encoding='utf-8') as stream:
    json.dump(result, stream, indent=2)
    stream.write('\n')
print(json.dumps(dict(status='passed', current_files=len(files),
                     parent_files=sum(len(p['verified_scientific_hashes']) for p in parents.values()),
                     pages=9, sha256=sha(OUT))))
