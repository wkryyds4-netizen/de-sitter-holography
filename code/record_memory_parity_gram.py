"""Record the completed parity Gram stage without rerunning its calculations."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import platform
import re
import subprocess

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'results/susskind-memory-parity-gram'
BUILD = ROOT / 'tmp/pdfs/susskind-memory-parity-gram-build'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    target = OUT / 'local-evidence.json'
    if target.exists():
        raise FileExistsError(target)
    records = {}
    for stem in ['signed-moments', 'car-gram', 'limit-checks']:
        data = json.loads((OUT / f'{stem}.json').read_text())
        assert data['status'] == 'passed', stem
        for name, digest in data['input_sha256'].items():
            assert sha(ROOT / name) == digest, name
        records[stem] = data
    assert records['signed-moments']['signed_moment_checks'] == 417
    assert records['signed-moments']['recurrence_checks'] == 139
    assert sum(r['coordinate_pairs'] for r in records['car-gram']['rows']) == 1468
    assert len(records['limit-checks']['records']) == 24
    parentpath = ROOT / 'results/susskind-memory-two-charge/local-evidence.json'
    parent = json.loads(parentpath.read_text())
    unchanged = []
    for name, digest in parent['hashes'].items():
        if name in ['README.md', 'PLAN.md', 'NEXT.md']:
            continue
        assert sha(ROOT / name) == digest, name
        unchanged.append(name)
    log = (OUT / 'main.log').read_text(errors='replace')
    assert re.search(r'5\s+pages', log)
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error', log)
    pdf = ROOT / 'output/pdf/de-sitter-susskind-memory-parity-gram.pdf'
    assert sha(pdf) == sha(BUILD / 'main.pdf')
    renders = [BUILD / f'page-{i}.png' for i in range(1, 6)]
    assert all(p.exists() for p in renders)
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    assert head == parent['current_git_head']
    names = [
        'code/memory_parity_moments.py', 'code/check_memory_parity_gram.py',
        'code/check_memory_parity_limit.py', 'code/build_memory_parity_gram.ps1',
        'code/record_memory_parity_gram.py',
        'manuscript/susskind-memory-parity-gram/main.tex',
        'output/pdf/de-sitter-susskind-memory-parity-gram.pdf',
        'notes/susskind-memory-parity-gram.md',
        'literature/susskind-memory-parity-gram-sources.md',
        'journal/2026-11-06-susskind-memory-parity-gram.md',
        'README.md', 'PLAN.md', 'NEXT.md',
    ]
    paths = [ROOT / n for n in names] + sorted(p for p in OUT.iterdir() if p.is_file())
    result = dict(
        created_at_utc=datetime.now(timezone.utc).isoformat(),
        status='Saved locally; no publication',
        scientific_status='Exact signed moments and the half-filled double-scaled (B, input deficit, output deficit) equal-time Gram. Output residual is positive after B projection and orthogonal to J_H F. Cross with A, weighted A tail, full old-span projection and dynamics remain open.',
        current_git_head=head, git_commit_created=False, git_push_performed=False,
        archive_created=False, parent_manifest_sha256=sha(parentpath),
        previous_stage_scientific_files_unchanged=unchanged,
        checks=dict(signed_moments=417, recurrences=139, full_car_coordinate_pairs=1468,
                    high_precision_cases=24),
        environment=dict(python=platform.python_version(),
                         calculations=[r['environment'] for r in records.values()]),
        build=dict(
            command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_parity_gram.ps1',
            compiler=log.splitlines()[0], passes=2, pages=5,
            render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-parity-gram.pdf tmp/pdfs/susskind-memory-parity-gram-build/page',
            visual_review='All five final page PNGs viewed at 95 dpi during this stage before recording. Equations, references and pagination readable without clipping or overlapping.'),
        failures=[dict(record='results/susskind-memory-parity-gram/verification-issues.txt',
                       cause='Negative integer exponents in parity introduced floats into Fraction sums.',
                       resolution='Integer conditional parity; formulas and tolerances unchanged; final exact checks passed.')],
        process_state='All computations, compilation and rendering terminal.',
        evidence_scope='Exact finite checks support the accompanying analytic derivation. Finite high-precision values concern Delta x/N, not S itself, and are not interval certificates. No full Poisson law or time-dependent spectrum asserted.',
        hashes={p.relative_to(ROOT).as_posix(): sha(p) for p in paths},
        render_hashes={p.relative_to(ROOT).as_posix(): sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt': sha(BUILD / f'pass{i}.txt') for i in [1, 2]},
    )
    target.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(status='saved locally', files=len(paths), pages=5,
                         prior_scientific_files_unchanged=len(unchanged), git_head=head)))


if __name__ == '__main__':
    main()
