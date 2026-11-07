"""Verify saved calculations and preserve the local ordering-metric stage."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, platform, re, subprocess
import mpmath as mp

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'results/susskind-ordering-metric'
BUILD = ROOT / 'tmp/pdfs/susskind-ordering-metric-build'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    target = OUT / 'local-evidence.json'
    if target.exists():
        raise FileExistsError(target)
    parent_path = ROOT / 'results/susskind-commutator-primary/local-evidence.json'
    parent = json.loads(parent_path.read_text())
    previous = []
    for stage in [p['stage'] for p in parent['previous_stages']] + ['susskind-commutator-primary']:
        path = ROOT / f'results/{stage}/local-evidence.json'
        record = json.loads(path.read_text())
        kept = []
        for name, value in record['hashes'].items():
            if name in ['README.md', 'PLAN.md', 'NEXT.md']:
                continue
            assert sha(ROOT / name) == value, name
            kept.append(name)
        previous.append(dict(stage=stage, manifest_sha256=sha(path), unchanged_verified=kept))
    exact = json.loads((OUT / 'exact.json').read_text())
    checks = json.loads((OUT / 'checks.json').read_text())
    for record in [exact, checks]:
        assert record['status'] == 'passed'
        for name, value in record['input_sha256'].items():
            assert sha(ROOT / name) == value, name
    assert len(exact['checks']) == 50 and all(exact['checks'].values())
    assert exact['relative_labels'] == [12, 22, 34, 48, 64, 82]
    assert exact['undressed_phase_norm']['quadratic'] == '181'
    assert len(checks['integrals']) == 5 and len(checks['phase_norm']) == 4
    assert sum(len(g['records']) for g in checks['finite_gram']) == 15
    mp.mp.dps = 120
    assert max(mp.mpf(r['absolute_error']) for r in checks['integrals']) < mp.mpf('5.15e-85')
    assert max(mp.mpf(r['difference']) for r in checks['phase_norm']) < mp.mpf('2.90e-70')
    log = (OUT / 'main.log').read_text(errors='replace')
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error', log)
    assert '(6 pages,' in log
    pdf = ROOT / 'output/pdf/de-sitter-susskind-ordering-metric.pdf'
    assert sha(pdf) == sha(BUILD / 'main.pdf')
    renders = [BUILD / f'page-{i}.png' for i in range(1, 7)]
    assert all(p.exists() for p in renders)
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    assert head == parent['current_git_head']
    names = ['code/ordering_metric.py', 'code/check_ordering_metric.py',
             'code/build_ordering_metric.ps1', 'code/record_ordering_metric.py',
             'manuscript/susskind-ordering-metric/main.tex',
             'output/pdf/de-sitter-susskind-ordering-metric.pdf',
             'notes/susskind-ordering-metric.md', 'literature/susskind-ordering-metric-sources.md',
             'journal/2026-11-06-susskind-ordering-metric.md', 'README.md', 'PLAN.md', 'NEXT.md']
    files = [ROOT / n for n in names] + sorted(p for p in OUT.iterdir() if p.is_file())
    result = dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
        status='Local stage saved; no publication or global freeze',
        scientific_status='Leading Casimir form and conditional free-translation/vertex incompatibility. Physical spatial dictionary remains open.',
        previous_stages=previous, current_git_head=head, git_commit_created=False,
        git_push_performed=False, archive_created=False, global_project_audit_performed=False,
        checks=dict(symbolic=50, kernel_integrals=5, finite_gram=15, phase_norm=4),
        environment=dict(python=platform.python_version(), calculation_environments=[exact['environment'], checks['environment']]),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_ordering_metric.ps1',
            compiler=log.splitlines()[0], passes=2, pages=6, visual_review='All six final rendered pages inspected',
            render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-ordering-metric.pdf tmp/pdfs/susskind-ordering-metric-build/page'),
        retained_failures=['Initial AA macro collision', 'Subsequent 5.98149pt and 0.54567pt line overflows; corrected'],
        error_scope=checks['error_scope'], process_state='All computation, compilation and rendering processes terminal',
        hashes={p.relative_to(ROOT).as_posix(): sha(p) for p in files},
        local_render_hashes={p.relative_to(ROOT).as_posix(): sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt': sha(BUILD / f'pass{i}.txt') for i in [1, 2]})
    target.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(status='saved locally', pages=6, files=len(files),
                         prior_scientific_files=sum(len(p['unchanged_verified']) for p in previous), git_head=head)))

if __name__ == '__main__':
    main()
