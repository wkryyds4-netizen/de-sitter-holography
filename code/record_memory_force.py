"""Bind the mean-memory stage without committing or changing prior records."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, re, subprocess, platform

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'results/susskind-memory-force'
BUILD = ROOT / 'tmp/pdfs/susskind-memory-force-build'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    target = OUT / 'local-evidence.json'
    if target.exists():
        raise FileExistsError(target)
    parent = json.loads((ROOT / 'results/susskind-bilocal-channel/local-evidence.json').read_text())
    previous = []
    for stage in [r['stage'] for r in parent['previous_stages']] + ['susskind-bilocal-channel']:
        path = ROOT / f'results/{stage}/local-evidence.json'
        record = json.loads(path.read_text())
        kept = []
        for name, digest in record['hashes'].items():
            if name in ['README.md', 'PLAN.md', 'NEXT.md']:
                continue
            assert sha(ROOT / name) == digest, name
            kept.append(name)
        previous.append(dict(stage=stage, manifest_sha256=sha(path), unchanged_verified=kept))
    records = {}
    for name in ['bridge', 'kernel-checks']:
        record = json.loads((OUT / f'{name}.json').read_text())
        assert record['status'] == 'passed'
        for path, digest in record['input_sha256'].items():
            assert sha(ROOT / path) == digest, path
        records[name] = record
    bridge, kernel = records['bridge'], records['kernel-checks']
    assert len(bridge['integer_renewal']) == 3
    assert sum(len(x['renewal_orders']) for x in bridge['integer_renewal']) == 18
    assert sum(len(x['values']) for x in bridge['comparison_bounds']) == 12
    assert len(kernel['tagged_car']['rows']) == 130 and kernel['tagged_car']['nonzero'] == 17
    assert len(kernel['chord_identities']) == 25 and len(kernel['mehler_density_limit']) == 9
    assert len(kernel['curves']) == 4
    log = (OUT / 'main.log').read_text(errors='replace')
    assert re.search(r'7\s+pages', log)
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error', log)
    pdf = ROOT / 'output/pdf/de-sitter-susskind-memory-force.pdf'
    assert sha(pdf) == sha(BUILD / 'main.pdf')
    renders = [BUILD / f'page-{i}.png' for i in range(1, 8)]
    assert all(p.exists() for p in renders)
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    assert head == parent['current_git_head']
    names = ['code/memory_force_bridge.py', 'code/check_force_kernel.py',
             'code/build_memory_force.ps1', 'code/record_memory_force.py',
             'manuscript/susskind-memory-force/main.tex',
             'output/pdf/de-sitter-susskind-memory-force.pdf',
             'notes/susskind-memory-force.md', 'literature/susskind-memory-force-sources.md',
             'journal/2026-11-06-susskind-memory-force.md', 'README.md', 'PLAN.md', 'NEXT.md']
    files = [ROOT / name for name in names] + sorted(p for p in OUT.iterdir() if p.is_file())
    result = dict(
        created_at_utc=datetime.now(timezone.utc).isoformat(),
        status='Local stage saved; no publication or global freeze',
        scientific_status='Compact-time scalar annealed centered-memory limit and elementary 1/p correction; exact renewal coefficients. Full operator fluctuation memory and bulk spatial dictionary remain open.',
        previous_stages=previous, current_git_head=head,
        git_commit_created=False, git_push_performed=False, archive_created=False,
        global_project_audit_performed=False,
        checks=dict(exact_renewal_orders=18, exact_first_coefficients=6,
                    finite_time_bounds=12, tagged_words=130, tagged_nonzero=17,
                    exact_chord_polynomials=25, mehler_limit_points=9, lambda_curves=4),
        environment=dict(python=platform.python_version(),
                         calculation_environments=[r['environment'] for r in records.values()]),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_force.ps1',
                   compiler=log.splitlines()[0], passes=2, pages=7,
                   visual_review='All seven final rendered pages inspected; no clipping, overlap or undefined citations',
                   render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-force.pdf tmp/pdfs/susskind-memory-force-build/page'),
        retained_failures=['Integer helper truncated half-integral coefficients; original source and failure saved; fixed exact dyadic conversion.',
                           'Random tagged-word nonzero coverage insufficient; initial source and failure saved; ten directed words added without lowering threshold.'],
        error_scope={name: record['error_scope'] for name, record in records.items()},
        process_state='All computation, compilation and rendering processes terminal',
        hashes={p.relative_to(ROOT).as_posix(): sha(p) for p in files},
        local_render_hashes={p.relative_to(ROOT).as_posix(): sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt': sha(BUILD / f'pass{i}.txt') for i in [1, 2]})
    target.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(status='saved locally', files=len(files), pages=7,
                          prior_scientific_files=sum(len(x['unchanged_verified']) for x in previous),
                          git_head=head)))


if __name__ == '__main__':
    main()
