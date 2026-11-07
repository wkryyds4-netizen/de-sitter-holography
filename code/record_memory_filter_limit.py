"""Preserve local sixth-moment, strong-filter, and PDF evidence without publishing."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, platform, re, subprocess

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'results/susskind-memory-filter-limit'
BUILD = ROOT / 'tmp/pdfs/susskind-memory-filter-limit-build'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    target = OUT / 'local-evidence.json'
    if target.exists():
        raise FileExistsError(target)
    records = {}
    for stem in ['sixth-moment-certificate', 'sixth-moment-checks', 'filter-certificate']:
        data = json.loads((OUT / f'{stem}.json').read_text())
        assert data['status'] == 'passed', stem
        for name, digest in data['input_sha256'].items():
            assert sha(ROOT / name) == digest, name
        records[stem] = data
    assert records['sixth-moment-certificate']['total_cancelled_negative_powers'] == 21
    assert len(records['sixth-moment-certificate']['symbolic_rows']) == 6
    assert records['sixth-moment-checks']['finite_moment_equalities'] == 28
    assert len(records['sixth-moment-checks']['records']) == 48
    assert records['filter-certificate']['finite_filter_identities'] == 130
    assert records['filter-certificate']['small_lambda_core_residual_over_Js_eight'] == '1024'
    parentpath = ROOT / 'results/susskind-memory-residual-norm/local-evidence.json'
    parent = json.loads(parentpath.read_text())
    unchanged = []
    for name, digest in parent['hashes'].items():
        if name in ['README.md', 'PLAN.md', 'NEXT.md']:
            continue
        assert sha(ROOT / name) == digest, name
        unchanged.append(name)
    log = (OUT / 'main.log').read_text(errors='replace')
    assert re.search(r'6\s+pages', log)
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error', log)
    pdf = ROOT / 'output/pdf/de-sitter-susskind-memory-filter-limit.pdf'
    assert sha(pdf) == sha(BUILD / 'main.pdf')
    renders = [BUILD / f'page-{i}.png' for i in range(1, 7)]
    assert all(p.exists() for p in renders)
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    assert head == parent['current_git_head']
    names = [
        'code/memory_sixth_moment.py', 'code/check_memory_sixth_moment.py',
        'code/memory_A_filter_uniform.py', 'code/build_memory_filter_limit.ps1',
        'code/record_memory_filter_limit.py',
        'manuscript/susskind-memory-filter-limit/main.tex',
        'output/pdf/de-sitter-susskind-memory-filter-limit.pdf',
        'notes/susskind-memory-filter-limit.md',
        'literature/susskind-memory-filter-limit-sources.md',
        'journal/2026-11-06-susskind-memory-filter-limit.md',
        'README.md', 'PLAN.md', 'NEXT.md',
    ]
    paths = [ROOT / name for name in names] + sorted(p for p in OUT.iterdir() if p.is_file())
    data = dict(
        created_at_utc=datetime.now(timezone.utc).isoformat(),
        status='Saved locally; no publication',
        scientific_status='Uniform ordinary sixth moments at charges m,m-1 and a growing output window prove the A-filter L2 replacement. The complete half-filled quadratic chaos of the quartic memory core has a four-direction expansion and a positive residual relative to B,O,J_hF. Raw Z2, central-charge window, energy-resolved spectra and the bulk dictionary remain open.',
        current_git_head=head, git_commit_created=False, git_push_performed=False, archive_created=False,
        parent_manifest_sha256=sha(parentpath), previous_stage_scientific_files_unchanged=unchanged,
        checks=dict(negative_laurent_coefficients=21, constant_moment_coefficients=6,
                    full_probability_tables=4, finite_moment_identities=28,
                    large_parameter_cases=48, finite_filter_identities=130),
        environment=dict(python=platform.python_version(), calculations=[r['environment'] for r in records.values()]),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_filter_limit.ps1',
                   compiler=log.splitlines()[0], passes=2, pages=6,
                   render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-filter-limit.pdf tmp/pdfs/susskind-memory-filter-limit-build/page',
                   visual_review='All six final pages viewed individually at 95 dpi; readable equations and references, no clipping or overlap.'),
        failures=[], process_state='All computations, compilation and rendering terminal; no background calculation.',
        evidence_scope='Exact polynomial-ring cancellations plus uniform rational analyticity, positive sixth-moment tails, prior coherent input-weighted norm and explicit growing-window proof. Finite decimals are X/N proxy moments, not actual large-system S moments. No full Poisson law, fixed-window A tightness, decimal interval certificate or new particle claim.',
        hashes={p.relative_to(ROOT).as_posix(): sha(p) for p in paths},
        render_hashes={p.relative_to(ROOT).as_posix(): sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt': sha(BUILD / f'pass{i}.txt') for i in [1, 2]},
    )
    target.write_text(json.dumps(data, indent=2) + '\n')
    print(json.dumps(dict(status='saved locally', pages=6, files=len(paths),
                         prior_scientific_files_unchanged=len(unchanged), git_head=head)))


if __name__ == '__main__':
    main()
