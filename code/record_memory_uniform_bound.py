"""Seal uniform memory second-moment evidence without rerunning science."""
from pathlib import Path
from decimal import Decimal
import datetime
import hashlib
import json
import platform
import re
import numpy
import scipy
import sympy
import mpmath

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'results/susskind-memory-uniform-bound'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


if __name__ == '__main__':
    destination = OUT / 'local-evidence.json'
    if destination.exists():
        raise FileExistsError(destination)
    mutable = {'README.md', 'PLAN.md', 'NEXT.md'}
    parents = {}
    for stage in ['susskind-memory-charged-transfer',
                  'susskind-memory-wick-rank',
                  'susskind-memory-return-background']:
        name = f'results/{stage}/local-evidence.json'
        path = ROOT / name
        record = json.loads(path.read_text(encoding='utf-8'))
        checked = []
        for item, digest in record['hashes'].items():
            if item in mutable:
                continue
            assert sha(ROOT / item) == digest, item
            checked.append(item)
        parents[name] = dict(sha256=sha(path), unchanged_scientific_files=checked)
    data = {name: json.loads((OUT / name).read_text(encoding='utf-8')) for name in
            ['jacobi-recurrence-checks.json', 'dyson-force-checks.json',
             'jacobi-entire-checks.json']}
    verified = {}
    for name, result in data.items():
        for source, digest in result['input_sha256'].items():
            assert sha(ROOT / source) == digest, (name, source)
        verified[name] = result['input_sha256']

    recurrence = data['jacobi-recurrence-checks.json']['systems']
    assert [(x['N'], x['p'], x['m'], x['max_degree']) for x in recurrence] == [
        (3, 1, 2, 5), (4, 1, 2, 4), (4, 2, 3, 3)]
    cases = [row for system in recurrence for row in system['cases']]
    assert len(cases) == 12
    assert all(row['annihilation_identity_error'] < 2e-9 for row in cases)
    assert all(row['cancelled_top_degree_error'] < 2e-9 for row in cases)
    assert all(row['residual_norm'] <= row['independent_tensor_bound'] + 2e-9
               for row in cases)
    dyson = data['dyson-force-checks.json']
    assert len(dyson['exact_degree_four_quadrature_words']) == 3
    assert all(row['complete_wick_word_error'] < 2e-12
               for row in dyson['exact_degree_four_quadrature_words'])
    rows = dyson['finite_time_rows']
    assert len(rows) == 9
    assert all(row['remainder_norm'] <= row['analytic_remainder_bound'] + 2e-11
               and row['Gaussian_order_difference'] < 2e-11
               and row['time_order_difference'] < 2e-12 for row in rows)
    entire = data['jacobi-entire-checks.json']
    assert entire['exact_geometric_cases'] == 105
    assert len(entire['records']) == 6
    assert all(row['outside_spectral_bound']
               and Decimal(row['absolute_error']) < Decimal('1e-40')
               and Decimal(row['precision_difference']) < Decimal('1e-65')
               for row in entire['records'])

    build = ROOT / 'tmp/pdfs/susskind-memory-uniform-bound-build'
    pdf = ROOT / 'output/pdf/de-sitter-susskind-memory-uniform-bound.pdf'
    log = (OUT / 'main.log').read_text(errors='replace')
    assert re.search(r'\(9 pages,', log)
    assert not re.search(r'undefined|Overfull|Rerun to get|Fatal error', log)
    assert sha(pdf) == sha(build / 'main.pdf')
    assert sha(OUT / 'main.log') == sha(build / 'main.log')
    pngs = [build / f'page-{i}.png' for i in range(1, 10)]
    passes = [build / f'pass{i}.txt' for i in (1, 2)]
    assert all(path.exists() for path in pngs + passes)
    assert all(path.stat().st_mtime >= pdf.stat().st_mtime for path in pngs)
    names = [
        'code/memory_jacobi_bridge.py', 'code/check_memory_dyson_force.py',
        'code/check_memory_jacobi_entire.py', 'code/build_memory_uniform_bound.ps1',
        'code/record_memory_uniform_bound.py',
        'manuscript/susskind-memory-uniform-bound/main.tex',
        'output/pdf/de-sitter-susskind-memory-uniform-bound.pdf',
        'notes/susskind-memory-uniform-bound.md',
        'notes/susskind-memory-joint-kernel-seed.md',
        'literature/susskind-memory-uniform-bound-sources.md',
        'journal/2026-11-07-susskind-memory-uniform-bound.md',
        'README.md', 'PLAN.md', 'NEXT.md']
    paths = [ROOT / name for name in names] + [OUT / name for name in data]
    paths += [OUT / 'main.log'] + sorted((OUT / 'diagnostics').glob('*'))
    result = dict(
        status='Local research progress; overall goal remains active',
        created_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        scientific_status='Uniform annealed L2 bounds for the actual finite-time pinched memory field and ordinary term X0. The two covariances are asymptotically equivalent uniformly on compact memory times and all real external times. The common limiting kernel and bulk dictionary remain open.',
        git_commit_created=False, git_push_performed=False,
        archive_created=False, external_reply_sent=False,
        previous_evidence=parents, verified_scientific_inputs=verified,
        checks=dict(complete_CAR_systems=3, complete_Hermite_cases=12,
                    full_fourth_moment_words=3, finite_time_Dyson_rows=9,
                    exact_geometric_inequalities=105, entire_series_rows=6,
                    max_annihilation_identity_error=max(x['annihilation_identity_error'] for x in cases),
                    max_Gaussian_order_difference=max(x['Gaussian_order_difference'] for x in rows),
                    max_entire_error=str(max(Decimal(x['absolute_error']) for x in entire['records'])),
                    max_entire_precision_difference=str(max(Decimal(x['precision_difference']) for x in entire['records']))),
        environment=dict(python=platform.python_version(), numpy=numpy.__version__,
                         scipy=scipy.__version__, sympy=sympy.__version__, mpmath=mpmath.__version__),
        parameters='Same J across all charge sectors; standard CAR overlapping complex p-body; beta=mu=0 normalized full trace; hbar=1 and string real time; even N, lambda_N in a fixed positive compact interval, fixed Js and S. Finite test parameters and dependency versions appear in the scientific JSONs.',
        evidence_scopes={name: value['scope'] for name, value in data.items()},
        proof_scope='Weighted Gaussian Sobolev Duhamel, even common-copy Dyson expansion, finite-parameter Jacobi chaos-to-energy comparison and whole-complex-plane expansion; binomial charge tail and Cauchy-Schwarz. Small tests do not replace the uniform proof. No path-space tightness or finite-G bulk dynamics claim.',
        continuation_seed_status='Read-only collaborator proposals for four-u leading strong regression and neutral low Taylor coefficients; q quadratic/cubic background algebra checked manually, but no finite-N strong-limit proof or new computation claimed.',
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_uniform_bound.ps1',
                   pages=9, passes=2, compiler_header=log.splitlines()[0],
                   render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-uniform-bound.pdf tmp/pdfs/susskind-memory-uniform-bound-build/page',
                   final_pages_visually_inspected=list(range(1, 10)),
                   visual_result='Final nine pages inspected after the last rebuild: legible equations and references, corrected differentials/central-sector sum, no clipping or overlap.'),
        diagnostics=[
            'All three scientific drivers passed their first execution; no scientific output rerun or overwritten.',
            'Initial compilation succeeded, but visual inspection found a literal qquad and malformed central subscript; corrected with the integration differential spacing. Original source, log and page 8 image retained in diagnostics.',
            'Final mathematical review added a compactness subsequence argument linking the previous convergent-lambda return theorem to the present positive-compact-lambda hypothesis.',
            'Final two-pass build and all-page rendering completed with exit 0; no undefined/overfull/fatal warnings.'],
        process_state='All scientific, build and rendering commands terminal; no live calculation handles.',
        hashes={path.relative_to(ROOT).as_posix(): sha(path) for path in paths},
        render_hashes={path.relative_to(ROOT).as_posix(): sha(path) for path in pngs},
        build_pass_hashes={path.relative_to(ROOT).as_posix(): sha(path) for path in passes})
    with destination.open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(dict(status='passed', pages=9, scientific_outputs=3,
                         previous_files_checked=sum(len(x['unchanged_scientific_files'])
                                                    for x in parents.values()),
                         bound_stage_files=len(paths), checks=result['checks'])))
