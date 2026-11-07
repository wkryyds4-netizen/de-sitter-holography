"""Seal the neutral marked-series proof, new checks and inspected PDF."""
from pathlib import Path
import datetime
import hashlib
import json
import platform
import re
import numpy
import mpmath

ROOT = Path(__file__).resolve().parent.parent
STAGE = 'susskind-memory-neutral-marked-series'
OUT = ROOT / 'results' / STAGE


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def comparisons(value):
    if isinstance(value, dict):
        if 'absolute_error' in value:
            yield value
        for child in value.values():
            yield from comparisons(child)
    elif isinstance(value, list):
        for child in value:
            yield from comparisons(child)


def main():
    destination = OUT / 'local-evidence.json'
    if destination.exists():
        raise FileExistsError(destination)
    mutable = {'README.md', 'PLAN.md', 'NEXT.md'}
    parents = {}
    for stage in ['susskind-memory-two-component', 'susskind-memory-uniform-bound',
                  'susskind-memory-wick-rank', 'susskind-memory-central-window']:
        path = ROOT / 'results' / stage / 'local-evidence.json'
        previous = json.loads(path.read_text(encoding='utf-8'))
        checked = []
        for name, digest in previous['hashes'].items():
            if name in mutable:
                continue
            assert sha(ROOT / name) == digest, (stage, name)
            checked.append(name)
        parents[path.relative_to(ROOT).as_posix()] = dict(
            sha256=sha(path), unchanged_scientific_files=checked)
    data = {name: json.loads((OUT / name).read_text(encoding='utf-8')) for name in
            ['coefficient-checks.json', 'rank-checks.json']}
    for name, record in data.items():
        assert record['status'] == 'passed'
        for filename, digest in record['input_sha256'].items():
            assert sha(ROOT / filename) == digest, (name, filename)
    coefficients = data['coefficient-checks.json']
    tests = list(comparisons(coefficients['systems']))
    assert len(tests) == 87
    assert all(test['absolute_error'] <= test['tolerance'] for test in tests)
    maximum_error = max(test['absolute_error'] for test in tests)
    assert maximum_error < 6e-14
    assert sum(len(row['polynomial_decomposition']) for row in coefficients['systems']) == 18
    assert sum(len(row['Wick_sources']) for row in coefficients['systems']) == 12
    assert all(row['low_orders']['eta4_chaos_zero_norm'] > 1e-5
               for row in coefficients['systems'])
    rank = data['rank-checks.json']
    expected = dict(finite_chain_sixth_cases=12,
                    finite_complete_conditional_distributions=38,
                    independent_product_mixtures=6, large_chain_proxy_bounds=54,
                    large_conditional_proxy_bounds=72, exact_slope_domination_cases=450)
    assert rank['counts'] == expected
    assert rank['finite']['raw_casimir_exact_equalities'] == 350
    with mpmath.workdps(120):
        assert mpmath.mpf(rank['slopes']['max_scaled_precision_change']) < mpmath.mpf('2e-78')
    pdf = ROOT / 'output/pdf' / f'de-sitter-{STAGE}.pdf'
    build = ROOT / 'tmp/pdfs' / f'{STAGE}-build'
    log = (OUT / 'main.log').read_text(errors='replace')
    assert re.search(r'\(9 pages,', log)
    assert not re.search(r'undefined|Rerun to get|Overfull|Fatal error', log)
    assert sha(pdf) == sha(build / 'main.pdf')
    assert sha(OUT / 'main.log') == sha(build / 'main.log')
    pages = [build / f'page-{i}.png' for i in range(1, 10)]
    passes = [build / f'pass{i}.txt' for i in [1, 2]]
    assert all(path.exists() for path in pages + passes)
    assert all(path.stat().st_mtime >= pdf.stat().st_mtime for path in pages)
    names = [
        'code/check_neutral_marked_series.py', 'code/check_neutral_marked_rank.py',
        'code/build_neutral_marked_series.ps1', 'code/record_neutral_marked_series.py',
        f'manuscript/{STAGE}/main.tex', f'output/pdf/de-sitter-{STAGE}.pdf',
        f'notes/{STAGE}.md', f'literature/{STAGE}-sources.md',
        f'journal/2026-11-07-{STAGE}.md', 'README.md', 'PLAN.md', 'NEXT.md']
    files = [ROOT / name for name in names] + [OUT / name for name in data] + [OUT / 'main.log']
    result = dict(
        status='Local research progress; long-term objective remains active',
        created_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        scientific_status='Actual neutral finite-time memory field has an absolutely summable highest-chaos rank-marked representation with forced lower-chaos corrections, strong full-charge L2 error on compact real memory times, and exact finite-N energy pinching retained. Joint limiting marked spectra, charged cross terms and the bulk dictionary remain open.',
        parameters='Standard CAR, overlapping complex Gaussian p-body, common original J across all charges, normalized annealed beta=mu=0 full trace, hbar=1/string real time, even N and lambda_N=p^2/N tending to positive lambda, fixed Js. Finite fixtures and all deviations are in the two scientific JSON records.',
        proof_scope='Fixed-degree sixth rank moments extend the old polynomial-ring certificate. A first-commutator graph-norm Taylor tail and fixed-degree double commutator give D^2-weight uniform integrability for the actual exponential. Finite-theta multiplier slope, Jacobi rearrangement, ordered Gaussian-Wick source bounds and factorial tails establish the representation. No uniform D^2 norm of the whole exponential is assumed.',
        continuation='Compute joint marked moments at every fixed degree, retain the charged reverse source, and control zero-frequency mass before passing exact pinching to a spectral limit; a fixed finite list of moments is insufficient.',
        external_scope='Issue 1 previous local review retained. No new GitHub comments checked, no external reply, no full 21143-cell reintegration; the historical 9890 endpoint replay failures remain unchanged.',
        git_commit_created=False, git_push_performed=False,
        archive_created=False, external_reply_sent=False,
        previous_evidence=parents,
        input_hashes_verified={name: record['input_sha256'] for name, record in data.items()},
        checks=dict(coefficient_comparisons=len(tests), maximum_coefficient_error=maximum_error,
                    rank_counts=expected, exact_Casimir_equalities=350,
                    largest_retained_relative_slope_deviation=max(
                        abs(float(row['relative_limit_deviation'])) for row in rank['slopes']['rows']
                        if row['relative_limit_deviation'] is not None),
                    precision_change=rank['slopes']['max_scaled_precision_change']),
        evidence_scopes={name: record['scope'] for name, record in data.items()},
        environment=dict(python=platform.python_version(), numpy=numpy.__version__,
                         mpmath=mpmath.__version__),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_neutral_marked_series.ps1',
                   pages=9, passes=2, compiler_header=log.splitlines()[0],
                   render_command=f'pdftoppm -r 95 -png output/pdf/de-sitter-{STAGE}.pdf tmp/pdfs/{STAGE}-build/page',
                   final_pages_visually_inspected=list(range(1, 10)),
                   visual_result='Final nine pages inspected after the last build. Equations, references, page breaks and margins are readable; no clipping or overlaps.'),
        diagnostics=[
            'Both new scientific drivers passed their first runs. Existing scientific outputs were not overwritten.',
            'First LaTeX build passed. Subsequent successful builds incorporated explicit tails, test statistics, a_j indexing, zero-frequency continuation and clearer reference pagination.',
            'One apply_patch text match failed without changing any file; corrected the target and applied successfully.',
            'Independent mathematical review verified the main proof and corrected the future-work claim from a fixed finite list to every fixed degree plus zero-frequency control.'],
        process_state='Scientific drivers, final compilation and rendering terminal; no live computation handles.',
        hashes={path.relative_to(ROOT).as_posix(): sha(path) for path in files},
        render_hashes={path.relative_to(ROOT).as_posix(): sha(path) for path in pages},
        build_pass_hashes={path.relative_to(ROOT).as_posix(): sha(path) for path in passes})
    with destination.open('x', encoding='utf-8') as handle:
        handle.write(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(dict(status='passed', pages=9, current_files=len(files),
                         previous_scientific_files_verified=sum(
                             len(parent['unchanged_scientific_files']) for parent in parents.values()),
                         coefficient_comparisons=len(tests), maximum_error=maximum_error)))


if __name__ == '__main__':
    main()
