"""Seal the Wick-descendant stage after the final PDF has been inspected.

The recorder executes no scientific driver, compiler, renderer or Git command.
Final page names follow Poppler's width determined by the final log page count.
"""
from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import platform
import re

import mpmath as mp
import numpy as np
import sympy as sp


ROOT = Path(__file__).resolve().parent.parent
STAGE = 'susskind-memory-wick-descendants'
OUT = ROOT/'results'/STAGE
SCIENCE = ['descendant-checks.json', 'source-regression-checks.json']
PARENTS = ['susskind-memory-size-bridge', 'susskind-memory-neutral-marked-series',
           'susskind-memory-two-component']
MUTABLE_PLANS = {'README.md', 'PLAN.md', 'NEXT.md'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text(encoding='utf-8'))


def text(path):
    raw = path.read_bytes()
    encoding = 'utf-16' if raw.startswith((b'\xff\xfe', b'\xfe\xff')) else 'utf-8-sig'
    return raw.decode(encoding, errors='replace')


def verify_inputs(record, name):
    assert record['status'] == 'passed', name
    for filename, digest in record['input_sha256'].items():
        path = ROOT/filename
        assert path.is_file(), (name, 'missing input', filename)
        assert sha(path) == digest, (name, 'changed input', filename)
    encoded = json.dumps(record['parameters'], sort_keys=True, separators=(',', ':')).encode()
    assert hashlib.sha256(encoded).hexdigest() == record['parameter_sha256'], (name, 'parameters')


def verify_parents():
    result = {}
    for stage in PARENTS:
        path = ROOT/'results'/stage/'local-evidence.json'
        record = load(path)
        checked = []
        for filename, digest in record['hashes'].items():
            if filename in MUTABLE_PLANS:
                continue
            assert (ROOT/filename).is_file(), (stage, 'missing predecessor', filename)
            assert sha(ROOT/filename) == digest, (stage, 'changed predecessor', filename)
            checked.append(filename)
        result[path.relative_to(ROOT).as_posix()] = dict(
            sha256=sha(path), unchanged_scientific_files=checked,
            excluded_mutable_plans=sorted(MUTABLE_PLANS.intersection(record['hashes'])))
    return result


def verify_descendants(record):
    finite = record['finite_comparisons']
    low = record['low_order_checks']
    source = record['composite_source_regression']['rows']
    centered = record['centered_zero_bridge']['rows']
    moments = record['moment_checks']['rows']
    chaos = record['moment_checks']['eta_chaos_norms']
    actual = dict(finite_OU_Hermite_comparisons=len(finite),
                  prior_V1_V2_calibrations=len(low['prior_V1_V2_moment_calibrations']),
                  composite_source_regressions=len(source),
                  nonvacuous_centered_no_bridge_cases=len(centered),
                  formal_moment_polynomials=len(moments), eta_chaos_norms=len(chaos))
    expected = dict(finite_OU_Hermite_comparisons=42, prior_V1_V2_calibrations=12,
                    composite_source_regressions=10, nonvacuous_centered_no_bridge_cases=6,
                    formal_moment_polynomials=59, eta_chaos_norms=9)
    assert actual == record['counts'] == expected
    parameters = record['parameters']
    assert parameters['finite_systems'] == [[2, 1, 1], [3, 1, 2]]
    assert parameters['finite_maximum_slots'] == 6 and parameters['formal_maximum_slots'] == 8
    assert parameters['maximum_descendant_degree'] == 5 and parameters['C'] == 1
    assert len(parameters['finite_layouts']) == 7
    assert all(row['absolute_error'] <= row['tolerance'] for row in finite)
    maximum_error = max(row['absolute_error'] for row in finite)
    assert maximum_error == record['maximum_finite_absolute_error'] and maximum_error < 1e-13
    assert {(row['N'], row['p'], row['m'], row['Gaussian_coordinates']) for row in finite} == {
        (2, 1, 1, 4), (3, 1, 2, 9)}
    assert all(row['raw_slots'] <= 6 and row['C_sector'] == 1 for row in finite)
    assert all(row['centered_residual'] == '0' and row['uncentered_no_bridge'] != '0'
               for row in centered)
    assert low['eta4_scalar'] == [dict(lambda_power=2, q_power=1, coefficient='-4')]
    assert low['eta4_scalar_squared_norm'] == [dict(lambda_power=4, q_power=2, coefficient='16')]
    assert low['eta4_explicit_difference_norm_squared'] == []
    negative_rows = negative_terms = 0
    for row in moments:
        assert row['maximum_raw_slots'] <= 8
        negative = sum(F(term['coefficient']) < 0 for term in row['polynomial'])
        assert negative == row['negative_coefficients']
        negative_rows += int(negative > 0)
        negative_terms += negative
    assert (negative_rows, negative_terms) == (19, 31)
    for family in record['fields'].values():
        for field in family:
            for term in field['terms']:
                assert 2*term['implicit_C_power']+term['raw_slots'] == field['physical_degree']
                assert all(any(a <= x < y <= b for a, b in term['wick_blocks'])
                           for x, y in term['mark_blocks'])
    return dict(counts=actual, maximum_finite_absolute_error=maximum_error,
                signed_polynomial_rows=negative_rows, retained_negative_coefficients=negative_terms,
                parameter_sha256=record['parameter_sha256'],
                finite_scope='Full-Fock halfsize before charge restriction; floating complete '
                             'Gaussian-Hermite comparisons, not finite rank=size or interval arithmetic.',
                exact_scope='Signed formal polynomial, regression and no-bridge cancellation checks '
                            'at the recorded fixed degrees. The manuscript supplies the all-degree proof.')


def verify_source_regression(record):
    assert record['schema_version'] == 1
    parameters = record['parameters']
    assert parameters['max_degree'] == 7 and parameters['b0_max_slots'] == 8
    assert (parameters['q_exact'], parameters['lambda_exact']) == ('1/3', '1/4')
    assert parameters['rational_point_scope'] == 'independent formal q and lambda, not q=exp(-4 lambda)'
    assert set(record['input_sha256']) == {'code/check_memory_source_regression.py'}
    assert record['environment']['dependencies'] == 'Python standard library only'
    rows, b0 = record['permutation_rows'], record['b0_rows']
    expected = dict(permutations=5913, weighted_mark_samples=80638,
                    regression_cases=56, negative_raw_mark_samples=12316,
                    b0_candidate_pairings=735, b0_pairings=168,
                    b0_negative_centered_terms=60, b0_nonzero_centered_terms=120,
                    b0_groups=38, b0_layouts=7)
    assert record['counts'] == expected
    assert {(row['d'], row['v'], row['orientation']) for row in rows} == {
        (d, v, orientation) for d in range(1, 8) for v in range(1, d+1)
        for orientation in ['prefix', 'suffix']}
    actual = dict(permutations=sum(math.factorial(d) for d in range(1, 8)),
                  weighted_mark_samples=sum(2*d*math.factorial(d) for d in range(1, 8)),
                  regression_cases=len(rows),
                  negative_raw_mark_samples=sum(row['negative_mark_samples'] for row in rows),
                  b0_layouts=len(b0), b0_groups=sum(len(row['groups']) for row in b0),
                  b0_pairings=sum(group['count'] for row in b0 for group in row['groups']),
                  b0_candidate_pairings=sum(math.prod(range(1, row['slots'], 2)) for row in b0))
    assert all(value == expected[key] for key, value in actual.items())

    def multiply(left, right):
        product = [0]*(len(left)+len(right)-1)
        for i, x in enumerate(left):
            for j, y in enumerate(right):
                product[i+j] += x*y
        return product

    qfactorials = {0: [1]}
    for degree in range(1, 8):
        qfactorials[degree] = multiply(qfactorials[degree-1], [1]*degree)
    for row in rows:
        d, v = row['d'], row['v']
        assert row['exact_identity'] is True
        assert (row['minimum_k'], row['maximum_k']) == (v*(2-v), v*v)
        assert 0 <= row['negative_mark_samples'] <= math.factorial(d)
        weighted = row['weighted_q_coefficients']
        assert len(weighted) == d*(d-1)//2+1
        marked = [(v*v-4*degree)*value for degree, value in enumerate(qfactorials[v])]
        assert multiply(weighted, qfactorials[v]) == row['cleared_polynomial']
        assert row['cleared_polynomial'] == multiply(qfactorials[d], marked)

    q, lam = F(parameters['q_exact']), F(parameters['lambda_exact'])
    assert [row['layout'] for row in b0] == parameters['b0_layouts']
    negative_examples = 0
    for row in b0:
        layout = row['layout']
        assert row['slots'] == sum(layout[key] for key in ['d1', 'd2', 'a', 'b']) == 8
        assert row['exact_factorization'] is True and row['exact_group_centering'] is True
        kappas = [F(row['kappa1']), F(row['kappa2'])]
        for key, kappa in zip(['v1', 'v2'], kappas):
            v = layout[key]
            factorial = qfactorials[v]
            denominator = sum(F(value)*q**degree for degree, value in enumerate(factorial))
            numerator = sum(F((v*v-4*degree)*value)*q**degree
                            for degree, value in enumerate(factorial))
            assert kappa == lam*numerator/denominator
        assert row['groups']
        for group in row['groups']:
            zero, one1, one2, two, centered = [F(group[key]) for key in
                                                ['zero', 'one1', 'one2', 'two', 'centered']]
            assert zero > 0 and group['count'] == math.factorial(layout['d1'])*math.factorial(layout['d2'])
            assert one1 == kappas[0]*zero and one2 == kappas[1]*zero
            assert two == kappas[0]*kappas[1]*zero
            assert centered == two-kappas[0]*one2-kappas[1]*one1+kappas[0]*kappas[1]*zero == 0
        for example in row['negative_centered_examples']:
            value = q**example['q_crossings']*(lam*example['k1']-kappas[0])*(lam*example['k2']-kappas[1])
            assert F(example['centered_weight']) == value < 0
            negative_examples += 1
    assert negative_examples == 10
    return dict(counts=record['counts'], counts_reconstructed_from_saved_rows=actual,
                negative_centered_examples_verified=negative_examples,
                parameter_sha256=record['parameter_sha256'],
                exact_scope='Independent standard-library permutation and no-bridge enumeration; '
                            'saved polynomial identities and rational group sums rechecked here. '
                            'Signed individual terms are retained before group cancellation.',
                rational_scope='q=1/3 and lambda=1/4 are independent formal evaluation values, '
                               'not a physical pair satisfying q=exp(-4 lambda).',
                limitations='Fixed degrees through 7 and no-bridge layouts with 8 slots; '
                            'no arbitrary-degree or asymptotic theorem follows from these counts.')


def page_count(log):
    match = re.search(r'Output written on .*?\((\d+) pages?,', log, re.S)
    assert match, 'Missing successful PDF output line'
    return int(match.group(1))


def main():
    destination = OUT/'local-evidence.json'
    if destination.exists():
        raise FileExistsError(destination)
    parents = verify_parents()
    data = {name: load(OUT/name) for name in SCIENCE}
    for name, record in data.items():
        verify_inputs(record, name)
    checks = dict(descendants=verify_descendants(data['descendant-checks.json']),
                  source_regression=verify_source_regression(data['source-regression-checks.json']))

    source = ROOT/'manuscript'/STAGE/'main.tex'
    pdf = ROOT/'output/pdf'/f'de-sitter-{STAGE}.pdf'
    build = ROOT/'tmp/pdfs'/f'{STAGE}-build'
    final_log = text(OUT/'main.log')
    pages_count = page_count(final_log)
    assert pages_count > 0
    assert not re.search(r'undefined|Rerun to get|Overfull|Fatal error', final_log)
    assert sha(pdf) == sha(build/'main.pdf')
    assert sha(OUT/'main.log') == sha(build/'main.log')
    assert (build/'main.pdf').stat().st_mtime_ns >= source.stat().st_mtime_ns
    passes = [build/f'pass{index}.txt' for index in [1, 2]]
    pass_records = []
    for index, path in enumerate(passes, 1):
        output = text(path)
        assert page_count(output) == pages_count
        assert 'Fatal error' not in output and '! Emergency stop' not in output
        if index == 2:
            assert not re.search(r'undefined|Rerun to get|Overfull', output)
        pass_records.append(dict(pass_number=index, path=path.relative_to(ROOT).as_posix(),
                                 sha256=sha(path), compiler_header=output.splitlines()[0]))
    assert passes[1].stat().st_mtime_ns >= passes[0].stat().st_mtime_ns >= source.stat().st_mtime_ns
    width = len(str(pages_count))
    pages = [build/f'page-{index:0{width}d}.png' for index in range(1, pages_count+1)]
    assert all(path.is_file() for path in pages), 'Missing final Poppler render'
    assert all(path.stat().st_mtime_ns >= pdf.stat().st_mtime_ns for path in pages), 'Stale render'

    names = [
        'code/memory_wick_descendant_chords.py', 'code/check_memory_wick_descendants.py',
        'code/build_memory_wick_descendants.ps1', 'code/record_memory_wick_descendants.py',
        f'manuscript/{STAGE}/main.tex', f'output/pdf/de-sitter-{STAGE}.pdf',
        f'notes/{STAGE}.md', f'literature/{STAGE}-sources.md',
        'notes/susskind-memory-oriented-transfer-seed.md',
        f'journal/2026-11-07-{STAGE}.md', 'README.md', 'PROJECT.md', 'PLAN.md', 'NEXT.md']
    # Bind the second producer's driver(s) from its recorded inputs once available.
    names += [name for name in data['source-regression-checks.json']['input_sha256']
              if name.startswith('code/') and name not in names
              and name not in data['descendant-checks.json']['input_sha256']]
    diagnostics = sorted(path for path in (OUT/'diagnostics').rglob('*') if path.is_file())
    files = [ROOT/name for name in names]+[OUT/name for name in SCIENCE]+[OUT/'main.log']+diagnostics
    assert all(path.is_file() for path in files), [str(path) for path in files if not path.is_file()]
    assert len(files) == len(set(files))
    result = dict(
        status='Local research progress; long-term objective remains active',
        created_at_utc=datetime.now(timezone.utc).isoformat(),
        scientific_status='Ordered common-coupling Gaussian-Wick descendants have a constructive '
                          'fixed-degree moment algorithm and deterministic energy regressions. '
                          'The cut-diagram no-bridge cancellation and all-external-degree bound '
                          'give absolutely continuous residual two-energy measures; exact '
                          'samplewise pinching is then controlled. The established time-series '
                          'tails yield a convergent neutral finite-time covariance construction. '
                          'The charged defect, AA/AB/BA and bulk dictionary remain open.',
        parameters='Standard CAR, overlapping complex Gaussian p-body, common original J across '
                   'charges, normalized annealed beta=mu=0 full trace, hbar=1/string real time, '
                   'even N and lambda_N=p^2/N tending to positive lambda, fixed Js. Polynomial '
                   'degrees and external real time are fixed before the double-scale limit.',
        proof_scope='Finite checks verify complete Gaussian bookkeeping and exact low-degree '
                    'diagram formulas. They do not by themselves prove the arbitrary-degree '
                    'spectral bound or energy-pinching limit. Those use the manuscript cut-field '
                    'argument, Hermite bridge suppression, positive-measure convergence and '
                    'the predecessor uniform summable tails. The series tail bound is uniform '
                    'in compact memory times and all external real times; this is not a claim '
                    'of uniform-in-external-time convergence of the finite-N spectral measures.',
        continuation='Compute charged forward/reverse sources and the AA, AB and BA kernels '
                     'in the same original-J diagram system, recovering the previously nonzero '
                     'low-order mixed spectral calibration. The saved oriented-transfer seed '
                     'starts from the exact full-Fock D identity; its proposed weighted '
                     'whole-exponential limit has not been established in this stage. '
                     'Then return to the Susskind bulk '
                     'mass/position dictionary and interaction observables.',
        external_scope='Earlier Issue 1 review retained; no new issue reply, Git publication '
                       'or reintegration of the historical audit is made by this stage.',
        git_commit_created=False, git_push_performed=False, archive_created=False,
        external_reply_sent=False,
        previous_evidence=parents,
        input_hashes_verified={name: record['input_sha256'] for name, record in data.items()},
        parameter_hashes_verified={name: record['parameter_sha256'] for name, record in data.items()},
        checks=checks, evidence_scopes={name: record['scope'] for name, record in data.items()},
        recorded_driver_environments={name: record['environment'] for name, record in data.items()},
        environment=dict(python=platform.python_version(), numpy=np.__version__,
                         sympy=sp.__version__, mpmath=mp.__version__),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_wick_descendants.ps1',
                   pages=pages_count, passes=2, pass_records=pass_records,
                   compiler_header=final_log.splitlines()[0],
                   render_command=f'pdftoppm -r 95 -png output/pdf/de-sitter-{STAGE}.pdf '
                                  f'tmp/pdfs/{STAGE}-build/page',
                   render_number_width=width,
                   final_render_filenames=[path.name for path in pages],
                   final_pages_visually_inspected=list(range(1, pages_count+1)),
                   visual_result='Coordinating agent inspected every final render after the '
                                 'last build, including equations, references and page boundaries. '
                                 'Only the naming width selected by the final page count is bound.'),
        diagnostics=dict(preserved_files=[path.relative_to(ROOT).as_posix() for path in diagnostics],
                         scope='Earlier source/build/render diagnostics, if any, are preserved '
                               'separately from the validated final PDF and renders.'),
        process_state='Scientific checks, both final LaTeX passes and rendering are terminal '
                      'when this recorder is run. The recorder executes none of them.',
        hashes={path.relative_to(ROOT).as_posix(): sha(path) for path in files},
        render_hashes={path.relative_to(ROOT).as_posix(): sha(path) for path in pages},
        build_pass_hashes={path.relative_to(ROOT).as_posix(): sha(path) for path in passes})
    with destination.open('x', encoding='utf-8') as handle:
        handle.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(dict(status='passed', pages=pages_count, current_files=len(files),
                         scientific_outputs=len(SCIENCE),
                         previous_scientific_files_verified=sum(len(parent['unchanged_scientific_files'])
                                                                for parent in parents.values()),
                         final_renders=len(pages), final_build_passes=len(passes))))


if __name__ == '__main__':
    main()
