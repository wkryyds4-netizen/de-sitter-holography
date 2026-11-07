"""Exclusive final evidence seal, run only after final PDF visual inspection.

Reads the three saved scientific outputs and their inputs. Does not run a
scientific driver, compiler, renderer, historical main or Git command.
"""
from collections import Counter
from datetime import datetime, timezone
from decimal import Decimal as Dec, localcontext
from fractions import Fraction as F
from itertools import product
from math import comb, isfinite
from pathlib import Path
import hashlib
import json
import platform
import re


ROOT = Path(__file__).resolve().parent.parent
STAGE = 'susskind-memory-joint-regression'
OUT = ROOT/'results'/STAGE
SCIENCE = ['joint-regression-checks.json', 'charge-tilt-checks.json', 'joint-series-checks.json']
PARENTS = ['susskind-memory-oriented-transfer', 'susskind-memory-wick-descendants',
           'susskind-memory-two-component', 'susskind-memory-uniform-bound']
MUTABLE_PLANS = {'README.md', 'PLAN.md', 'NEXT.md'}
REQUIRED_DOCS = [f'notes/{STAGE}.md', f'notes/{STAGE}-proof.md',
                 'notes/susskind-memory-charge-tilt-proof.md',
                 f'literature/{STAGE}-sources.md', 'literature/susskind-memory-to-bilocal-sources.md',
                 f'journal/2026-11-07-{STAGE}.md', 'README.md', 'PROJECT.md', 'PLAN.md', 'NEXT.md']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text(encoding='utf-8'))


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def text(path):
    raw = path.read_bytes()
    encoding = 'utf-16' if raw.startswith((b'\xff\xfe', b'\xfe\xff')) else 'utf-8-sig'
    return raw.decode(encoding, errors='replace')


def verify_inputs(record, name):
    assert record['schema_version'] == 1 and record['status'] == 'passed', name
    assert record['environment'] and record['scope'] and record['input_sha256']
    assert digest(record['parameters']) == record['parameter_sha256'], (name, 'parameters')
    for filename, expected in record['input_sha256'].items():
        assert (ROOT/filename).is_file(), (name, 'missing input', filename)
        assert sha(ROOT/filename) == expected, (name, 'changed input', filename)
    for filename, parent in record.get('parents_verified', {}).items():
        assert parent['input_hashes_verified'] is True
        assert parent['sha256'] == record['input_sha256'][filename] == sha(ROOT/filename)
    for filename, expected in record.get('parents_sha256', {}).items():
        assert expected == record['input_sha256'][filename] == sha(ROOT/filename)


def verify_parents():
    records = {}
    for stage in PARENTS:
        path = ROOT/'results'/stage/'local-evidence.json'
        parent = load(path)
        checked = []
        for filename, expected in parent['hashes'].items():
            if filename in MUTABLE_PLANS:
                continue
            assert (ROOT/filename).is_file(), (stage, 'missing predecessor', filename)
            assert sha(ROOT/filename) == expected, (stage, 'changed predecessor', filename)
            checked.append(filename)
        records[path.relative_to(ROOT).as_posix()] = dict(
            sha256=sha(path), unchanged_scientific_files=checked,
            excluded_mutable_plans=sorted(MUTABLE_PLANS.intersection(parent['hashes'])))
    return records


def verify_regression(record):
    expected = dict(joint_projection_tests=68, nonzero_joint_projection_tests=68,
                    residual_Gram_layouts=39, residual_Gram_equalities=117,
                    nonzero_residual_Gram_layouts=39, no_bridge_layouts=30,
                    nonzero_raw_no_bridge_layouts=28, mixed_neutral_layouts=44,
                    nonzero_mixed_neutral_layouts=32, neutral_algorithm_crosschecks=16,
                    nonzero_A0B2_Mehler_checks=18, zero_A0B1_checks=12,
                    charge_static_diagnostic=2, G2_joint_norm_calibration=2)
    total, statistics = Counter(), Counter()
    assert [model['parameters'] for model in record['models']] == record['parameters']['models']
    assert record['parameters']['C'] == 1 and record['parameters']['maximum_raw_slots'] == 10
    for model in record['models']:
        projection, grams, b0 = model['joint_projection'], model['residual_Grams'], model['no_bridge']
        mixed, neutral = model['mixed_neutral'], model['neutral_algorithm_crosschecks']
        calibration = model['low_order_spectral_calibration']
        actual = dict(joint_projection_tests=len(projection),
            nonzero_joint_projection_tests=sum(F(row['raw_joint_moment']) != 0 for row in projection),
            residual_Gram_layouts=len(grams), residual_Gram_equalities=3*len(grams),
            nonzero_residual_Gram_layouts=sum(F(row['RR']) != 0 for row in grams),
            no_bridge_layouts=len(b0), nonzero_raw_no_bridge_layouts=sum(F(row['GG']) != 0 for row in b0),
            mixed_neutral_layouts=len(mixed),
            nonzero_mixed_neutral_layouts=sum(F(row['whole_residual_V']) != 0 for row in mixed),
            neutral_algorithm_crosschecks=len(neutral),
            nonzero_A0B2_Mehler_checks=sum(row['channel'] == 'A0B2' for row in calibration),
            zero_A0B1_checks=sum(row['channel'] == 'A0B1' for row in calibration),
            charge_static_diagnostic=1, G2_joint_norm_calibration=1)
        assert actual == model['counts']
        assert all(row['exact'] is True for collection in [projection, grams, b0, mixed, neutral, calibration]
                   for row in collection)
        for row in projection:
            assert F(row['raw_joint_moment']) == F(row['regression_joint_moment'])
        for row in grams:
            assert F(row['RR']) == F(row['RS']) == F(row['SR']) == F(row['SS'])
        for row in b0:
            assert F(row['GG'])-F(row['gG'])-F(row['Gg'])+F(row['gg']) == F(row['centered_group_sum']) == 0
        for row in mixed:
            assert F(row['whole_residual_V']) == F(row['resummed_J_V'])
        for row in neutral:
            assert F(row['sealed_algorithm']) == F(row['new_algorithm'])
        lam, q = F(model['parameters']['lambda']), F(model['parameters']['q'])
        w = 4*lam+16*lam*lam*q/(1+q)
        eta = -4*lam*q+16*lam*lam*q/(1+q)
        for row in calibration:
            if row['channel'] == 'A0B2':
                assert F(row['known_prefactor']) == 4*lam*q*q*(w+eta)
                assert F(row['direct']) == F(row['reference']) == \
                    F(row['known_prefactor'])*F(row['Mehler_moment_q_squared']) > 0
            else:
                assert row['channel'] == 'A0B1' and F(row['direct']) == 0
        diagnostic = model['charge_static_diagnostic']
        assert F(diagnostic['energy_only_b0_value']) == F(diagnostic['known_charge_variance']) == 4*lam
        assert F(diagnostic['joint_centered_b0_value']) == 0
        assert F(diagnostic['joint_G2_norm']) == F(diagnostic['old_w']) == w
        assert len(model['candidate_fields']) == model['parameters']['maximum_G_degree']-1
        for field in model['candidate_fields']:
            assert all(term['charge_power'] in [0, 1] for term in field['g_n'])
            assert all(term['charge_power'] == 0 for term in field['resummed_J'])
        total.update(actual)
        statistics.update(model['diagram_statistics'])
    assert dict(total) == record['counts'] == expected
    assert dict(statistics) == record['diagram_statistics'] == dict(
        accepted_term_diagrams=93277, nonzero_charge_connected_diagrams=5203,
        negative_diagram_weights=34370, term_moment_cache_misses=2305)
    return dict(counts=expected, signed_diagram_statistics=dict(statistics),
                scope='Exact formal rational common-J moments at recorded degrees. The '
                      'nonzero A0/B2 calibration and charge-only static remainder are retained. '
                      'These tests do not themselves prove arbitrary-degree regression or pinching.')


def dc(value):
    return (Dec(value['real']), Dec(value['imag']))


def cadd(left, right):
    return (left[0]+right[0], left[1]+right[1])


def csub(left, right):
    return (left[0]-right[0], left[1]-right[1])


def cmul(left, right):
    return (left[0]*right[0]-left[1]*right[1], left[0]*right[1]+left[1]*right[0])


def norm2(value):
    return value[0]**2+value[1]**2


def fraction_decimal(value):
    value = F(value)
    return Dec(value.numerator)/Dec(value.denominator)


def verify_tilt(record):
    finite, singleton = record['finite'], record['singleton']
    scalar, background = record['scalar'], record['G2_background']
    assert [(row['N'], row['p']) for row in finite['systems']] == [(2, 1), (3, 1), (4, 2)]
    actual = {key:sum(row['counts'][key] for row in finite['systems']) for key in finite['counts']}
    assert actual == finite['counts'] == dict(support_pairing_layout_cases=19502,
        charge_source_coefficients=315192, nonzero_coefficients=21258, negative_coefficients=6929)
    for row in finite['systems']:
        assignments = comb(row['N'], row['p'])**4
        assert row['assignments'] == assignments
        assert row['counts']['support_pairing_layout_cases'] == 14*assignments
        assert row['counts']['charge_source_coefficients'] == 46*(row['N']+1)*assignments
        assert row['all_charge_and_source_coefficients_equal'] is True
    assert len(finite['examples']) == 18
    for example in finite['examples']:
        N, p = example['N'], example['p']
        for row in example['moments']:
            coefficients = list(map(F, example['charge_coefficients'][str(row['mask'])]))
            assert len(coefficients) == N+1
            for r, expected in enumerate(row['X_moments']):
                assert sum(coefficients[m]*F(2*m-N, 2*p)**r for m in range(N+1))/2**N == F(expected)
    assert len(singleton['singleton_rows']) == singleton['singleton_cases'] == 16
    assert 4*len(singleton['collision_rows']) == singleton['collision_ratios'] == 16
    for row in singleton['singleton_rows']:
        rho, p = F(row['rho']), row['p']
        ratio = rho/(1-rho)
        assert F(row['exp_t_over_p']) == ratio
        assert F(row['one'].get('1', '0')) == p*(ratio-1)
        assert F(row['two'].get('3', '0')) == p*p*(ratio-1)**2
    for row in singleton['collision_rows']:
        rho, ratios = F(row['rho']), row['ratios']
        assert F(ratios['n,n']['0']) == 1/rho
        assert F(ratios['n,bar']['2']) == F(ratios['bar,n']['1']) == 1/(1-rho)
        assert F(ratios['bar,bar']['3']) == rho/(1-rho)**3
    assert len(scalar['rows']) == scalar['cases'] == 216
    assert len(background['rows']) == background['cases'] == 20
    assert scalar['precisions'] == background['precisions'] == [70, 110]
    with localcontext() as context:
        context.prec = 100
        tolerance = Dec('1e-52')
        for row in scalar['rows']:
            assert row['all_three_bounds'] is True
            assert Dec(row['maximum_precision_difference']) < Dec('1e-58')
            values = row['values']
            t = tuple(fraction_decimal(x) for x in row['t'])
            assert norm2(csub(csub(dc(values['source_coefficient']), t), dc(values['source_error']))) < tolerance**2
            for error, bound in [('source_error', 'source_bound'), ('singleton_log', 'singleton_log_bound'),
                                  ('charge_log_error', 'charge_log_bound')]:
                assert norm2(dc(values[error])) <= Dec(values[bound]['real'])**2+tolerance
                assert Dec(values[bound]['imag']) == 0 and Dec(values[bound]['real']) >= 0
        for row in background['rows']:
            assert F(row['lambda_N']) == F(row['p']**2, row['N'])
            assert row['all_charge_sectors'] == row['N']+1 and row['exactly_zero_low_sectors'] == row['p']
            assert Dec(row['maximum_precision_difference']) < Dec('1e-57')
            lam = fraction_decimal(row['lambda_N'])
            t = tuple(fraction_decimal(x) for x in row['t'])
            shifted = cadd(t, (-2*lam, Dec(0)))
            values = row['values']
            M = dc(values['target_charge_MGF'])
            assert norm2(csub(dc(values['target_mean']), cmul(shifted, M))) < tolerance**2
            second = cmul(cadd(cmul(shifted, shifted), (4*lam, Dec(0))), M)
            assert norm2(csub(dc(values['target_second']), second)) < tolerance**2
            for key in ['mean', 'second']:
                difference = norm2(csub(dc(values[key]), dc(values['target_'+key])))
                assert abs(difference-Dec(row[key+'_absolute_error'])**2) < Dec('1e-36')
    counts = dict(**actual, exact_singletons=16, exact_collision_ratios=16,
                  scalar_complex_tilt_cases=216, G2_full_charge_tilt_cases=20)
    assert counts == record['counts']
    assert record['maximum_absolute_integer'] == 256 < record['integer_operation_guard'] == 2**60
    return dict(counts=counts, recorded_exact_charge_moment_examples=18,
                scope='All listed support assignments and charge/source coefficients are covered. '
                      'Unmarked singleton finite corrections, negative coefficients and actual '
                      'finite-N tilt errors remain. Complex precision comparisons are not interval '
                      'certificates; uniform-integrability and Cauchy arguments are in the proof note.')


def verify_series(record):
    parameters, exact, numerical = record['parameters'], record['exact'], record['numerical']
    actual = dict(ordinary_Hermite_identities=sum(len(row['a_rows']) for row in exact['rows']),
                  centered_force_coefficient_identities=sum(len(row['force_checks'])*(parameters['force_degree']+1)
                                                           for row in exact['rows']),
                  quadrature_coefficients=len(numerical['rows']))
    assert actual == record['counts'] == dict(ordinary_Hermite_identities=90,
        centered_force_coefficient_identities=1782, quadrature_coefficients=405)
    assert exact['ordinary_Hermite_identities'] == 90 and exact['centered_force_coefficient_identities'] == 1782
    assert {(row['q'], row['C']) for row in exact['rows']} == set(product(parameters['rational_q'], parameters['rational_C']))
    for row in exact['rows']:
        q, C = F(row['q']), F(row['C'])
        assert row['degree_identities'] == len(row['a_rows']) == parameters['maximum_degree']+1 == 15
        for n, a in enumerate(row['a_rows']):
            assert F(a[str(n)]) == 1
            assert all((n-int(d)) % 2 == 0 and F(value) > 0 for d, value in a.items())
        for entry in row['force_checks']:
            d = entry['d']
            if d < 2:
                continue
            x, y = map(F, entry['energies'])
            e = (x+y)/2
            expected = [F(d == 2), F(d == 3)-3*e*(d == 2),
                        F(d == 4)-4*e*(d == 3)+(C*(3+2*q+q*q)+x*y+5*e*e)*(d == 2)]
            assert list(map(F, entry['centered_first_three'])) == expected
    assert parameters['quadrature_sizes'] == [48, 72]
    assert numerical['cases'] == 405
    keys = {(row['lambda_value'], row['time'], row['degree'], row['derivative']) for row in numerical['rows']}
    assert keys == set(product(parameters['physical_lambda'], parameters['times'],
                               parameters['spectral_degrees'], parameters['derivative_orders']))
    maximum = max(row['absolute_difference'] for row in numerical['rows'])
    assert maximum == numerical['maximum_absolute_difference'] == 1.2800225526324438e-14
    assert maximum < parameters['floating_comparison_tolerance'] == 2e-10
    for row in numerical['rows']:
        assert all(isfinite(value) for key in ['value_48', 'value_72'] for value in row[key])
        assert 0 <= row['absolute_difference'] <= maximum
        assert Dec(row['analytic_Taylor_quadrature_tail_48']) >= 0
    return dict(counts=actual, maximum_quadrature_difference=maximum,
                scope='Formal rational coefficient identities and physical q=exp(-4 lambda) '
                      'limiting coefficient quadratures are distinct. The analytic tail covers '
                      'Taylor/quadrature truncation, not floating roundoff or finite-N dynamics.')


def page_count(log):
    match = re.search(r'Output written on .*?\((\d+) pages?,', log, re.S)
    assert match, 'Missing successful PDF output line'
    return int(match.group(1))


def main():
    destination = OUT/'local-evidence.json'
    if destination.exists():
        raise FileExistsError(destination)
    parents = verify_parents()
    science = {name:load(OUT/name) for name in SCIENCE}
    for name, record in science.items():
        verify_inputs(record, name)
    checks = {name:function(science[name]) for name,function in zip(
        SCIENCE, [verify_regression, verify_tilt, verify_series])}

    source = ROOT/'manuscript'/STAGE/'main.tex'
    pdf = ROOT/'output/pdf'/f'de-sitter-{STAGE}.pdf'
    build = ROOT/'tmp/pdfs'/f'{STAGE}-build'
    render = ROOT/'tmp/pdfs'/f'{STAGE}-render'
    log = text(OUT/'main.log')
    pages_count = page_count(log)
    assert pages_count > 0 and not re.search(r'undefined|Rerun to get|Overfull|Fatal error', log)
    assert sha(pdf) == sha(build/'main.pdf')
    assert sha(OUT/'main.log') == sha(build/'main.log')
    assert (build/'main.pdf').stat().st_mtime_ns >= source.stat().st_mtime_ns
    passes = [build/f'pass{i}.txt' for i in [1,2]]
    pass_records = []
    for i, path in enumerate(passes,1):
        output = text(path)
        assert page_count(output) == pages_count
        assert 'Fatal error' not in output and '! Emergency stop' not in output
        if i == 2:
            assert not re.search(r'undefined|Rerun to get|Overfull', output)
        pass_records.append(dict(pass_number=i,path=path.relative_to(ROOT).as_posix(),
                                 sha256=sha(path),compiler_header=output.splitlines()[0]))
    assert passes[1].stat().st_mtime_ns >= passes[0].stat().st_mtime_ns >= source.stat().st_mtime_ns
    width = len(str(pages_count))
    pages = [render/f'page-{i:0{width}d}.png' for i in range(1,pages_count+1)]
    assert set(render.glob('page-*.png')) == set(pages), 'Final render filenames/count mismatch'
    assert all(path.stat().st_mtime_ns >= pdf.stat().st_mtime_ns for path in pages), 'Stale final render'

    names = set(REQUIRED_DOCS) | {
        'code/build_memory_joint_regression.ps1', 'code/record_memory_joint_regression.py',
        f'manuscript/{STAGE}/main.tex', f'output/pdf/de-sitter-{STAGE}.pdf', f'results/{STAGE}/main.log'}
    for record in science.values():
        names.update(record['input_sha256'])
    names.update(f'results/{STAGE}/{name}' for name in SCIENCE)
    diagnostics = sorted(path for path in (OUT/'diagnostics').rglob('*') if path.is_file())
    files = [ROOT/name for name in sorted(names)]+diagnostics
    assert all(path.is_file() for path in files), [str(path) for path in files if not path.is_file()]
    assert len(files) == len(set(files))
    result = dict(
        status='Local research milestone; the long-term de Sitter objective remains active',
        created_at_utc=datetime.now(timezone.utc).isoformat(),
        scientific_status='The manuscript constructs the complete charged/neutral common-J '
                          'finite-time memory covariance AA+AB+BA+BB in the stated model and limit. '
                          'Complex charge tilt identifies the joint charge-energy regression; '
                          'same-sample mixed Grams give strong residual reduction, fixed-field '
                          'spectral absolute continuity controls exact energy pinching, and '
                          'summable time tails give the convergent full kernel. The bulk de Sitter '
                          'dictionary remains the next physical problem.',
        parameters='Standard CAR, overlapping complex Gaussian p-body model, same J across '
                   'charge sectors, normalized annealed beta=mu=0 full-Fock trace, hbar=1/string '
                   'time; even N, p to infinity, lambda_N=p^2/N to positive lambda, fixed Js. '
                   'Polynomial degrees and real s,u,T are fixed before the double-scale limit.',
        proof_scope='Finite exact graph/CAR checks, complex-source comparisons and coefficient '
                    'quadratures validate their recorded interfaces. They do not substitute for '
                    'the manuscript arbitrary-degree joint regression, positive-measure limit, '
                    'pinching and time-resummation proofs. The infinite spectral density is L1; '
                    'fixed-field L2 density bounds are not claimed uniform in field degree. '
                    'Uniformity for all external real T applies to series tails, not finite-N '
                    'spectral convergence. No growing-time or simultaneous small-lambda limit is taken.',
        continuation='Use the constructed spectral matrix to test the Susskind de Sitter '
                     'mass/position and bilocal dictionary, while checking the remaining A1 naming '
                     'calibration. Finite-temperature, uniform small-lambda and finite-G four-point '
                     'interpretations are separate questions, not consequences of this covariance proof.',
        external_scope='Earlier Issue 1 review is retained; this stage makes no new issue reply or historical reintegration.',
        git_commit_created=False,git_push_performed=False,archive_created=False,external_reply_sent=False,
        previous_evidence=parents,checks=checks,
        input_hashes_verified={name:record['input_sha256'] for name,record in science.items()},
        parameter_hashes_verified={name:record['parameter_sha256'] for name,record in science.items()},
        evidence_scopes={name:record['scope'] for name,record in science.items()},
        recorded_driver_environments={name:record['environment'] for name,record in science.items()},
        environment=dict(python=platform.python_version(),platform=platform.platform(),
                         recorder_dependencies='Python standard library only'),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_joint_regression.ps1',
                   pages=pages_count,passes=2,pass_records=pass_records,compiler_header=log.splitlines()[0],
                   render_command=f'pdftoppm -r 95 -png output/pdf/de-sitter-{STAGE}.pdf tmp/pdfs/{STAGE}-render/page',
                   render_number_width=width,final_render_filenames=[path.name for path in pages],
                   final_pages_visually_inspected=list(range(1,pages_count+1)),
                   visual_result='Root inspected every final 95dpi page after the last source/build '
                                 'revision; equations, page boundaries and references show no clipping or overlap.'),
        diagnostics=dict(preserved_files=[path.relative_to(ROOT).as_posix() for path in diagnostics],
                         scope='Earlier diagnostics remain separate from the verified final artifacts.'),
        process_state='All three science outputs, both final LaTeX passes and final rendering are '
                      'terminal when this recorder runs. It does not execute any of those processes.',
        hashes={path.relative_to(ROOT).as_posix():sha(path) for path in files},
        render_hashes={path.relative_to(ROOT).as_posix():sha(path) for path in pages},
        build_pass_hashes={path.relative_to(ROOT).as_posix():sha(path) for path in passes})
    with destination.open('x',encoding='utf-8') as handle:
        handle.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(status='passed',pages=pages_count,scientific_outputs=len(SCIENCE),
        current_files=len(files),final_renders=len(pages),final_build_passes=2,
        previous_scientific_files_verified=sum(len(parent['unchanged_scientific_files']) for parent in parents.values()))))


if __name__ == '__main__':
    main()
