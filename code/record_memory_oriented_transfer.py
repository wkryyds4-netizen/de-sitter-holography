"""Seal oriented-transfer evidence after root inspects every final PDF page.

No scientific driver, compiler, renderer, historical main or Git is executed.
This file is deliberately a separate, exclusive-create final recording step.
"""
from datetime import datetime, timezone
from decimal import Decimal as Dec, localcontext
from fractions import Fraction as F
from math import comb
from pathlib import Path
import hashlib
import json
import platform
import re


ROOT = Path(__file__).resolve().parent.parent
STAGE = 'susskind-memory-oriented-transfer'
OUT = ROOT/'results'/STAGE
SCIENCE = ['finite-transfer-checks.json', 'oriented-mark-checks.json',
           'body-source-checks.json', 'charge-background-checks.json',
           'body-geometry-checks.json']
PARENTS = ['susskind-memory-wick-descendants', 'susskind-memory-uniform-bound',
           'susskind-memory-charged-transfer', 'susskind-memory-two-component']
MUTABLE_PLANS = {'README.md', 'PLAN.md', 'NEXT.md'}


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
    assert record['input_sha256'] and record['environment'] and record['scope']
    for filename, expected in record['input_sha256'].items():
        path = ROOT/filename
        assert path.is_file(), (name, 'missing input', filename)
        assert sha(path) == expected, (name, 'changed input', filename)
    assert digest(record['parameters']) == record['parameter_sha256'], (name, 'parameters')


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


def verify_transfer(record):
    systems = record['systems']
    assert [[row['N'], row['p'], row['maximum_degree']] for row in systems] == \
        record['parameters']['systems'] == [[2, 1, 4], [3, 1, 3], [3, 2, 3], [4, 2, 2]]
    expected = dict(matrix_identities=17465, integer_component_equalities=468408,
                    Hermite_transfer_coefficient_cases=4368, complete_OU_polynomial_cases=46)
    assert record['counts'] == expected
    setup_matrices = setup_components = 0
    transfer_coefficients = OU_cases = 0
    for system in systems:
        N, p, degree = system['N'], system['p'], system['maximum_degree']
        coordinates = comb(N, p)**2
        assert system['real_Gaussian_coordinates'] == system['full_basis_DH_equals_pH'] == coordinates
        assert system['full_Fock_dimension'] == 2**N
        assert system['all_sectors'] == list(range(N+1))
        assert system['neutral_matrix_units'] == comb(2*N, N)
        assert system['matrix_unit_sector_checks'] == N*comb(2*N, N)
        assert system['variance_check_all_sectors'] is True
        assert {(row['degree'], row['m']) for row in system['rows']} == {
            (j, m) for j in range(degree+1) for m in range(N+1)}
        setup_matrices += N*N+1+3*coordinates+(N+1)
        setup_components += (N*N+1+3*coordinates)*2*(2**N)**2+2*comb(2*N, N)
        for row in system['rows']:
            if row['m'] == 0:
                assert row['zero_charge_convention'] is True and row['E0_zero'] is True
                continue
            m = row['m']
            assert F(row['alpha']) == F(max(m-p, 0), m)
            assert row['transfer_and_number_identity'] is True
            assert row['conditional_OU_identity'] is True
            if row['degree'] <= 1 or m < p:
                assert row['nonzero_E_coefficients'] == 0
            transfer_coefficients += row['Hermite_coefficients']
            OU_cases += 1
    assert (transfer_coefficients, OU_cases) == (4368, 46)
    assert sum(row['counts']['matrix_identities'] for row in systems)+setup_matrices == 17465
    assert sum(row['counts']['integer_component_equalities'] for row in systems)+setup_components == 468408
    assert record['maximum_absolute_integer'] == 612 < record['pre_operation_integer_limit'] == 2**60
    taylor = record['taylor']
    rows = taylor['rows']
    assert digest(rows) == taylor['rows_sha256']
    actual = dict(cases=len(rows), endpoint_m_equals_p=sum(row['m'] == row['p'] for row in rows),
                  k_zero_to_three=sum(row['k'] <= 3 for row in rows),
                  zero_remainder=sum(F(row['remainder']) == 0 for row in rows),
                  upper_equality=sum(F(row['remainder']) == F(row['upper_bound']) for row in rows))
    assert actual == taylor['counts'] == dict(cases=12582, endpoint_m_equals_p=1026,
        k_zero_to_three=1272, zero_remainder=954, upper_equality=1272)
    assert all(0 <= F(row['remainder']) <= F(row['upper_bound']) for row in rows)
    return dict(counts=expected, taylor=actual, complete_neutral_matrix_units=116,
                setup_matrix_identities=setup_matrices, setup_integer_components=setup_components,
                scope='Complete small-system Gaussian coordinates and Hermite coefficients. '
                      'Variance 2/1 integer normalization is distinct from physical Js units.')


def verify_oriented_marks(record):
    counts = record['counts']
    geometry, finite = record['geometry_rows'], record['finite_car_rows']
    expected = dict(local_deficit_comparisons=381, local_words=127, local_nonzero_odd_words=6,
                    geometry_pairings=85, geometry_coefficient_comparisons=6120,
                    geometry_cases=765, negative_single_coefficients=312,
                    negative_connected_coefficients=184, legal_two_slot_cases=446,
                    finite_neutral_loss_comparisons=480, finite_words=80, finite_type_cases=720,
                    finite_trace_comparisons=5760, finite_negative_double_marks=4,
                    finite_zero_base_nonzero_mark=8, pair_site_table_rows=27,
                    p_minus_D_H_calibrations=45)
    assert counts == expected
    assert len(geometry) == 765 and len(finite) == 720
    assert len({digest({key: row[key] for key in ['length', 'blocks', 'pairing']})
                for row in geometry}) == 85
    assert sum(F(row['kappa'].get(str(key), '0')) < 0 for row in geometry for key in [1, 2]) == 312
    assert sum(F(row['kappa'].get('3', '0')) < 0 for row in geometry) == 184
    assert all(row['both_loss_forms_exact'] is True and row['geometry_exact'] is True for row in geometry)
    assert all(row['both_loss_forms_exact'] is True and len(row['all_mark_subsets']) == 4 for row in finite)
    assert sum(F(row['all_mark_subsets'][3]) < 0 for row in finite) == 4
    assert sum(F(row['all_mark_subsets'][0]) == 0 and any(F(x) for x in row['all_mark_subsets'][1:])
               for row in finite) == 8
    assert sum(row['legal_two_slot_cases'] for row in record['conditional_two_slot_rows']) == 446
    assert len(record['pair_site_table']) == 27 and len(record['calibration_rows']) == 45
    assert all(row['p_minus_D_H_exactly_zero'] is True for row in record['calibration_rows'])
    assert counts['local_words'] == sum(2**j for j in range(record['parameters']['local_word_max_length']+1))
    assert sum(len(row['word']) % 2 for row in record['local_word_rows']) == 6
    for row in record['pair_site_table']:
        assert F(row['connected']) == sum(F(item['kappa'].get('3', '0')) for item in row['signed_pair_site_terms'])
    diagnostic = record['raw_nonadjacent_shortcut_diagnostic']
    assert diagnostic['status'] == 'expected_failure_retained' and diagnostic['old_assertion_rejected'] is True
    assert diagnostic['direct_minor_alpha'] == -2 and diagnostic['invalid_shortcut_alpha'] == 2
    return dict(counts=counts, invalid_shortcut_preserved=True,
                scope='Coupling-Wick blocks without internal pairs; seeded selected finite words '
                      'are recorded explicitly and are not an exhaustive index-set enumeration.')


def verify_body_sources(record):
    finite, diagrams = record['finite'], record['no_bridge']
    actual = dict(finite_scalar_coefficients=sum(row['scalar_coefficients'] for row in finite),
                  no_bridge_diagrams=sum(len(row['diagrams']) for row in diagrams),
                  nonzero_no_bridge_connected=sum(row['nonzero_connected'] for row in diagrams))
    assert actual == record['counts'] == dict(finite_scalar_coefficients=6408,
        no_bridge_diagrams=294, nonzero_no_bridge_connected=80)
    assert len(finite) == 4 and len(diagrams) == 10
    for row in finite:
        assert row['exact'] is True
        assert row['scalar_coefficients'] == 4*row['pairings']*row['assignments']*row['layouts']
        assert 0 < row['negative_components'] < row['scalar_coefficients']
    for row in diagrams:
        assert row['cases'] == len(row['diagrams'])
        assert row['nonzero_connected'] == sum(item['internal'][0]*item['internal'][1] > 0
                                              for item in row['diagrams'])
        for item in row['diagrams']:
            assert F(item['source'].get('3', '0')) == 4*item['internal'][0]*item['internal'][1]
    assert record['maximum_absolute_integer'] == 4
    assert record['static_G2']['exponent'] == {'1': '-2', '2': '-2', '3': '4'}
    return dict(counts=actual,
                scope='Fixed CAR index assignments include internal contractions; the larger '
                      'cases use listed deterministic assignments. The nonzero static connected '
                      'terms are retained, not treated as a dynamic residual kernel.')


def verify_charge_background(record):
    exact, pair, numerical = record['exact'], record['pair_source'], record['high_precision']
    actual = dict(exact_systems=len(exact['rows']),
                  exact_sector_identities=sum(row['N']+1 for row in exact['rows']),
                  cross_collision_types=len(pair['cross_collision_ratios']),
                  own_collision_types=len(pair['own_collision_ratios']),
                  exponential_checks=len(pair['evaluated_exponentials']),
                  high_precision_systems=len(numerical['rows']),
                  central_cases=sum(len(row['central']) for row in numerical['rows']))
    assert actual == record['counts'] == dict(exact_systems=81, exact_sector_identities=1643,
        cross_collision_types=4, own_collision_types=2, exponential_checks=4,
        high_precision_systems=12, central_cases=72)
    assert exact['systems'] == 81 and exact['sector_identities'] == 1643
    assert numerical['systems'] == 12 and numerical['central_cases'] == 72
    for row in exact['rows']:
        assert F(row['lambda_N']) == F(row['p']**2, row['N'])
        assert all(row[key] is True for key in ['exact_difference_identity', 'exact_trace_shift_identity',
                                               'independent_support_overlap_identity'])
        assert F(row['variance_G0_over_C0_squared']) == \
            F(row['second_G0_over_C0_squared'])-F(row['mean_G0_over_C0'])**2 > 0
    assert pair['kappa'] == {'1': '-2', '2': '-2', '3': '4'}
    assert pair['cross_connected_exponent'] == {'3': '4'}
    for row in pair['evaluated_exponentials']:
        lam = F(row['lambda_value'])
        assert F(row['exponential']['1']) == F(row['exponential']['2']) == -2*lam
        assert F(row['exponential']['3']) == 4*lam+4*lam*lam
        assert F(row['connected_coefficient']) == 4*lam
    assert numerical['precisions'] == record['parameters']['precisions'] == [70, 110]
    with localcontext() as context:
        context.prec = 90
        for row in numerical['rows']:
            assert F(row['lambda_N']) == F(row['p']**2, row['N'])
            assert row['all_charge_sectors'] == row['N']+1
            assert row['identically_zero_low_sectors'] == row['p']
            assert Dec(row['maximum_precision_difference']) < Dec('1e-58')
            mean, second, variance = [Dec(row[key]) for key in
                ['mean_G0_over_C0', 'second_G0_over_C0_squared', 'variance_G0_over_C0_squared']]
            assert abs(second-mean*mean-variance) < Dec('1e-54')
            assert len(row['central']) == 6
            for point in row['central']:
                assert point['m'] == F(row['N'], 2)+F(point['x'])*row['p']
                assert abs(abs(Dec(point['observed'])-Dec(point['target']))-
                           Dec(point['absolute_error'])) < Dec('1e-54')
    return dict(counts=actual,
                maximum_precision_difference=max((row['maximum_precision_difference'] for row in numerical['rows']), key=Dec),
                scope='Every charge sector is included; 70/110-digit agreement is not interval '
                      'certification. The charge-static limit uses the separate analytic L2 proof.')


def verify_body_geometry(record):
    rows, calibration = record['geometry_rows'], record['G2_calibration']
    actual = dict(
        coefficient_comparisons=8*len(rows), geometry_cases=len(rows),
        cases_with_internal_pairs=sum(any(row['internal_pairs']) for row in rows),
        mixed_D_L_with_D_internal=sum(row['kinds'][0] != row['kinds'][1] and
                                     any(row['internal_pairs']) for row in rows),
        negative_linear_coefficients=sum(F(row['closed_kappa'][str(bit)]) < 0
                                         for row in rows for bit in [1, 2]),
        negative_connected_coefficients=sum(F(row['closed_kappa']['3']) < 0 for row in rows),
        negative_G_linear_coefficients=sum(F(row['G_kappa'][str(bit)]) < 0
                                           for row in rows for bit in [1, 2]),
        no_bridge_positive_charge_connected=sum(not row['shared_open_chords'] and
            F(row['closed_kappa']['3']) > 0 for row in rows),
        pure_Wick_cases=sum(row['pure_Wick_match_checked'] for row in rows),
        pure_Wick_comparisons=4*sum(row['pure_Wick_match_checked'] for row in rows),
        G2_diagrams=len(calibration['diagrams']),
        G2_formal_rational_checks=len(calibration['formal_rational_evaluations']))
    actual.update({'cases_'+kind: sum(''.join(row['kinds']) == kind for row in rows)
                   for kind in ['DD', 'DL', 'LD', 'LL']})
    expected = dict(coefficient_comparisons=18944, geometry_cases=2368, cases_DD=858,
                    cases_with_internal_pairs=724, mixed_D_L_with_D_internal=277,
                    negative_linear_coefficients=1454, negative_connected_coefficients=178,
                    negative_G_linear_coefficients=1946, no_bridge_positive_charge_connected=80,
                    cases_DL=572, cases_LD=527, pure_Wick_comparisons=6576,
                    pure_Wick_cases=1644, cases_LL=411, G2_diagrams=3,
                    G2_formal_rational_checks=3)
    assert actual == record['counts'] == expected
    for row in rows:
        assert row['closed_geometry_exact'] is True and row['G_shift_exact'] is True
        assert row['local_kappa'] == row['closed_kappa']
        assert row['pure_Wick_match_checked'] == (not any(row['internal_pairs']))
        assert row['G_kappa']['0'] == row['closed_kappa']['0']
        assert row['G_kappa']['3'] == row['closed_kappa']['3']
        for index, kind in enumerate(row['kinds']):
            if kind == 'L':
                assert not row['internal_pairs'][index]
            degree = len(row['open_chords'][index])
            shift = degree*(degree-1) if kind == 'D' else 0
            assert F(row['G_kappa'][str(1 << index)]) == F(row['closed_kappa'][str(1 << index)])-shift
    assert len(record['parents_verified']) == 3
    for filename, parent in record['parents_verified'].items():
        assert parent['input_hashes_verified'] is True
        assert sha(ROOT/filename) == parent['sha256'] == record['input_sha256'][filename]
    assert calibration['old_w_recovered_after_charge_subtraction'] is True
    assert calibration['energy_only_equals_old_w'] is False
    combined = {}
    for diagram in calibration['diagrams']:
        for term in diagram['norm_contribution_over_C_squared']:
            key = (term['lambda_power'], term['q_power'])
            combined[key] = combined.get(key, F())+F(term['coefficient'])
    assert combined == {(term['lambda_power'], term['q_power']): F(term['coefficient'])
                        for term in calibration['raw_norm_polynomial']}
    for row in calibration['formal_rational_evaluations']:
        lam, q = F(row['lambda']), F(row['q'])
        raw = F(row['raw_norm_over_C_squared'])
        energy = F(row['energy_only_residual_over_C_squared'])
        joint = F(row['charge_and_energy_residual_over_C_squared'])
        assert raw == 8*lam+4*lam*lam+16*lam*lam*q
        assert energy == 8*lam+16*lam*lam*q/(1+q)
        assert energy-joint == 4*lam and joint == F(row['old_w_over_C_squared'])
    return dict(counts=actual, same_stage_input_links=record['parents_verified'],
                G2_energy_only_equals_old_w=False, G2_joint_subtraction_recovers_old_w=True,
                scope='Signed finite diagrams compare whole-field D with Wick-restricted L. '
                      'Formal rational lambda/q checks are independent algebraic substitutions, '
                      'not physical q=exp(-4 lambda) parameter pairs. The G2 subtraction is a '
                      'low-degree calibration and not a general pinching theorem.')


def page_count(log):
    match = re.search(r'Output written on .*?\((\d+) pages?,', log, re.S)
    assert match, 'Missing successful PDF output line'
    return int(match.group(1))


def main():
    destination = OUT/'local-evidence.json'
    if destination.exists():
        raise FileExistsError(destination)
    parents = verify_parents()
    science = {name: load(OUT/name) for name in SCIENCE}
    for name, record in science.items():
        verify_inputs(record, name)
    verifiers = [verify_transfer, verify_oriented_marks, verify_body_sources,
                 verify_charge_background, verify_body_geometry]
    checks = {name: function(science[name]) for name, function in zip(SCIENCE, verifiers)}
    assert science['body-source-checks.json']['static_G2']['exponent'] == \
        science['charge-background-checks.json']['pair_source']['kappa']

    source = ROOT/'manuscript'/STAGE/'main.tex'
    pdf = ROOT/'output/pdf'/f'de-sitter-{STAGE}.pdf'
    build = ROOT/'tmp/pdfs'/f'{STAGE}-build'
    render = ROOT/'tmp/pdfs'/f'{STAGE}-render'
    final_log = text(OUT/'main.log')
    number = page_count(final_log)
    assert number > 0
    assert not re.search(r'undefined|Rerun to get|Overfull|Fatal error', final_log)
    assert sha(pdf) == sha(build/'main.pdf')
    assert sha(OUT/'main.log') == sha(build/'main.log')
    assert (build/'main.pdf').stat().st_mtime_ns >= source.stat().st_mtime_ns
    passes = [build/f'pass{i}.txt' for i in [1, 2]]
    pass_records = []
    for i, path in enumerate(passes, 1):
        output = text(path)
        assert page_count(output) == number
        assert 'Fatal error' not in output and '! Emergency stop' not in output
        if i == 2:
            assert not re.search(r'undefined|Rerun to get|Overfull', output)
        pass_records.append(dict(pass_number=i, path=path.relative_to(ROOT).as_posix(),
                                 sha256=sha(path), compiler_header=output.splitlines()[0]))
    assert passes[1].stat().st_mtime_ns >= passes[0].stat().st_mtime_ns >= source.stat().st_mtime_ns
    width = len(str(number))
    pages = [render/f'page-{i:0{width}d}.png' for i in range(1, number+1)]
    assert set(render.glob('page-*.png')) == set(pages), 'Final render names/count differ from final log'
    assert all(path.stat().st_mtime_ns >= pdf.stat().st_mtime_ns for path in pages), 'Stale final render'

    names = {
        'code/build_memory_oriented_transfer.ps1', 'code/record_memory_oriented_transfer.py',
        f'manuscript/{STAGE}/main.tex', f'output/pdf/de-sitter-{STAGE}.pdf',
        f'notes/{STAGE}.md', f'notes/{STAGE}-proof-draft.md', f'notes/{STAGE}-seed.md',
        'notes/susskind-memory-charge-background.md', f'literature/{STAGE}-sources.md',
        f'journal/2026-11-07-{STAGE}.md', 'README.md', 'PROJECT.md', 'PLAN.md', 'NEXT.md'}
    for record in science.values():
        names.update(record['input_sha256'])
    names.update(f'results/{STAGE}/{name}' for name in SCIENCE)
    names.add(f'results/{STAGE}/main.log')
    diagnostics = sorted(path for path in (OUT/'diagnostics').rglob('*') if path.is_file())
    files = [ROOT/name for name in sorted(names)]+diagnostics
    assert all(path.is_file() for path in files), [str(path) for path in files if not path.is_file()]
    assert len(files) == len(set(files))
    result = dict(
        status='Local research progress; long-term objective remains active',
        created_at_utc=datetime.now(timezone.utc).isoformat(),
        scientific_status='Exact full-Fock adjacent-charge transfer is represented by the body-loss '
                          'operator D. Weighted compact-time replacement by the whole-field G source '
                          'and summable ordinary-power/force tails are established in the manuscript. '
                          'Nilpotent source normalization retains internal Gaussian contractions '
                          'and identifies the charge-static background already present in G2. '
                          'Full joint charge-energy regression and charged AA/AB/BA spectra remain open.',
        parameters='Standard CAR, overlapping complex Gaussian p-body model, same original J '
                   'across charges, annealed normalized beta=mu=0 full trace, hbar=1/string time; '
                   'even N, p to infinity, lambda_N=p^2/N tending to positive lambda, fixed Js. '
                   'Memory times are in fixed compact intervals; polynomial degrees are fixed '
                   'before the double-scale limit. Finite integer checks use a separate 2/1 variance.',
        proof_scope='Saved finite identities, selected CAR words, signed local graph coefficients '
                    'and full-binomial numerical sums check the recorded assertions at their '
                    'specified parameters. They do not by themselves prove weighted all-time '
                    'bounds, arbitrary-degree Poisson geometry or a static-subtraction theorem. '
                    'Those manuscript arguments and the unresolved joint charge-energy regression '
                    'are distinguished explicitly. No-Gaussian-bridge does not remove charge-static '
                    'correlations in whole-field oriented sources.',
        continuation='Construct the joint charge-energy regression of G_n and K_n, recover the '
                     'known low-order A0 and nonzero A0/B2 mixed calibration after proper static '
                     'subtraction, and determine AA/AB/BA. Then combine with the neutral kernel '
                     'before returning to the Susskind bulk dictionary and interactions.',
        external_scope='Earlier Issue 1 audit retained; no new issue reply or historical reintegration.',
        git_commit_created=False, git_push_performed=False, archive_created=False, external_reply_sent=False,
        previous_evidence=parents, checks=checks,
        input_hashes_verified={name: record['input_sha256'] for name, record in science.items()},
        parameter_hashes_verified={name: record['parameter_sha256'] for name, record in science.items()},
        evidence_scopes={name: record['scope'] for name, record in science.items()},
        recorded_driver_environments={name: record['environment'] for name, record in science.items()},
        environment=dict(python=platform.python_version(), platform=platform.platform(),
                         recorder_dependencies='Python standard library only'),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_oriented_transfer.ps1',
                   pages=number, passes=2, pass_records=pass_records,
                   compiler_header=final_log.splitlines()[0],
                   render_command=f'pdftoppm -r 95 -png output/pdf/de-sitter-{STAGE}.pdf '
                                  f'tmp/pdfs/{STAGE}-render/page',
                   render_number_width=width, final_render_filenames=[path.name for path in pages],
                   final_pages_visually_inspected=list(range(1, number+1)),
                   visual_result='Root inspected every final 95dpi render after the last build, '
                                 'including formula line lengths, reference text and page boundaries.'),
        diagnostics=dict(preserved_files=[path.relative_to(ROOT).as_posix() for path in diagnostics],
                         scope='Earlier failures/layouts remain separate from the verified final build.'),
        process_state='The five scientific checks, both final LaTeX passes and final rendering '
                      'are terminal when this recorder is run; it executes none of them.',
        hashes={path.relative_to(ROOT).as_posix(): sha(path) for path in files},
        render_hashes={path.relative_to(ROOT).as_posix(): sha(path) for path in pages},
        build_pass_hashes={path.relative_to(ROOT).as_posix(): sha(path) for path in passes})
    with destination.open('x', encoding='utf-8') as handle:
        handle.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(dict(status='passed', pages=number, scientific_outputs=len(SCIENCE),
                         current_files=len(files), final_renders=len(pages), final_build_passes=2,
                         previous_scientific_files_verified=sum(len(parent['unchanged_scientific_files'])
                                                                for parent in parents.values()))))


if __name__ == '__main__':
    main()
