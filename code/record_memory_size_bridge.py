"""Seal the size-bridge stage after the final ten pages have been inspected.

This recorder never runs scientific drivers, LaTeX, plotting or rendering.
Run only after the final manuscript, two build passes and zero-padded renders
have been completed and visually inspected by the coordinating agent.
"""
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path
import hashlib
import json
import platform
import re

import matplotlib
import mpmath as mp
import numpy as np
import sympy as sp


ROOT = Path(__file__).resolve().parent.parent
STAGE = 'susskind-memory-size-bridge'
OUT = ROOT / 'results' / STAGE
EXPECTED_PAGES = 10
SCIENCE = ['finite-bridge.json', 'marked-chord-checks.json',
           'geometry-checks.json', 'rank-law-checks.json']
PARENTS = ['susskind-memory-neutral-marked-series',
           'susskind-memory-central-spectrum',
           'susskind-memory-oriented-kernel',
           'susskind-memory-two-component']
MUTABLE_PLANS = {'README.md', 'PLAN.md', 'NEXT.md'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text(encoding='utf-8'))


def read_build_text(path):
    raw = path.read_bytes()
    encoding = 'utf-16' if raw.startswith((b'\xff\xfe', b'\xfe\xff')) else 'utf-8-sig'
    return raw.decode(encoding, errors='replace')


def verify_hash_map(mapping, label):
    for name, digest in mapping.items():
        path = ROOT / name
        assert path.is_file(), (label, 'missing input', name)
        assert sha(path) == digest, (label, 'changed input', name)


def verify_parents():
    parents = {}
    for stage in PARENTS:
        path = ROOT / 'results' / stage / 'local-evidence.json'
        record = load(path)
        checked = []
        for name, digest in record['hashes'].items():
            if name in MUTABLE_PLANS:
                continue
            assert (ROOT/name).is_file(), (stage, 'missing predecessor', name)
            assert sha(ROOT/name) == digest, (stage, 'changed predecessor', name)
            checked.append(name)
        parents[path.relative_to(ROOT).as_posix()] = dict(
            sha256=sha(path), unchanged_scientific_files=checked,
            excluded_mutable_plans=sorted(MUTABLE_PLANS.intersection(record['hashes'])))
    return parents


def check_finite(record):
    fixtures = record['fixture_checks']['rows']
    chaos = record['highest_chaos_checks']['rows']
    clifford = record['Clifford_size']['rows']
    counts = dict(fixture_sector_word_cases=len(fixtures),
                  fixture_identities=sum(row['identities'] for row in fixtures),
                  highest_chaos_sector_coefficient_cases=sum(row['homogeneous_coefficients'] for row in chaos),
                  highest_chaos_identities=sum(row['identities'] for row in chaos),
                  even_clifford_strings=sum(row['total_even_strings'] for row in clifford))
    expected = dict(fixture_sector_word_cases=96, fixture_identities=384,
                    highest_chaos_sector_coefficient_cases=1428,
                    highest_chaos_identities=5712, even_clifford_strings=168)
    assert counts == record['counts'] == expected
    assert all(row['identities'] == 4 for row in fixtures)
    assert sum(row['resonance_factor'] == 0 for row in fixtures) == 6
    assert sum(row['nonzero_finite_bridge_error'] for row in fixtures) == 96
    assert {(row['j'], row['m']) for row in chaos} == {(j, m) for j in range(1, 5) for m in [1, 2]}
    for row in chaos:
        assert row['coordinates'] == 9 and row['N'] == 3 and row['p'] == 1
        assert row['homogeneous_coefficients'] == comb(9+row['j']-1, row['j'])
        assert row['identities'] == 4*row['homogeneous_coefficients']
    for row in clifford:
        assert sum(row['even_strings_by_length'].values()) == row['total_even_strings']
        assert row['total_even_strings'] == 2**(2*row['N']-1)
    assert record['arithmetic']['largest_product_bound'] < 2**48
    return dict(counts=counts, finite_resonances=6, nonzero_fixture_bridge_errors=96,
                maximum_integer_product_bound=record['arithmetic']['largest_product_bound'])


def check_marked(record):
    rows, moments = record['finite_records'], record['moment_rows']
    assert len(rows) == record['finite_words'] == 160
    comparisons = sum(len(row['all_mark_subsets']) for row in rows)
    assert comparisons == record['finite_comparisons'] == 640
    assert all(len(row['all_mark_subsets']) == 4 for row in rows)
    for row in rows:
        for value in row['all_mark_subsets']:
            F(value)
    assert len(moments) == 54
    keys = [(row['j'], row['k'], row['left_energy'], row['right_energy'], row['mark_mask'])
            for row in moments]
    assert len(set(keys)) == len(keys)
    for row in moments:
        assert row['mark_mask'] in [0, 1, 3] and row['allowed_pairings'] >= 0
        for term in row['polynomial']:
            assert term['lambda_power'] >= 0 and term['q_power'] >= 0
            F(term['coefficient'])
    return dict(finite_words=len(rows), finite_comparisons=comparisons,
                limiting_moment_polynomials=len(moments))


def check_geometry(record):
    rows = record['geometry_layouts']
    totals = Counter()
    for row in rows:
        totals.update(row['counts'])
    expected = dict(all_pairings=8038, allowed_pairings=5716,
                    exact_coefficient_comparisons=25884, single_positive=9450,
                    double_positive=2598, double_zero=4429, single_zero=2274,
                    connected_degree_at_least_three_zero=1085,
                    single_negative=298, double_negative=34)
    assert dict(totals) == record['geometry_totals'] == expected
    assert len(rows) == 34
    assert sum(row['counts'].get('allowed_pairings', 0) > 0 for row in rows) == 33
    encoded = json.dumps(record['parameters'], sort_keys=True, separators=(',', ':')).encode()
    assert hashlib.sha256(encoded).hexdigest() == record['parameter_sha256']
    assert record['predecessor_output_sha256'] == sha(OUT/'marked-chord-checks.json')
    bridge, u1 = record['all_bridge_moments'], record['U1_two_energy_moments']
    assert [row['j'] for row in bridge] == list(range(1, 6))
    assert all(row['all_bridge_pairings'] == factorial(row['j']) for row in bridge)
    assert sum(row['all_bridge_pairings'] for row in bridge) == 153
    assert all(row['closed_log_derivative_form_verified'] for row in bridge)
    assert len(u1) == 16
    assert {(row['left_energy'], row['right_energy']) for row in u1} == {
        (left, total-left) for total in [0, 2, 4, 6] for left in range(total+1)}
    assert all(row['independent_vacuum_and_Mehler_recurrences_verified'] for row in u1)
    return dict(layouts=34, nonempty_layouts=33, geometry_totals=expected,
                all_bridge_degree_cases=5, all_bridge_pairings=153,
                U1_two_energy_cases=16, parameter_sha256=record['parameter_sha256'],
                predecessor_output_sha256=record['predecessor_output_sha256'])


def check_law(record):
    algebra = record['neutral_algebra']['rows']
    integers = record['integer_weights']['rows']
    probabilities = record['probabilities']['rows']
    counts = dict(neutral_matrix_units=sum(row['complete_neutral_matrix_units'] for row in algebra),
                  superoperator_systems=len(algebra),
                  exact_integer_weight_powers=sum(row['tested_power_range'][1]-
                                                 row['tested_power_range'][0]+1 for row in integers),
                  probability_cases=len(probabilities),
                  MGF_cases=sum(len(row['MGF']) for row in probabilities))
    expected = dict(neutral_matrix_units=96, superoperator_systems=3,
                    exact_integer_weight_powers=1542, probability_cases=18, MGF_cases=54)
    assert counts == record['counts'] == expected
    for row in algebra:
        assert row['complete_neutral_matrix_units'] == comb(2*row['N'], row['N'])
        assert row['exact_superoperator_identities'] == 5 and row['rank_less_equal_half_size']
    assert [row['j'] for row in integers] == list(range(1, 7))
    assert all(row['all_power_sums_nonnegative'] and row['tested_power_range'] == [0, 256]
               and row['permutation_count'] == factorial(row['j']) for row in integers)
    assert {(row['j'], row['lambda_value']) for row in probabilities} == {
        (j, lam) for j in range(1, 7) for lam in ['1/10', '1/4', '1']}
    cutoff_values, tail_bounds = [], []
    with mp.workdps(120):
        assert mp.mpf(record['probabilities']['max_scaled_precision_change']) < mp.mpf('4e-78')
        missing = []
        for row in probabilities:
            normalization = row['normalization']
            assert len(normalization['pmf']) == normalization['cutoff']+1
            assert all(mp.mpf(value) >= 0 for value in normalization['pmf'])
            assert {item['z'] for item in row['MGF']} == {'1/2', '1', '3/2'}
            for item in row['MGF']:
                bound = F(item['exact_rational_tail_upper'])
                assert 0 < bound <= F(1, 10**42)
                cutoff_values.append(item['cutoff'])
                tail_bounds.append(bound)
                missing.append(mp.mpf(item['numerical_missing_tail']))
                if item['z'] == '1':
                    assert item['cutoff'] == normalization['cutoff']
                    assert mp.mpf(item['normalized_MGF']) == 1
            assert F(normalization['exact_rational_first_moment_tail_upper']) > 0
            assert F(normalization['exact_rational_second_moment_tail_upper']) > 0
        largest_missing = mp.nstr(max(missing), 45)
    assert [min(cutoff_values), max(cutoff_values)] == [19, 227]
    return dict(counts=counts, cutoff_range=[19, 227],
                maximum_exact_tail_upper=str(max(tail_bounds)),
                maximum_numerical_missing_tail=largest_missing,
                precision_change=record['probabilities']['max_scaled_precision_change'],
                largest_integer_product_bound=record['largest_integer_product_bound'])


def check_plot(law):
    metadata = load(OUT/'deficit-law-plot.json')
    # The plot schema deliberately has no status field and has one output hash.
    verify_hash_map(metadata['input_sha256'], 'deficit-law-plot.json')
    assert metadata['output_sha256'] == sha(OUT/'deficit-law.png')
    rows = metadata['rows']
    assert len(rows) == 4 and [row['j'] for row in rows] == [1, 2, 4, 6]
    for row in rows:
        saved = next(item for item in law['probabilities']['rows']
                     if item['j'] == row['j'] and item['lambda_value'] == row['lambda_value'])
        assert row['lambda_value'] == '1/4' and row['display_cutoff'] == 28
        assert row['saved_summation_cutoff'] == saved['normalization']['cutoff']
        actual = sum(float(value) for value in saved['normalization']['pmf'][29:])
        assert actual == row['saved_probability_beyond_display']
    return metadata


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
        assert record['status'] == 'passed', name
        verify_hash_map(record['input_sha256'], name)
    checks = dict(finite_bridge=check_finite(data['finite-bridge.json']),
                  marked_chords=check_marked(data['marked-chord-checks.json']),
                  geometry=check_geometry(data['geometry-checks.json']),
                  rank_law=check_law(data['rank-law-checks.json']))
    plot = check_plot(data['rank-law-checks.json'])

    source = ROOT/'manuscript'/STAGE/'main.tex'
    pdf = ROOT/'output/pdf'/f'de-sitter-{STAGE}.pdf'
    build = ROOT/'tmp/pdfs'/f'{STAGE}-build'
    final_log = read_build_text(OUT/'main.log')
    assert page_count(final_log) == EXPECTED_PAGES
    assert not re.search(r'undefined|Rerun to get|Overfull|Fatal error', final_log)
    assert sha(pdf) == sha(build/'main.pdf')
    assert sha(OUT/'main.log') == sha(build/'main.log')
    assert (build/'main.pdf').stat().st_mtime_ns >= source.stat().st_mtime_ns, 'Source newer than final PDF'
    passes = [build/f'pass{number}.txt' for number in [1, 2]]
    pass_records = []
    for number, path in enumerate(passes, 1):
        text = read_build_text(path)
        assert page_count(text) == EXPECTED_PAGES
        assert 'Fatal error' not in text and '! Emergency stop' not in text
        if number == 2:
            assert not re.search(r'undefined|Rerun to get|Overfull', text)
        pass_records.append(dict(pass_number=number, path=path.relative_to(ROOT).as_posix(),
                                 sha256=sha(path), compiler_header=text.splitlines()[0]))
    assert passes[1].stat().st_mtime_ns >= passes[0].stat().st_mtime_ns >= source.stat().st_mtime_ns
    # Bind only final zero-padded pages; old page-1.png..page-9.png are stale.
    pages = [build/f'page-{number:02d}.png' for number in range(1, EXPECTED_PAGES+1)]
    assert all(path.is_file() for path in pages)
    assert all(path.stat().st_mtime_ns >= pdf.stat().st_mtime_ns for path in pages), 'Stale final render'

    names = [
        'code/check_neutral_size_bridge.py', 'code/memory_size_marked_chords.py',
        'code/check_size_chord_geometry.py', 'code/check_size_rank_law.py',
        'code/plot_size_deficit_law.py', 'code/build_memory_size_bridge.ps1',
        'code/record_memory_size_bridge.py',
        f'manuscript/{STAGE}/main.tex', f'output/pdf/de-sitter-{STAGE}.pdf',
        f'notes/{STAGE}.md', f'literature/{STAGE}-sources.md',
        f'journal/2026-11-07-{STAGE}.md', 'README.md', 'PROJECT.md', 'PLAN.md', 'NEXT.md']
    diagnostics = sorted(path for path in (OUT/'diagnostics').rglob('*') if path.is_file())
    assert diagnostics, 'Expected preserved initial-build diagnostics'
    files = [ROOT/name for name in names]+[OUT/name for name in SCIENCE]
    files += [OUT/'main.log', OUT/'deficit-law-plot.json', OUT/'deficit-law.png']+diagnostics
    assert all(path.is_file() for path in files), [str(path) for path in files if not path.is_file()]
    assert len(files) == len(set(files))
    result = dict(
        status='Local research progress; long-term objective remains active',
        created_at_utc=datetime.now(timezone.utc).isoformat(),
        scientific_status='Fixed-degree rank marks have a strong central/full-trace bridge '
                          'to half-Clifford-size defects. Common-coupling marked chord geometry '
                          'gives their fixed-energy moments, fixed-j deficit law and determinate '
                          'positive two-energy measures. The j=1 static regression and Mehler '
                          'remainder are explicit. The forced neutral descendants, their static '
                          'part and the full charged-memory cross kernel remain unresolved.',
        parameters='Standard CAR, overlapping complex Gaussian p-body, common original J '
                   'across charges, normalized annealed beta=mu=0 full trace, hbar=1 and string '
                   'real time. Even N, lambda_N=p^2/N tends to positive lambda, fixed Js; '
                   'all chaos and energy degrees are fixed before this limit.',
        proof_scope='Finite Casimir composition is performed before charge restriction; '
                    'overlap marks and rank tails control the strong bridge. The auxiliary '
                    'full-neutral spin algebra distinguishes rank from Clifford length. '
                    'The scientific JSONs verify finite algebra, local diagrams and analytic '
                    'candidate-law formulas; the uniform asymptotic statements use the '
                    'separate manuscript proofs. No finite list of numerical checks establishes '
                    'an actual rank-MGF or finite-time kernel limit.',
        continuation='Apply the local mark rule to the forced Gaussian-Wick descendants '
                     'while retaining every forbidden internal contraction; identify zero-frequency '
                     'mass before passing exact samplewise pinching to the limit, then add '
                     'the charged defect with its reverse source in the same original-J calculation.',
        external_scope='Earlier Issue 1 review retained. This stage makes no new GitHub '
                       'reply or publication and no full reintegration of the historical audit.',
        git_commit_created=False, git_push_performed=False, archive_created=False,
        external_reply_sent=False,
        previous_evidence=parents,
        input_hashes_verified={name: record['input_sha256'] for name, record in data.items()},
        geometry_predecessor_verified=checks['geometry']['predecessor_output_sha256'],
        plot_inputs_verified=plot['input_sha256'], plot_output_verified=plot['output_sha256'],
        checks=checks,
        evidence_scopes={name: record['scope'] for name, record in data.items()},
        recorded_driver_environments={name: record['environment'] for name, record in data.items()},
        recorded_plot_matplotlib_version=plot['matplotlib_version'],
        environment=dict(python=platform.python_version(), numpy=np.__version__,
                         sympy=sp.__version__, mpmath=mp.__version__, matplotlib=matplotlib.__version__),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_size_bridge.ps1',
                   pages=EXPECTED_PAGES, passes=2, pass_records=pass_records,
                   compiler_header=final_log.splitlines()[0],
                   render_command=f'pdftoppm -r 95 -png output/pdf/de-sitter-{STAGE}.pdf '
                                  f'tmp/pdfs/{STAGE}-build/page',
                   render_naming='Only final page-01.png through page-10.png are bound; '
                                 'un-padded page-1.png through page-9.png are excluded.',
                   final_pages_visually_inspected=list(range(1, EXPECTED_PAGES+1)),
                   visual_result='Coordinating agent inspected all final zero-padded rendered '
                                 'pages after the last build, including the corrected summation '
                                 'indices, equations, figure, references and page boundaries.'),
        diagnostics=dict(preserved_files=[path.relative_to(ROOT).as_posix() for path in diagnostics],
                         scope='Earlier build source/log/pass records are retained separately '
                               'from the validated final compilation and zero-padded renders.'),
        process_state='Recorder is run only after scientific drivers, both final build passes '
                      'and rendering have terminated; it executes none of them.',
        hashes={path.relative_to(ROOT).as_posix(): sha(path) for path in files},
        render_hashes={path.relative_to(ROOT).as_posix(): sha(path) for path in pages},
        build_pass_hashes={path.relative_to(ROOT).as_posix(): sha(path) for path in passes})
    with destination.open('x', encoding='utf-8') as handle:
        handle.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(dict(status='passed', pages=EXPECTED_PAGES, current_files=len(files),
                         scientific_outputs=len(SCIENCE),
                         previous_scientific_files_verified=sum(len(parent['unchanged_scientific_files'])
                                                                for parent in parents.values()),
                         final_zero_padded_renders=len(pages), final_build_passes=len(passes))))


if __name__ == '__main__':
    main()
