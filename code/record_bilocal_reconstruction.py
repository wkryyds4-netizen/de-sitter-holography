"""Seal saved bilocal evidence after both builds and actual final-page review.

This recorder reads artifacts only: no scientific main, compiler, renderer or
Git command runs. Use --visual-review-complete only after every final page has
actually been inspected. Successful evidence is exclusive-create.
"""
import argparse
from datetime import datetime, timezone
from fractions import Fraction
from itertools import product
from math import isfinite, isclose
from pathlib import Path
import hashlib
import json
import platform
import re
import struct

ROOT = Path(__file__).resolve().parent.parent
STAGE = 'susskind-bilocal-reconstruction'
OUT = ROOT/'results'/STAGE
SCIENCE = ['finite-volterra-checks.json', 'spectral-filter-checks.json']
PARENTS = ['susskind-memory-joint-regression', 'susskind-memory-return-expansion',
           'susskind-memory-uniform-bound', 'susskind-bilocal-channel']
MUTABLE_PLANS = {'README.md', 'PLAN.md', 'NEXT.md'}
REQUIRED_DOCS = [f'notes/{STAGE}.md', f'notes/{STAGE}-proof.md',
                 'notes/susskind-bilocal-multiplier-proof.md',
                 f'literature/{STAGE}-sources.md', f'journal/2026-11-07-{STAGE}.md',
                 'README.md', 'PROJECT.md', 'PLAN.md', 'NEXT.md']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text(encoding='utf-8'))


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def read_text(path):
    raw = path.read_bytes()
    encoding = 'utf-16' if raw.startswith((b'\xff\xfe', b'\xfe\xff')) else 'utf-8-sig'
    return raw.decode(encoding, errors='replace')


def finite_numbers(value):
    if isinstance(value, dict):
        for item in value.values():
            finite_numbers(item)
    elif isinstance(value, list):
        for item in value:
            finite_numbers(item)
    elif isinstance(value, float):
        assert isfinite(value), 'Nonfinite saved numeric value'


def verify_inputs(record, name):
    assert record['status'] == 'passed' and record['scope'] and record['environment'], name
    assert record['input_sha256'] and record['parameters']
    finite_numbers(record)
    for filename, expected in record['input_sha256'].items():
        assert (ROOT/filename).is_file(), (name, 'missing input', filename)
        assert sha(ROOT/filename) == expected, (name, 'changed input', filename)
    parameter_hash = digest(record['parameters'])
    if name == 'spectral-filter-checks.json':
        assert record['schema_version'] == 1
        assert record['parameter_sha256'] == parameter_hash
    else:
        assert name == 'finite-volterra-checks.json'
        assert 'parameter_sha256' not in record and 'schema_version' not in record
    return {'sha256': parameter_hash,
            'origin': 'Verified against saved digest' if 'parameter_sha256' in record else
                      'Computed now from the original finite driver parameter record; no digest was stored there'}


def verify_parents():
    result = {}
    for stage in PARENTS:
        path = ROOT/'results'/stage/'local-evidence.json'
        record = load(path)
        checked = []
        for filename, expected in record['hashes'].items():
            if filename in MUTABLE_PLANS:
                continue
            assert (ROOT/filename).is_file(), (stage, 'missing predecessor', filename)
            assert sha(ROOT/filename) == expected, (stage, 'changed predecessor', filename)
            checked.append(filename)
        result[path.relative_to(ROOT).as_posix()] = {
            'sha256': sha(path), 'unchanged_scientific_files': checked,
            'excluded_mutable_plans': sorted(MUTABLE_PLANS.intersection(record['hashes']))}
    return result


def verify_finite(record):
    p, rows = record['parameters'], record['records']
    assert [(r['N'], r['p'], r['m'], r['seed']) for r in p['CAR']] == [
        (4, 2, 3, 110701), (5, 2, 3, 110702), (4, 2, 2, 110703)]
    assert p['quadrature_orders'] == [8, 14] and p['time_tests'] == [.37, .91]
    assert p['laplace_z'] == [[.7, .9], [1.3, -.4]]
    assert p['small_time_steps'] == [.04, .02, .01]
    assert p['generic'] == dict(seed=110704, h_eigenvalues=[-.4, -.4, .7],
        bath_dimension=4, alpha=.64, m=3, coupling_scale=.25)
    assert p['matrix_tolerance'] == 2e-11 and p['quadrature_tolerance'] == 3e-10
    counts = dict(systems=len(rows), actual_CAR_systems=sum(r['parameters']['kind']=='actual_complex_CAR' for r in rows),
        degenerate_bath_systems=sum(r['parameters']['kind']=='independent_degenerate_Hermitian_bath' for r in rows),
        Laplace_points=sum(len(r['resolvent_checks']) for r in rows),
        time_quadrature_runs=sum(len(t['quadrature']) for r in rows for t in r['time_checks']),
        contraction_times=sum(len(r['contraction_checks']) for r in rows),
        small_time_curvatures=sum(len(r['small_time_curvature']) for r in rows))
    assert counts == record['counts'] == dict(systems=4, actual_CAR_systems=3, degenerate_bath_systems=1,
        Laplace_points=8, time_quadrature_runs=16, contraction_times=24, small_time_curvatures=12)
    maximum_matrix = maximum_quadrature = 0.
    finite_remainders = []
    for index, r in enumerate(rows):
        spec = r['parameters']; d = r['system_dimension']; k = r['bath_dimension']
        if index < 3:
            assert all(spec[key] == value for key, value in p['CAR'][index].items())
            assert spec['force_initial_identity_error'] < p['matrix_tolerance']
        else:
            assert all(spec[key] == value for key, value in p['generic'].items())
        assert sum(r['energy_block_dimensions']) == d
        assert {name: value['shape'] for name, value in r['input_matrices'].items()} == {
            'h': [d, d], 'a': [d, d], 'bath': [k, k], 'R': [k, d]}
        for matrix in r['input_matrices'].values():
            a, b = matrix['shape']
            assert all(len(matrix[part]) == a and all(len(row)==b for row in matrix[part])
                       for part in ['real', 'imag'])
        assert [x['time'] for x in r['contraction_checks']] == p['contraction_times']
        for item in r['contraction_checks']:
            assert 0 <= item['F0_operator_norm'] <= 1+p['matrix_tolerance']
            assert item['block_dilation_dimensions'] == [k+n for n in r['energy_block_dimensions']]
        assert [x['z'] for x in r['resolvent_checks']] == p['laplace_z']
        for item in r['resolvent_checks']:
            errors = [value for name, value in item.items() if name.endswith('_error')]
            assert len(errors) == 4 and min(errors) >= 0 and max(errors) < p['matrix_tolerance']
            maximum_matrix = max(maximum_matrix, *errors)
        assert [x['time'] for x in r['time_checks']] == p['time_tests']
        for item in r['time_checks']:
            errors = [value for name, value in item.items() if name.endswith('_error')]
            assert min(errors) >= 0 and max(errors) < p['matrix_tolerance']
            assert item['background_after_pinching_HS'] < p['matrix_tolerance']
            maximum_matrix = max(maximum_matrix, *errors)
            assert [x['order'] for x in item['quadrature']] == p['quadrature_orders']
            for scan in item['quadrature']:
                errors = [value for name, value in scan.items() if name.endswith('_error')]
                assert len(errors) == 5 and min(errors) >= 0 and max(errors) < p['quadrature_tolerance']
                maximum_quadrature = max(maximum_quadrature, *errors)
            changes = item['resolution_differences']
            assert set(changes) == {'_linear', '_remainder', '_filter'}
            assert min(changes.values()) >= 0 and max(changes.values()) < p['quadrature_tolerance']
        assert [x['time'] for x in r['small_time_curvature']] == p['small_time_steps']
        if r['dynamic_initial_memory_HS'] > p['finite_nonzero_threshold']:
            errors = [x['even_curvature_relative_error'] for x in r['small_time_curvature']]
            assert errors[1] < errors[0]/3 and errors[2] < errors[1]/3
            assert max(x['incorrect_order_difference_HS'] for x in r['resolvent_checks']) > p['finite_nonzero_threshold']
            remainder = r['time_checks'][-1]['quadrature'][-1]['finite_projected_scaled_remainder_HS']
            assert remainder > p['finite_nonzero_threshold']
            finite_remainders.append(remainder)
        else:
            assert spec['zero_memory_control'] is True and spec['m'] == spec['p']
        if index == 3:
            assert r['energy_block_dimensions'] == [2, 1]
            assert r['degenerate_background_commutator_HS'] > 1e-5
            assert r['incorrect_diagonal_only_pinching_loss_HS'] > 1e-3
    return dict(counts=counts, maximum_matrix_identity_error=maximum_matrix,
        maximum_time_identity_error=maximum_quadrature,
        finite_projected_scaled_remainders_at_t091=finite_remainders,
        scope='Finite deterministic checks retain nonzero remainders and a genuinely noncommuting '
              'degenerate background. They do not measure a DS convergence rate.')


def z(pair):
    assert len(pair) == 2
    return complex(*pair)


def verify_spectral(record):
    p, exact, numerical = record['parameters'], record['exact'], record['numerical']
    assert p['exact_degree'] == 14 and p['field_degrees'] == list(range(1, 9))
    assert p['rational_q'] == ['1/3', '3/4'] and p['rational_C'] == ['1/2', '2']
    assert p['rational_energies'] == [['-2', '3'], ['1/3', '1/3'], ['0', '0']]
    assert p['quadrature_sizes'] == [[40, 24], [64, 36]] and p['floating_tolerance'] == 3e-10
    exact_keys = {(r['q'], r['C'], tuple(r['energies']), r['degree']) for r in exact['rows']}
    assert exact_keys == set(product(p['rational_q'], p['rational_C'], map(tuple, p['rational_energies']), p['field_degrees']))
    assert len(exact_keys) == len(exact['rows']) == 96
    charged = neutral = 0
    for r in exact['rows']:
        assert r['exact'] is True and r['charged_checked'] == (r['degree'] >= 2)
        n = list(map(Fraction, r['neutral']))
        assert len(n) == 15 and n[:3] == [0, 0, 0]
        assert n[3] == 2*Fraction(r['C'])*int(r['degree']==1)
        neutral += len(n)
        if r['charged_checked']:
            c = list(map(Fraction, r['charged']))
            assert len(c) == 15 and c[:2] == [0, 0] and c[2] == int(r['degree']==2)
            charged += len(c)
        else:
            assert r['charged'] is None
    keys = {(r['lambda_value'], r['time'], tuple(r['relative_energies']), r['degree']) for r in numerical['rows']}
    assert keys == set(product(p['physical_lambda'], p['times'], map(tuple, p['relative_energy_pairs']), p['numerical_field_degrees']))
    assert len(keys) == len(numerical['rows']) == numerical['cases'] == 144
    for r in numerical['rows']:
        assert r['charged_checked'] == (r['degree'] >= 2)
        fine, coarse = r['fine'], r['coarse']
        identity = [abs(z(fine['neutral'])-z(fine['neutral_target']))]
        if r['charged_checked']:
            identity.append(abs(z(fine['charged'])-z(fine['charged_target'])))
        refinement = max(abs(z(fine[k])-z(coarse[k])) for k in fine)
        assert isclose(max(identity), r['identity_error'], rel_tol=1e-12, abs_tol=1e-25)
        assert isclose(refinement, r['refinement_error'], rel_tol=1e-12, abs_tol=1e-25)
        assert max(identity+[refinement]) < p['floating_tolerance']
    special = numerical['repeated_nodes']
    assert {(tuple(r['nodes']), r['time']) for r in special} == set(product(map(tuple, p['scalar_repeated_nodes']), p['times']))
    assert len(special) == 20
    for r in special:
        assert abs(z(r['expm'])-z(r['simplex'])) < p['floating_tolerance']
        assert isclose(abs(z(r['expm'])-z(r['simplex'])), r['error'], rel_tol=1e-12, abs_tol=1e-25)
        assert abs(z(r['expm'])) <= r['time']**2/2+1e-13
    for name, values in [('maximum_identity_error', [r['identity_error'] for r in numerical['rows']]),
                         ('maximum_refinement_error', [r['refinement_error'] for r in numerical['rows']]),
                         ('maximum_repeated_node_error', [r['error'] for r in special])]:
        assert numerical[name] == max(values) < p['floating_tolerance']
    assert charged == exact['charged_coefficient_identities'] == 1260
    assert neutral == exact['neutral_coefficient_identities'] == 1440
    counts = dict(charged_coefficient_identities=charged, neutral_coefficient_identities=neutral,
                  floating_field_cases=len(numerical['rows']), repeated_scalar_cases=len(special))
    assert counts == record['counts']
    return dict(counts=counts, maximum_identity_error=numerical['maximum_identity_error'],
        maximum_refinement_error=numerical['maximum_refinement_error'],
        maximum_repeated_node_error=numerical['maximum_repeated_node_error'],
        scope='Saved rational coefficient checks and two-resolution quadratures are inspected. '
              'The charged d=1 case is deliberately excluded from the cancellation identity; '
              'coincident and nearly coincident energies remain in the record.')


def page_count(log):
    match = re.search(r'Output written on .*?\((\d+) pages?,', log, re.S)
    assert match, 'Missing successful PDF output line'
    return int(match.group(1))


def png_metadata(path):
    raw = path.read_bytes()
    assert raw[:8] == b'\x89PNG\r\n\x1a\n'
    offset = 8; size = density = None
    while offset < len(raw):
        length = struct.unpack('>I', raw[offset:offset+4])[0]
        kind = raw[offset+4:offset+8]; data = raw[offset+8:offset+8+length]
        if kind == b'IHDR':
            size = struct.unpack('>II', data[:8])
        elif kind == b'pHYs':
            density = struct.unpack('>IIB', data)
        offset += length+12
        if kind == b'IEND':
            break
    assert size and density == (3740, 3740, 1), (path, 'not a 95 dpi render', density)
    return dict(width_pixels=size[0], height_pixels=size[1], pixels_per_metre=list(density[:2]),
                density_units='metres', nominal_dpi=95)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--visual-review-complete', action='store_true',
                        help='Attest that all final rendered pages were actually inspected')
    args = parser.parse_args()
    assert args.visual_review_complete, 'Run only after final-page review; supply --visual-review-complete'
    destination = OUT/'local-evidence.json'
    if destination.exists():
        raise FileExistsError(destination)
    parents = verify_parents()
    science = {name: load(OUT/name) for name in SCIENCE}
    parameter_hashes = {name: verify_inputs(record, name) for name, record in science.items()}
    checks = {SCIENCE[0]: verify_finite(science[SCIENCE[0]]), SCIENCE[1]: verify_spectral(science[SCIENCE[1]])}
    source = ROOT/'manuscript'/STAGE/'main.tex'
    pdf = ROOT/'output/pdf'/f'de-sitter-{STAGE}.pdf'
    build = ROOT/'tmp/pdfs'/f'{STAGE}-build'
    render = ROOT/'tmp/pdfs'/f'{STAGE}-render'
    log = read_text(OUT/'main.log'); pages_count = page_count(log)
    assert pages_count > 0 and not re.search(r'undefined|Rerun to get|Overfull|Fatal error', log)
    assert sha(pdf) == sha(build/'main.pdf') and sha(OUT/'main.log') == sha(build/'main.log')
    assert (build/'main.pdf').stat().st_mtime_ns >= source.stat().st_mtime_ns
    passes = [build/f'pass{i}.txt' for i in [1, 2]]
    pass_records = []
    for i, path in enumerate(passes, 1):
        content = read_text(path)
        assert page_count(content) == pages_count
        assert 'Fatal error' not in content and '! Emergency stop' not in content
        if i == 2:
            assert not re.search(r'undefined|Rerun to get|Overfull', content)
        pass_records.append(dict(pass_number=i, path=path.relative_to(ROOT).as_posix(),
                                 sha256=sha(path), compiler_header=content.splitlines()[0]))
    assert passes[1].stat().st_mtime_ns >= passes[0].stat().st_mtime_ns >= source.stat().st_mtime_ns
    width = len(str(pages_count))
    pages = [render/f'page-{i:0{width}d}.png' for i in range(1, pages_count+1)]
    assert set(render.glob('page-*.png')) == set(pages), 'Final render filenames/count mismatch'
    assert all(path.stat().st_mtime_ns >= pdf.stat().st_mtime_ns for path in pages), 'Stale render'
    render_metadata = {path.relative_to(ROOT).as_posix(): png_metadata(path) for path in pages}
    names = set(REQUIRED_DOCS) | {'code/build_bilocal_reconstruction.ps1',
        'code/record_bilocal_reconstruction.py', f'manuscript/{STAGE}/main.tex',
        f'output/pdf/de-sitter-{STAGE}.pdf', f'results/{STAGE}/main.log'}
    for record in science.values():
        names.update(record['input_sha256'])
    names.update(f'results/{STAGE}/{name}' for name in SCIENCE)
    diagnostics = sorted(path for path in (OUT/'diagnostics').rglob('*') if path.is_file())
    files = [ROOT/name for name in sorted(names)]+diagnostics
    assert all(path.is_file() for path in files), [str(path) for path in files if not path.is_file()]
    assert len(files) == len(set(files))
    result = dict(status='Local research milestone; long-term de Sitter research remains open',
        created_at_utc=datetime.now(timezone.utc).isoformat(),
        scientific_status='The stage reconstructs the background-subtracted bilocal covariance from '
            'the completed common-J memory kernel. Contractive background dilation and L2/L8 '
            'interpolation control the strong return error by O(p^(-1/3)) on compact real-time '
            'windows. Charged coefficients simplify and neutral coefficients have a regular '
            'divided-difference representation.',
        proof_scope='Finite matrix checks and limiting coefficient quadratures verify the recorded '
            'interfaces; they are not DS extrapolations. Uniform full-charge error and the '
            'multiplier/series passage are justified in the proof notes. No canonical boson '
            'algebra, microscopic P-minus, bulk mass identification or finite-G interaction '
            'is inferred from covariance alone.',
        parameters='Standard CAR, common complex Gaussian p-body couplings, normalized annealed '
            'beta=mu=0 full-Fock trace, hbar=1, string time; even N, p^2/N in a positive '
            'compact set approaching lambda>0, fixed Js. No growing-time/small-lambda exchange.',
        external_scope='No new issue reply, publication, historical reintegration, Git operation or archive.',
        git_commit_created=False, git_push_performed=False, archive_created=False, external_reply_sent=False,
        previous_evidence=parents, checks=checks,
        input_hashes_verified={name: record['input_sha256'] for name, record in science.items()},
        parameter_hashes_verified=parameter_hashes,
        evidence_scopes={name: record['scope'] for name, record in science.items()},
        recorded_driver_environments={name: record['environment'] for name, record in science.items()},
        environment=dict(python=platform.python_version(), platform=platform.platform(),
                         recorder_dependencies='Python standard library only'),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_bilocal_reconstruction.ps1',
            pages=pages_count, passes=2, pass_records=pass_records, compiler_header=log.splitlines()[0],
            render_command=f'pdftoppm -r 95 -png output/pdf/de-sitter-{STAGE}.pdf tmp/pdfs/{STAGE}-render/page',
            final_render_filenames=[path.name for path in pages], render_metadata=render_metadata,
            final_pages_visually_inspected=list(range(1, pages_count+1)),
            visual_review_attestation='Caller supplied --visual-review-complete after inspecting all '
                'final rendered pages. The recorder verifies artifact hashes, freshness, counts and '
                'PNG resolution; it does not itself perform visual inspection.'),
        diagnostics=dict(preserved_files=[path.relative_to(ROOT).as_posix() for path in diagnostics]),
        process_state='Both science outputs, final two-pass build and rendering are saved terminal '
                      'artifacts. This recorder executes none of those processes.',
        hashes={path.relative_to(ROOT).as_posix(): sha(path) for path in files},
        render_hashes={path.relative_to(ROOT).as_posix(): sha(path) for path in pages},
        build_pass_hashes={path.relative_to(ROOT).as_posix(): sha(path) for path in passes})
    with destination.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(status='passed', pages=pages_count, scientific_outputs=len(SCIENCE),
        current_files=len(files), final_renders=len(pages), final_build_passes=2,
        previous_scientific_files_verified=sum(len(r['unchanged_scientific_files']) for r in parents.values()))))


if __name__ == '__main__':
    main()
