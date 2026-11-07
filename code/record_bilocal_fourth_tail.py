"""Read-only evidence seal; run only after final build/render/visual review.

The three saved science schemas are validated separately. Only saved scalar
arithmetic and combinatorial bookkeeping are replayed: no CAR matrices,
scientific drivers, compiler, renderer, or Git command are executed.
"""
import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, factorial, isclose, prod, sqrt
from pathlib import Path
import json
import platform
import re

from record_bilocal_reconstruction import (
    sha, load, digest, read_text, finite_numbers, page_count, png_metadata,
)

ROOT = Path(__file__).resolve().parent.parent
STAGE = 'susskind-bilocal-fourth-tail'
OUT = ROOT/'results'/STAGE
SCIENCE = ['tagged-covariance-checks.json', 'coherence-channel-checks.json',
           'charged-reassembly-checks.json']
SCIENCE_SHA256 = {
    'tagged-covariance-checks.json': '9169292b18de235437ebf22bd69b40b532da15e8d45e2e82040734a18d723985',
    'coherence-channel-checks.json': '99b2e3d852c5e706853c703310ee8fa2494aa3ac611bef59e5fb4dd89f3e430a',
    'charged-reassembly-checks.json': 'bb76a4b5b195ecb203612afbf9b285e9d91fd59bb5c57ace47018d65aa18a0b4',
}
PARENTS = ['susskind-bilocal-four-point', 'susskind-memory-neutral-marked-series',
           'susskind-memory-oriented-transfer']
MUTABLE = {'README.md', 'PLAN.md', 'NEXT.md'}
DOCS = [f'notes/{STAGE}.md', f'notes/{STAGE}-diagram-draft.md',
        f'notes/{STAGE}-operator-draft.md', f'notes/{STAGE}-charged-proof.md',
        'notes/susskind-bilocal-charged-tail-audit.md',
        'notes/susskind-bilocal-full-four-point-transfer.md',
        f'literature/{STAGE}-sources.md', 'literature/susskind-bilocal-clifford-riesz-sources.md',
        f'journal/2026-11-07-{STAGE}.md', 'README.md', 'PROJECT.md', 'PLAN.md', 'NEXT.md']


def near(a, b):
    assert isclose(a, b, rel_tol=2e-12, abs_tol=1e-25), (a, b)


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def inputs(record, name):
    assert sha(OUT/name) == SCIENCE_SHA256[name], (name, 'sealed output changed')
    finite_numbers(record)
    assert record['scope']
    if name == SCIENCE[0]:
        assert record['success'] is True and record['failures'] == []
        assert record['stage'] == STAGE and record['versions']
        mapping = record['input_sha256']
        parameters = dict(conventions=record['conventions'],
            covariance_systems=[dict(N=r['N'], p=r['p']) for r in record['covariance_systems']],
            network_configurations=record['network_configurations'])
        origin = 'Original schema has no parameter digest; computed at sealing from conventions and configurations'
    elif name == SCIENCE[1]:
        assert record['status'] == 'passed' and record['environment']
        mapping = {record['source']['path']: record['source']['sha256']}
        parameters = record['parameters']
        origin = 'Original schema has no parameter digest; computed at sealing from saved parameters'
    else:
        assert name == SCIENCE[2] and record['status'] == 'passed' and record['environment']
        assert record['arithmetic'] == 'exact fractions.Fraction'
        mapping = record['source_hashes']
        assert mapping[record['source']['path']] == record['source']['sha256']
        parameters = record['parameters']
        origin = 'Original schema has no parameter digest; computed at sealing from saved parameters'
    for filename, expected in mapping.items():
        assert (ROOT/filename).is_file() and sha(ROOT/filename) == expected, (name, filename)
    return mapping, dict(origin=origin, sha256=digest(parameters), parameters=parameters)


def parents():
    output = {}
    for stage in PARENTS:
        path = ROOT/'results'/stage/'local-evidence.json'
        record = load(path)
        assert isinstance(record['hashes'], dict) and record['hashes']
        checked = []
        for name, expected in record['hashes'].items():
            if name in MUTABLE:
                continue
            assert (ROOT/name).is_file() and sha(ROOT/name) == expected, (stage, name)
            checked.append(name)
        output[path.relative_to(ROOT).as_posix()] = dict(sha256=sha(path),
            unchanged_scientific_files=checked,
            excluded_mutable_plans=sorted(MUTABLE.intersection(record['hashes'])))
    return output


def covariance_bound(n, p, m, ti, tj, side):
    if not 0 <= m <= n or m < p:
        return 0
    if side == 'row':
        ti, tj = tj, ti
    return choose(m-len(tj), p-len(tj))*choose(n-m+p-len(ti), p-len(ti))


def covariance_formula(n, p, occupied, ti, tj, side):
    if side == 'row':
        ti, tj = tj, ti
    m = len(occupied)
    if m < p or not tj <= occupied:
        return 0
    required = tj | (ti & occupied)
    return choose(m-len(required), p-len(required))*choose(n-m+p-len(ti), p-len(ti))


def tagged(record):
    systems = record['covariance_systems']
    assert [(r['N'], r['p']) for r in systems] == [(n, p) for n in range(2, 6) for p in (1, 2)]
    counts = Counter()
    for r in systems:
        n, p = r['N'], r['p']
        tags = [frozenset(t) for size in range(3) for t in combinations(range(n), size)]
        assert r['labels'] == comb(n, p)**2 and r['tags_per_side'] == len(tags)
        assert r['comparisons'] == 2*(1 << n)*len(tags)**2
        assert r['exact_formula_max_abs_error'] == r['upper_bound_max_positive_excess'] == r['individual_covariance_max_offdiagonal'] == 0
        near(r['lambda_N'], p*p/n)
        near(r['physical_w_N_for_Js_1'], 4**p/(p*p/n)*(factorial(p)/n**p)**2)
        for ti, tj in product(tags, repeat=2):
            for side in ('column', 'row'):
                for state in range(1 << n):
                    occupied = frozenset(i for i in range(n) if state & (1 << i))
                    value = covariance_formula(n, p, occupied, ti, tj, side)
                    bound = covariance_bound(n, p, len(occupied), ti, tj, side)
                    assert 0 <= value <= bound
                    counts.update(covariance_comparisons=1,
                        overlapping_TI_TJ_comparisons=int(bool(ti & tj)),
                        empty_restricted_family_comparisons=int(len(ti) > p or len(tj) > p),
                        strict_covariance_bound_comparisons=int(value < bound))
    modes = ('plain', 'sandwich', 'tagged', 'union')
    configs = [(2, 1, k, mode) for k in (2, 3, 4) for mode in modes]
    configs += [(3, 1, 3, mode) for mode in modes]+[(4, 2, 2, mode) for mode in modes]
    assert record['network_configurations'] == [list(c) for c in configs]
    seen = defaultdict(set)
    maximum_ratio = 0.0
    for r in record['networks']:
        n, p, k, mode = config = (r['N'], r['p'], r['pairs'], r['mode'])
        assert config in configs
        pairing = tuple(tuple(edge) for edge in r['pairing'])
        assert len(pairing) == k and all(a < b for a, b in pairing)
        assert sorted(i for edge in pairing for i in edge) == list(range(2*k))
        assert pairing == tuple(sorted(pairing)) and pairing not in seen[config]
        seen[config].add(pairing)
        crossings = sum(a < c < b < d or c < a < d < b for (a,b),(c,d) in combinations(pairing, 2))
        assert r['crossings'] == crossings
        assert len(r['endpoint_constraints']) == k
        endpoint, tags_by_edge, sizes = {}, {}, []
        for index, (edge, constraint) in enumerate(zip(pairing, r['endpoint_constraints'])):
            left_transpose = bool(index % 2)
            assert constraint['edge'] == list(edge) and constraint['left_transpose'] == left_transpose
            expected_tags = [[[], []], [[], []]]
            if mode == 'tagged' and index == 0:
                expected_tags = [[[0], [0]], [[0], []]]
            if mode == 'union' and index == 0:
                expected_tags = [[[0], []], [[], [n-1]]]
            if mode == 'union' and index == 1:
                expected_tags = [[[], [0]], [[0], [0]]]
            assert constraint['endpoint_tags'] == expected_tags
            ti, tj = set(), set()
            for position, transposed, restriction in zip(edge, (left_transpose, not left_transpose), expected_tags):
                a, b = map(set, restriction)
                if transposed:
                    a, b = b, a
                ti.update(a); tj.update(b)
                endpoint[position] = (index, transposed)
            assert constraint['base_union_TI'] == sorted(ti) and constraint['base_union_TJ'] == sorted(tj)
            size = choose(n-len(ti), p-len(ti))*choose(n-len(tj), p-len(tj))
            assert size == constraint['domain_size']
            tags_by_edge[index] = (ti, tj); sizes.append(size)
        assert r['assignments'] == prod(sizes)
        sandwiches = ['I']*(2*k+1)
        if mode != 'plain':
            sandwiches[0], sandwiches[k] = 'c0_dagger', 'c0'
        assert r['sandwiches'] == sandwiches
        shifts = {'I': 0, 'c0': -1, 'c0_dagger': 1}
        active, permutations, peak = [], 0, 4**n
        if all(sizes):
            for position in reversed(range(2*k)):
                edge, _ = endpoint[position]
                if edge not in active:
                    active.append(edge)
                else:
                    permutations += int(active[-1] != edge)
                    active.remove(edge)
                peak = max(peak, 4**n*prod(sizes[e] for e in active))
            assert not active
        else:
            permutations, peak = 0, 0
        assert r['ancilla_permutations'] == permutations and r['peak_ancilla_tensor_entries'] == peak
        assert r['direct_vs_ancilla_max_abs_error'] == 0
        assert [s['m'] for s in r['sectors']] == list(range(n+1))
        for s in r['sectors']:
            m = s['m']; charge = m; active = set(); factors = []
            assert s['input_dimension'] == comb(n, m)
            for position in reversed(range(2*k)):
                edge, transposed = endpoint[position]
                operation = 'column' if edge not in active else 'row'
                side = ('row' if operation == 'column' else 'column') if transposed else operation
                ti, tj = tags_by_edge[edge]
                bound = covariance_bound(n, p, charge, ti, tj, side)
                factors.append(dict(position=position, charge=charge, operation=operation,
                                    base_covariance=side, bound=bound))
                active.symmetric_difference_update({edge})
                charge += shifts[sandwiches[position]]
            assert s['endpoint_factors'] == factors and charge == m and not active
            squared = prod(f['bound'] for f in factors)
            assert s['bound_squared_exact'] == squared
            assert s['exact_psd_pass'] is True and F(s['minimum_schur_pivot']) >= 0
            norm = s['operator_norm_float']
            assert norm >= 0 and norm <= sqrt(squared)*(1+2e-12)
            near(s['norm_over_bound'], norm/sqrt(squared) if squared else 0)
            assert abs(float(F(s['normalized_trace_exact']))) <= norm+2e-12
            maximum_ratio = max(maximum_ratio, s['norm_over_bound'])
        assert r['nonzero'] == any(s['operator_norm_float'] > 0 for s in r['sectors'])
        counts.update(networks=1, crossed_networks=int(bool(crossings)), assignments=r['assignments'],
            ancilla_permutations=permutations, sector_psd_checks=len(r['sectors']),
            nonzero_networks=int(r['nonzero']),
            nonzero_crossed_sandwich_networks=int(r['nonzero'] and bool(crossings) and mode != 'plain'))
    for config in configs:
        k = config[2]
        assert len(seen[config]) == factorial(2*k)//(2**k*factorial(k))
    counts.update(max_direct_ancilla_abs_error=0.0,
                  max_operator_norm_over_bound_float=maximum_ratio, failure_count=0)
    expected = dict(covariance_comparisons=42336, overlapping_TI_TJ_comparisons=19312,
        empty_restricted_family_comparisons=17736, strict_covariance_bound_comparisons=14080,
        networks=564, crossed_networks=452, assignments=94743, ancilla_permutations=565,
        sector_psd_checks=1776, nonzero_networks=282, nonzero_crossed_sandwich_networks=113,
        max_direct_ancilla_abs_error=0.0, max_operator_norm_over_bound_float=1.0, failure_count=0)
    assert dict(counts) == record['summary'] == expected
    return dict(counts=dict(counts),
        replay='Scalar binomial formulas, complete pairing coverage, tag unions, charge paths, saved bounds and ratios',
        saved_certification='Exact PSD flags and nonnegative Schur pivots are validated as saved; matrices are not saved and PSD is not rerun',
        scope='Finite CAR fixtures; uniform double-scaling and Schatten claims require the separately bound analytic documents')


def coherence(record):
    p = record['parameters']
    assert p == dict(N=[2,3,4,5,6], p=[1,2], a=[0,1,2],
        m='every a <= m <= N-a; a <= floor(N/2)', covariance_weight_w_N=1,
        arithmetic='fractions.Fraction', CAR_site_numbering='zero-based')
    assert len(record['identities']) == 3 and all(record['identities'])
    rows = record['fixtures']; counts = Counter(); expected_keys = set(); seen = set(); visits = 0
    for n in p['N']:
        for degree in p['p']:
            for a in range(min(2, n//2)+1):
                length = n-2*a
                names = ['constant']+(['one_site'] if length else [])
                if length >= 2:
                    names += ['mixed_projector', 'signed_polynomial']
                for m in range(a, n-a+1):
                    expected_keys.update((n, degree, a, m, name) for name in names)
                    visits += len(names)*comb(n, degree)**2*comb(length, m-a)
                    kappa = choose(m-a, degree)*choose(n-m-a+degree, degree)
                    if kappa:
                        dimension = comb(length, m-a)
                        counts.update(kernel_rows=dimension, kernel_entries=dimension**2,
                                      tagged_flip_bounds=dimension*length)
    for r in rows:
        n, degree, a, m, name = key = tuple(r[k] for k in ['N','p','a','m','function'])
        assert key in expected_keys and key not in seen
        seen.add(key); length = n-2*a; dimension = comb(length, m-a)
        cm = choose(m, degree)*choose(n-m+degree, degree)
        kappa = choose(m-a, degree)*choose(n-m-a+degree, degree)
        assert r['C_m'] == cm and r['kappa_a'] == kappa
        values = []
        for subset in combinations(range(length), m-a):
            x, y = int(0 in subset), int(1 in subset)
            values.append(dict(constant=1, one_site=x, mixed_projector=x*(1-y),
                               signed_polynomial=2*x-y+3*x*y-1)[name])
        assert r['input_nonzero_entries'] == sum(value != 0 for value in values)
        output = r['output_nonzero_entries']
        assert isinstance(output, int) and 0 <= output <= dimension
        if not cm or not kappa or not r['input_nonzero_entries']:
            assert output == 0
        if name == 'constant' and kappa:
            assert output == dimension
        counts.update(coherence_identities=1, charge_path_identities=int(bool(cm)),
            inactive_zero_cases=int(not cm), sparse_entries_compared=output*(1+bool(cm)),
            nonconstant_input_fixtures=int(len(set(values)) > 1), nonzero_outputs=int(output > 0))
    assert seen == expected_keys
    # Full gauge tables were not saved; bind the reported count and inspect its witness.
    assert record['counts']['varying_gauge_sectors'] == 14
    counts['varying_gauge_sectors'] = record['counts']['varying_gauge_sectors']
    expected = dict(coherence_identities=340, charge_path_identities=265, inactive_zero_cases=75,
        kernel_rows=265, kernel_entries=2609, tagged_flip_bounds=1305, sparse_entries_compared=2044,
        nonconstant_input_fixtures=132, varying_gauge_sectors=14, nonzero_outputs=196)
    assert dict(counts) == record['counts'] == expected
    witness = record['gauge_witness']
    assert witness == dict(N=4,p=1,a=1,m=2,A=[1],B=[3],common_gauges={'1':1,'4':-1})
    occupied_support = sum(1 << site for site in witness['A']+witness['B'])
    assert not set(witness['A']) & set(witness['B'])
    assert all(int(mask).bit_count() == witness['m']-witness['a'] and
               not int(mask) & occupied_support for mask in witness['common_gauges'])
    assert set(witness['common_gauges'].values()) == {-1, 1}
    assert record['estimate']['coherence_fixtures'] == len(rows)
    assert record['estimate']['direct_kraus_entry_visits_upper_bound'] == visits == 113141
    return dict(counts=dict(counts), direct_entry_visit_upper_bound=visits,
        replay='Complete fixture coverage, scalar coefficients, input support counts and saved output-count arithmetic',
        saved_certification='CAR channel equality, positive kernels and flip inequalities are original driver certifications; no CAR or kernel main rerun. Full gauge-variation count is bound as reported, with its saved witness checked.',
        scope=record['scope'])


def charged(record):
    p = record['parameters']; rows = record['records']; counts = record['counts']
    assert p['polynomial_fixtures'] == [[2,1],[3,1],[3,2],[4,2]]
    assert p['maximum_gaussian_degree'] == 6
    neutral_units = sum(comb(2*n, n) for n in range(2,6))
    assert p['neutral_basis_units'] == counts['neutral_matrix_units'] == neutral_units == 348
    assert p['superoperator_identities'] == 9*neutral_units == 3132
    assert [(r['N'],r['p'],r['degree']) for r in rows] == [
        (n,degree,j) for n,degree in p['polynomial_fixtures'] for j in range(7)]
    for r in rows:
        assert isinstance(r['charge_difference_nonzero'], bool)
        assert r['reconstruction_residual_nonzero_entries'] == r['charged_residual_nonzero_entries'] == 0
        if r['degree'] < 2:
            assert r['charge_difference_nonzero'] is False
    assert counts['polynomial_identities'] == 3*len(rows) == 84
    assert counts['nonzero_charge_difference_fixtures'] == sum(r['charge_difference_nonzero'] for r in rows) == 16
    assert counts == dict(neutral_matrix_units=348, matrix_identities=3312,
        entry_comparisons=9949, polynomial_identities=84, nonzero_charge_difference_fixtures=16)
    assert len(record['checks']) == 7 and all(record['checks'])
    assert 'one standard real Gaussian' in record['gaussian_fixture']
    assert 'without an ensemble covariance assumption' in record['jacobi_fixture']
    return dict(counts=counts, polynomial_fixtures=len(rows),
        replay='Neutral basis count, superoperator count, complete polynomial fixture coverage, zero saved residuals and nonzero charge-difference count',
        saved_certification='The 3312 matrix identities and 9949 sparse entry comparisons are source-bound original driver certificates; matrices and Gaussian polynomials are not rerun.',
        scope=record['scope'], gaussian_fixture=record['gaussian_fixture'],
        jacobi_fixture=record['jacobi_fixture'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--visual-review-complete', action='store_true')
    args = parser.parse_args()
    assert args.visual_review_complete, 'Finish final-page review before sealing'
    target = OUT/'local-evidence.json'
    if target.exists():
        raise FileExistsError(target)
    previous = parents()
    science = {name:load(OUT/name) for name in SCIENCE}
    verified = {name:inputs(record,name) for name,record in science.items()}
    checks = {name:check(science[name]) for name,check in zip(SCIENCE,[tagged,coherence,charged])}
    source = ROOT/'manuscript'/STAGE/'main.tex'
    pdf = ROOT/'output/pdf'/f'de-sitter-{STAGE}.pdf'
    build = ROOT/'tmp/pdfs'/f'{STAGE}-build'
    render = ROOT/'tmp/pdfs'/f'{STAGE}-render'
    log = read_text(OUT/'main.log'); count = page_count(log)
    assert count > 0 and not re.search(r'undefined|Rerun to get|Overfull|Fatal error',log)
    assert sha(pdf) == sha(build/'main.pdf') and sha(OUT/'main.log') == sha(build/'main.log')
    assert (build/'main.pdf').stat().st_mtime_ns >= source.stat().st_mtime_ns
    passes = [build/f'pass{i}.txt' for i in [1,2]]; pass_records = []
    for number,path in enumerate(passes,1):
        content = read_text(path)
        assert page_count(content) == count
        assert not re.search(r'Fatal error|! Emergency stop',content)
        if number == 2:
            assert not re.search(r'undefined|Rerun to get|Overfull',content)
        pass_records.append(dict(pass_number=number,path=path.relative_to(ROOT).as_posix(),sha256=sha(path)))
    assert passes[1].stat().st_mtime_ns >= passes[0].stat().st_mtime_ns >= source.stat().st_mtime_ns
    width = len(str(count))
    pages = [render/f'page-{i:0{width}d}.png' for i in range(1,count+1)]
    assert set(render.glob('page-*.png')) == set(pages)
    assert all(path.stat().st_mtime_ns >= pdf.stat().st_mtime_ns for path in pages)
    metadata = {path.relative_to(ROOT).as_posix():png_metadata(path) for path in pages}
    names = set(DOCS) | {'code/record_bilocal_fourth_tail.py','code/record_bilocal_reconstruction.py',
        'code/build_bilocal_fourth_tail.ps1',f'manuscript/{STAGE}/main.tex',
        f'output/pdf/de-sitter-{STAGE}.pdf',f'results/{STAGE}/main.log'}
    for mapping,_ in verified.values():
        names.update(mapping)
    names.update(f'results/{STAGE}/{name}' for name in SCIENCE)
    diagnostics = sorted(path for path in (OUT/'diagnostics').rglob('*') if path.is_file())
    files = [ROOT/name for name in sorted(names)]+diagnostics
    assert len(files) == len(set(files)) and all(path.is_file() for path in files)
    result = dict(status='Local research milestone; long-term de Sitter research remains open',
        created_at_utc=datetime.now(timezone.utc).isoformat(),
        scientific_status='Saved finite tagged-covariance, coherence-channel and charged-reassembly fixtures accompany the final analytic stage documents.',
        proof_scope='The scientific conclusion and its hypotheses are in the hash-bound main note, charged proof, diagram/operator notes, and LaTeX manuscript. This recorder verifies provenance and saved scalar arithmetic, not an analytic theorem or a double-scaling limit.',
        proof_documents=[name for name in DOCS if name.startswith('notes/')]+[f'manuscript/{STAGE}/main.tex'],
        parameters='Finite fixtures use standard CAR, unit covariance weight, full finite Fock space with occupation-sector restrictions, and the separately saved parameters. Analytic full-Fock/annealed norms and limits must be read with their stated hypotheses in the final proof.',
        previous_evidence=previous,checks=checks,
        sealed_scientific_output_hashes_verified=SCIENCE_SHA256,
        input_hashes_verified={name:value[0] for name,value in verified.items()},
        parameter_hashes_recorded={name:value[1] for name,value in verified.items()},
        evidence_scopes={name:record['scope'] for name,record in science.items()},
        recorded_driver_environments={name:record.get('environment',record.get('versions')) for name,record in science.items()},
        git_commit_created=False,git_push_performed=False,archive_created=False,external_reply_sent=False,
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_bilocal_fourth_tail.ps1',
            pages=count,passes=2,pass_records=pass_records,compiler_header=log.splitlines()[0],
            render_command=f'pdftoppm -r 95 -png output/pdf/de-sitter-{STAGE}.pdf tmp/pdfs/{STAGE}-render/page',
            final_render_filenames=[path.name for path in pages],render_metadata=metadata,
            final_pages_visually_inspected=list(range(1,count+1)),
            visual_review_attestation='Caller supplied --visual-review-complete after final-page inspection. Recorder checks hashes, page counts, 95dpi PNG metadata and freshness; it does not see the pages.'),
        diagnostics=dict(preserved_files=[path.relative_to(ROOT).as_posix() for path in diagnostics]),
        environment=dict(python=platform.python_version(),platform=platform.platform(),dependencies='Python standard library only'),
        process_state='Saved science/build/render artifacts only; recorder executes no scientific main, compiler, renderer or Git.',
        hashes={path.relative_to(ROOT).as_posix():sha(path) for path in files},
        render_hashes={path.relative_to(ROOT).as_posix():sha(path) for path in pages},
        build_pass_hashes={path.relative_to(ROOT).as_posix():sha(path) for path in passes})
    with target.open('x',encoding='utf-8') as stream:
        json.dump(result,stream,ensure_ascii=False,indent=2,allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(status='passed',pages=count,scientific_outputs=len(SCIENCE),current_files=len(files),
        final_renders=len(pages),final_build_passes=2,
        previous_scientific_files_verified=sum(len(r['unchanged_scientific_files']) for r in previous.values()))))


if __name__ == '__main__':
    main()
