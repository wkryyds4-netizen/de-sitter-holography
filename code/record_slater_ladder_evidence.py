"""Seal saved Slater-ladder evidence after final build, render and visual review.

This recorder verifies the saved science hash and inputs, then replays only
scalar/binomial/Dyck bookkeeping. It does not run a CAR or Wick driver, compiler,
renderer, Git command, or scientific main. Successful evidence is exclusive-create.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, exp, factorial, isclose, prod
from pathlib import Path
import json
import platform
import re

from record_bilocal_reconstruction import (
    sha, load, digest, read_text, finite_numbers, page_count, png_metadata,
)

ROOT = Path(__file__).resolve().parent.parent
STAGE = 'susskind-slater-ladder'
OUT = ROOT/'results'/STAGE
SCIENCE = 'exact-checks.json'
SCIENCE_SHA256 = 'c3ffab60baf7b8302dfe73312bd5a386d9a5b942a893bc34d4dc5ba956b18ac1'
PARAMETER_SHA256 = 'babd1b557801241b185ab2632ee270137dae0a2ee33d96a90eb8e7706c288d7e'
INPUT_SHA256 = {
    'code/check_slater_ladder.py': '78d9a7a7ae45b636de7efe53e175dedca97844d473cab675fba24677523b0787',
    'results/susskind-bilocal-fourth-tail/local-evidence.json':
        '9055c2f1249cec6e4b2347a639532dd2f39ea8bcb381a880d77dd748a58692c8',
    'results/susskind-bilocal-four-point/local-evidence.json':
        '85701857a8e6f0040f837168dbf32077716505b11d36389ef48f9be0efb8ed21',
}
PARENTS = ['susskind-bilocal-fourth-tail', 'susskind-bilocal-four-point']
MUTABLE = {'README.md', 'PLAN.md', 'NEXT.md'}
# Artifact list and nine-page final build confirmed by the stage owner.
CONFIG_FINALIZED = True
EXPECTED_PAGES = 9
DOCS = [
    f'notes/{STAGE}.md', f'notes/{STAGE}-proof.md',
    f'notes/{STAGE}-literature-interface.md',
    'notes/susskind-bilocal-polarization-interface.md',
    'notes/susskind-bilocal-kms-polarization.md',
    'literature/susskind-bilocal-polarization-sources.md',
    'literature/susskind-bilocal-canonical-sources.md',
    f'journal/2026-11-07-{STAGE}.md',
    'README.md', 'PROJECT.md', 'PLAN.md', 'NEXT.md',
]


def near(a, b):
    assert isclose(a, b, rel_tol=2e-12, abs_tol=2e-14), (a, b)


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def falling(n, k):
    return prod(n-i for i in range(k))


def inputs(record):
    assert sha(OUT/SCIENCE) == SCIENCE_SHA256, 'Saved scientific output changed'
    finite_numbers(record)
    assert record['status'] == 'passed' and record['failures'] == []
    assert record['scope'] and record['environment'] and record['created_utc']
    assert record['elapsed_seconds'] >= 0
    assert record['input_sha256'] == INPUT_SHA256
    for filename, expected in INPUT_SHA256.items():
        assert (ROOT/filename).is_file() and sha(ROOT/filename) == expected, filename
    assert record['parameter_sha256'] == digest(record['parameters']) == PARAMETER_SHA256
    return dict(origin='Verified against original saved parameter digest',
                sha256=PARAMETER_SHA256, parameters=record['parameters'])


def parents():
    result = {}
    for stage in PARENTS:
        path = ROOT/'results'/stage/'local-evidence.json'
        record = load(path)
        assert isinstance(record['hashes'], dict) and record['hashes']
        checked = []
        for filename, expected in record['hashes'].items():
            if filename in MUTABLE:
                continue
            assert (ROOT/filename).is_file() and sha(ROOT/filename) == expected, (stage, filename)
            checked.append(filename)
        result[path.relative_to(ROOT).as_posix()] = dict(
            sha256=sha(path), unchanged_scientific_files=checked,
            excluded_mutable_plans=sorted(MUTABLE.intersection(record['hashes'])))
    return result


def dyck_data(word):
    height, multiplicity, births = 0, 1, []
    for sign in reversed(word):
        if sign == '+':
            births.append(height)
            height += 1
        else:
            assert sign == '-'
            if height == 0:
                return None
            multiplicity *= height
            height -= 1
    assert height == 0
    return births, multiplicity


def fixture_checks(record):
    p = record['parameters']
    configurations = [[4, 1, 4], [6, 1, 3], [6, 2, 3], [8, 2, 3]]
    assert p['fixtures_N_p_maximum_pairs'] == configurations
    assert p['filling'] == 'm=N/2'
    assert p['CAR'] == '{c_i,c_j^dagger}=delta_ij'
    assert p['normalization'] == 'unit Gaussian variance, sigma_squared=L=binom(m,p)binom(N-m,p)'
    assert p['creators'] == 'ascending operator order'
    assert p['word_order'] == 'left to right product; rightmost acts first; + creates a particle-hole block'
    assert p['expected_maximum_label_assignments'] == 1600000
    assert p['scope'] == record['scope']
    fixtures = record['fixtures']
    assert [(r['N'], r['p']) for r in fixtures] == [(n, degree) for n, degree, _ in configurations]
    totals, witness_count = Counter(), 0
    for fixture, (n, degree, maximum) in zip(fixtures, configurations):
        m = n//2
        a = n-m
        labels = choose(m, degree)*choose(a, degree)
        assert fixture['m'] == m and fixture['L'] == labels
        expected_words = {}
        for k in range(1, maximum+1):
            for plus_positions in combinations(range(2*k), k):
                plus = set(plus_positions)
                word = ''.join('+' if i in plus else '-' for i in range(2*k))
                data = dyck_data(word)
                if data is not None:
                    expected_words[word] = data
        rows = fixture['records']
        assert len(rows) == len(expected_words)
        assert {r['word'] for r in rows} == set(expected_words)
        counts = Counter(
            non_Dyck_words_zero_by_right_action=sum(
                comb(2*k, k)-comb(2*k, k)//(k+1) for k in range(1, maximum+1)),
            Gaussian_label_assignments=sum(
                comb(2*k, k)//(k+1)*factorial(k)*labels**k for k in range(1, maximum+1)),
            forward_CAR_interval_comparisons=0,
            word_moments=len(rows), nonzero_backward_word_corrections=0)
        for row in rows:
            word = row['word']
            births, multiplicity = expected_words[word]
            k = len(word)//2
            assert row['creation_heights'] == births
            assert row['forward_pairings'] == multiplicity
            assert row['backward_pairings'] == factorial(k)-multiplicity
            forward = F(multiplicity)*prod(
                F(choose(a-h*degree, degree)*choose(m-h*degree, degree), labels)
                for h in births)
            assert F(row['forward_moment']) == forward
            backward = F(row['backward_correction'])
            assert F(row['full_moment']) == forward+backward
            assert F(row['full_moment']) >= 0
            if row['backward_pairings'] == 0:
                assert backward == 0
            assert row['limiting_shift_moment'] == f'{multiplicity} * q^{sum(births)}'
            counts.update(forward_CAR_interval_comparisons=multiplicity*labels**k,
                          nonzero_backward_word_corrections=int(backward != 0))
        assert dict(counts) == fixture['counts']
        totals.update(counts)
        witness = fixture['nonzero_backward_pair_witness']
        assert witness == dict(word='-+-+', pairing=[[1, 2], [3, 0]],
                               label_assignment=[0, 0], CAR_amplitude=1)
        assert any(minus > plus for plus, minus in witness['pairing'])
        assert sorted(pos for pair in witness['pairing'] for pos in pair) == list(range(4))
        assert all(witness['word'][plus] == '+' and witness['word'][minus] == '-'
                   for plus, minus in witness['pairing'])
        witness_count += 1
        first = fixture['full_H_first_step']
        cm = choose(m, degree)*choose(a+degree, degree)
        weights = {str(k):choose(m, degree)*choose(a, k)*choose(degree, k)
                   for k in range(degree+1)}
        captured = F(labels, cm)
        assert first['C_m_at_w_1'] == cm and first['sigma_squared_at_w_1'] == labels
        assert first['weights_by_created_pairs'] == weights
        assert sum(weights.values()) == cm and weights[str(degree)] == labels
        assert F(first['captured_weight']) == F(first['normalized_overlap_squared']) == captured
        assert F(first['leakage_fraction']) == 1-captured
        assert first['H_Vdagger_overlap_at_w_1'] == labels
    totals['finite_systems'] = len(fixtures)
    assert totals['Gaussian_label_assignments'] <= p['expected_maximum_label_assignments']
    return totals, dict(
        replay='Complete right-action Dyck-word coverage, scalar forward moments, pairing/count arithmetic and full-H first-step hypergeometric weights.',
        saved_certification='Full CAR/Wick moments, finite backward corrections and the nonzero backward-pair CAR witness are bound original driver outputs; CAR actions and label assignments are not rerun.',
        saved_backward_witnesses=witness_count)


def asymptotic_checks(record):
    p = record['parameters']['asymptotic_parameters']
    assert p == dict(lambda_N='1/2', p_values=[2, 4, 8, 16, 32, 64], maximum_r=8)
    rows = record['hypergeometric_asymptotic']
    assert [r['p'] for r in rows] == p['p_values']
    for row in rows:
        degree = row['p']
        n = 2*degree**2
        a = n//2
        assert row['N'] == n and F(row['lambda_N']) == F(degree**2, n) == F(1, 2)
        denominator = choose(a+degree, degree)
        distribution = [F(choose(degree, r)*choose(a, degree-r), denominator)
                        for r in range(degree+1)]
        assert sum(distribution) == 1
        assert [F(value) for value in row['saved_r_probabilities']] == distribution[:p['maximum_r']+1]
        assert F(row['captured_weight']) == distribution[0]
        near(row['captured_weight_decimal'], float(distribution[0]))
        near(row['limit_captured_weight'], exp(-1))
        assert F(row['exact_mean_deficit']) == F(degree**2, a+degree)
        assert [F(value) for value in row['factorial_moments']] == [
            F(falling(degree, j)**2, falling(a+degree, j)) for j in range(1, min(4, degree)+1)]
        poisson = [exp(-1)/factorial(r) for r in range(degree+1)]
        distance = .5*(sum(abs(float(value)-poisson[r]) for r, value in enumerate(distribution))
                       +1-sum(poisson))
        near(row['total_variation_to_Poisson1'], distance)
        assert 0 <= row['total_variation_to_Poisson1'] <= 1
    return len(rows)


def cumulant_checks(record):
    expected = {(q, coefficients) for q in [F(1, 4), F(1, 2), F(3, 4)]
                for coefficients in product([-1, 0, 1, 2], repeat=3)
                if coefficients != (0, 0, 0)}
    seen = set()
    for row in record['linear_cumulants']:
        q = F(row['q'])
        coefficients = tuple(F(value) for value in row['coefficients'])
        key = q, coefficients
        assert key in expected and key not in seen
        seen.add(key)
        a, b, c = coefficients
        aa, bb = a*a+b*b, c*c
        variance = aa+bb
        cumulant = -(1-q)*aa**2-2*(1-q*q)*aa*bb-(1-q**4)*bb**2
        assert F(row['variance']) == variance > 0
        assert F(row['fourth_cumulant']) == cumulant < 0
        assert F(row['fourth_moment']) == 3*variance**2+cumulant
        assert (q**4-1)*variance**2 <= cumulant <= (q-1)*variance**2
    assert seen == expected
    return len(seen)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--visual-review-complete', action='store_true')
    args = parser.parse_args()
    assert CONFIG_FINALIZED and isinstance(EXPECTED_PAGES, int) and EXPECTED_PAGES > 0, 'Final artifact list/page count pending'
    assert args.visual_review_complete, 'Finish final-page visual review before sealing'
    target = OUT/'local-evidence.json'
    if target.exists():
        raise FileExistsError(target)
    science = load(OUT/SCIENCE)
    parameter_hash = inputs(science)
    previous = parents()
    counts, checks = fixture_checks(science)
    counts['asymptotic_rows'] = asymptotic_checks(science)
    counts['linear_cumulant_checks'] = cumulant_checks(science)
    expected_counts = dict(non_Dyck_words_zero_by_right_action=136,
        Gaussian_label_assignments=1537310, forward_CAR_interval_comparisons=754030,
        word_moments=46, nonzero_backward_word_corrections=25, finite_systems=4,
        asymptotic_rows=6, linear_cumulant_checks=189)
    assert dict(counts) == science['counts'] == expected_counts
    checks['counts'] = dict(counts)
    checks['asymptotic_replay'] = 'Six saved scalar hypergeometric distributions, factorial moments and finite total-variation diagnostics; not a proof of convergence.'
    checks['cumulant_replay'] = 'All 189 saved rational linear-combination variances and fourth cumulants, including both analytic bounds.'
    source = ROOT/'manuscript'/STAGE/'main.tex'
    pdf = ROOT/'output/pdf'/f'de-sitter-{STAGE}.pdf'
    build = ROOT/'tmp/pdfs'/f'{STAGE}-build'
    render = ROOT/'tmp/pdfs'/f'{STAGE}-render'
    log = read_text(OUT/'main.log')
    count = page_count(log)
    assert count == EXPECTED_PAGES
    assert not re.search(r'undefined|Rerun to get|Overfull|Fatal error', log)
    assert sha(pdf) == sha(build/'main.pdf')
    assert sha(OUT/'main.log') == sha(build/'main.log')
    assert (build/'main.pdf').stat().st_mtime_ns >= source.stat().st_mtime_ns
    passes = [build/f'pass{i}.txt' for i in [1, 2]]
    pass_records = []
    for number, path in enumerate(passes, 1):
        content = read_text(path)
        assert page_count(content) == count
        assert not re.search(r'Fatal error|! Emergency stop', content)
        if number == 2:
            assert not re.search(r'undefined|Rerun to get|Overfull', content)
        pass_records.append(dict(pass_number=number, path=path.relative_to(ROOT).as_posix(), sha256=sha(path)))
    assert passes[1].stat().st_mtime_ns >= passes[0].stat().st_mtime_ns >= source.stat().st_mtime_ns
    width = len(str(count))
    pages = [render/f'page-{i:0{width}d}.png' for i in range(1, count+1)]
    assert set(render.glob('page-*.png')) == set(pages)
    assert all(path.stat().st_mtime_ns >= pdf.stat().st_mtime_ns for path in pages)
    metadata = {path.relative_to(ROOT).as_posix():png_metadata(path) for path in pages}
    names = set(DOCS) | set(INPUT_SHA256) | {
        'code/record_slater_ladder_evidence.py', 'code/record_bilocal_reconstruction.py',
        'code/build_slater_ladder.ps1', f'manuscript/{STAGE}/main.tex',
        f'output/pdf/de-sitter-{STAGE}.pdf', f'results/{STAGE}/main.log',
        f'results/{STAGE}/{SCIENCE}'}
    diagnostics = sorted(path for path in (OUT/'diagnostics').rglob('*') if path.is_file())
    files = [ROOT/name for name in sorted(names)]+diagnostics
    assert len(files) == len(set(files)) and all(path.is_file() for path in files)
    result = dict(
        status='Local research milestone; long-term de Sitter research remains open',
        created_at_utc=datetime.now(timezone.utc).isoformat(),
        scientific_status='Saved finite CAR/Wick fixtures accompany the fixed-word pure Slater-ladder proof, its full-H first-step return decomposition, and the literature interface.',
        proof_scope='Analytic conclusions and hypotheses are in the hash-bound stage documents. This recorder verifies provenance, saved scalar arithmetic and finite coverage; it does not prove an analytic limit or rerun CAR/Wick science.',
        proof_documents=[name for name in DOCS if name.startswith('notes/')]+[f'manuscript/{STAGE}/main.tex'],
        parameters='Standard complex CAR, m=N/2, unit Gaussian variance in finite fixtures, reference Slater state and saved normalization. Analytic fixed-word double-scaling and self-averaging scopes are given separately in the proof.',
        previous_evidence=previous, checks=checks,
        sealed_scientific_output_hashes_verified={SCIENCE:SCIENCE_SHA256},
        input_hashes_verified={SCIENCE:INPUT_SHA256},
        parameter_hashes_recorded={SCIENCE:parameter_hash},
        evidence_scopes={SCIENCE:science['scope']},
        recorded_driver_environments={SCIENCE:science['environment']},
        git_commit_created=False, git_push_performed=False,
        archive_created=False, external_reply_sent=False,
        build=dict(
            command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_slater_ladder.ps1',
            pages=count, passes=2, pass_records=pass_records, compiler_header=log.splitlines()[0],
            render_command=f'pdftoppm -r 95 -png output/pdf/de-sitter-{STAGE}.pdf tmp/pdfs/{STAGE}-render/page',
            final_render_filenames=[path.name for path in pages], render_metadata=metadata,
            final_pages_visually_inspected=list(range(1, count+1)),
            visual_review_attestation='Caller supplied --visual-review-complete after final-page inspection. Recorder checks hashes, page counts, 95dpi PNG metadata and freshness; it does not see the pages.'),
        diagnostics=dict(preserved_files=[path.relative_to(ROOT).as_posix() for path in diagnostics]),
        environment=dict(python=platform.python_version(), platform=platform.platform(), dependencies='Python standard library only'),
        process_state='Saved science/build/render artifacts only; recorder executes no scientific main, compiler, renderer or Git.',
        hashes={path.relative_to(ROOT).as_posix():sha(path) for path in files},
        render_hashes={path.relative_to(ROOT).as_posix():sha(path) for path in pages},
        build_pass_hashes={path.relative_to(ROOT).as_posix():sha(path) for path in passes})
    with target.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(status='passed', pages=count, scientific_outputs=1,
        current_files=len(files), final_renders=len(pages), final_build_passes=2,
        previous_scientific_files_verified=sum(len(r['unchanged_scientific_files']) for r in previous.values()))))


if __name__ == '__main__':
    main()
