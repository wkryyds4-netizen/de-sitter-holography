"""Read-only evidence seal; run only after final build/render/visual review.

No science driver, compiler, renderer or Git command is executed. The four
saved schemas are deliberately handled separately. Output is exclusive-create.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction as F
from math import comb, exp, factorial, isclose, prod, sqrt
from pathlib import Path
import json
import platform
import re
import struct

from record_bilocal_reconstruction import sha, load, digest, read_text, finite_numbers, page_count, png_metadata

ROOT = Path(__file__).resolve().parent.parent
STAGE = 'susskind-bilocal-four-point'
OUT = ROOT/'results'/STAGE
SCIENCE = ['finite-four-source-checks.json', 'joint-moment-checks.json',
           'transfer-checks.json', 'mixed-otoc-checks.json']
PARENTS = ['susskind-bilocal-reconstruction', 'susskind-singlet-multimark', 'susskind-memory-size-bridge']
MUTABLE = {'README.md', 'PLAN.md', 'NEXT.md'}
DOCS = [f'notes/{STAGE}.md', 'notes/susskind-bilocal-multisource-proof.md',
        'notes/susskind-bilocal-four-point-transfer.md', f'literature/{STAGE}-sources.md',
        f'journal/2026-11-07-{STAGE}.md', 'README.md', 'PROJECT.md', 'PLAN.md', 'NEXT.md']


def near(a, b):
    assert isclose(a, b, rel_tol=2e-12, abs_tol=1e-25), (a, b)


def inputs(record, name):
    assert record['status'] == 'passed' and record['scope']; finite_numbers(record)
    if name == 'transfer-checks.json':
        assert record['schema'] == 1 and record['arithmetic'] and record['python']
        mapping = {'code/bilocal_four_point_transfer_check.py': record['source_sha256']}
        parameter_record = dict(origin='Original schema records deterministic parameters within cases',
                                sha256=digest(record['cases']))
    else:
        assert record['schema_version'] == 1 and record['environment']
        assert digest(record['parameters']) == record['parameter_sha256']
        mapping = record['input_sha256']
        parameter_record = dict(origin='Verified saved parameter digest', sha256=record['parameter_sha256'])
    for filename, expected in mapping.items():
        assert (ROOT/filename).is_file() and sha(ROOT/filename) == expected, (name, filename)
    return mapping, parameter_record


def parents():
    output = {}
    for stage in PARENTS:
        path = ROOT/'results'/stage/'local-evidence.json'; record = load(path); checked = []
        assert isinstance(record['hashes'], dict) and record['hashes']
        for name, expected in record['hashes'].items():
            if name in MUTABLE:
                continue
            assert (ROOT/name).is_file() and sha(ROOT/name) == expected, (stage, name)
            checked.append(name)
        output[path.relative_to(ROOT).as_posix()] = dict(sha256=sha(path), unchanged_scientific_files=checked,
            excluded_mutable_plans=sorted(MUTABLE.intersection(record['hashes'])))
    return output


def four_source(record):
    p = record['parameters']; totals = Counter(); geometry = Counter()
    assert p['systems'] == [[2, 1], [3, 1]] and p['all_source_masks'] == list(range(16))
    assert p['double_slot_kinds'] == ['DDDD', 'LLLL', 'DLDL', 'LDLD']
    assert p['formal_lambda_for_partition_check'] == '1/4'
    assert [(r['N'], r['p'], r['mode']) for r in record['finite']] == [(2,1,'double'),(2,1,'single'),(3,1,'double'),(3,1,'single')]
    for r in record['finite']:
        counts = r['counts']; double = r['mode'] == 'double'
        assert r['assignments'] == (4 if double else r['N']**4)
        assert counts['pairings'] == (60 if double else 3)
        assert counts['raw_cases'] == r['assignments']*counts['pairings']*(4 if double else 1)
        assert counts['raw_coefficients'] == 16*counts['raw_cases']
        if not double:
            assert counts['W1_cases'] == counts['raw_cases'] and counts['W1_coefficients'] == counts['raw_coefficients']
        for kind in (['raw'] if double else ['raw', 'W1']):
            assert 0 <= counts[kind+'_negative'] <= counts[kind+'_coefficients']-counts[kind+'_zero']
            assert 0 < counts[kind+'_four_mark_nonzero'] <= counts[kind+'_cases']
        totals.update(counts)
    assert totals['raw_cases'] == 2211 and totals['raw_coefficients'] == 35376
    assert totals['W1_cases'] == 291 and totals['W1_coefficients'] == 4656
    assert totals['raw_negative'] == 68 and totals['W1_negative'] == 992
    assert len(record['finite_nonzero_examples']) == 6
    for r in record['finite_nonzero_examples']:
        assert F(r['coefficients']['15']) != 0
        assert all(0 <= int(mask) < 16 for mask in r['coefficients'])
    lam = F(p['formal_lambda_for_partition_check'])
    for r in record['geometry']:
        kappa = {int(mask): F(v) for mask, v in r['kappa'].items()}
        assert all(kappa.get(mask, 0) == 0 for mask in range(16) if mask.bit_count() >= 3)
        # Replay the saved nilpotent exponent, not the CAR/Poisson computation.
        value = {0: F(1)}; power = {0: F(1)}
        reduced = {mask: lam*v for mask, v in kappa.items() if mask}
        for degree in range(1, 5):
            next_power = Counter()
            for a, av in power.items():
                for b, bv in reduced.items():
                    if not a & b: next_power[a | b] += av*bv
            power = dict(next_power)
            for mask, coefficient in power.items():
                value[mask] = value.get(mask, F(0))+coefficient/factorial(degree)
        assert value.get(15, F(0)) == F(r['four_mark_reduced_moment'])
        geometry.update(exponent_cases=1, all_mask_comparisons=16, higher_connected_zero_coefficients=5,
                        partition_comparisons=16, four_mark_moment_nonzero=int(value.get(15, 0) != 0))
        if r['mode'] == 'mixed':
            assert kappa.get(1, 0) == kappa.get(4, 0) == 1
            assert all(kappa.get(mask, 0) == 0 for mask in [3, 9, 6, 12])
            assert kappa.get(5, 0) == int([0, 3] in r['pairing'])
            geometry.update(single_L_double_D_zero_checks=4, single_L_pair_checks=1)
    assert dict(geometry) == record['geometry_counts'] == dict(exponent_cases=298, all_mask_comparisons=4768,
        higher_connected_zero_coefficients=1490, partition_comparisons=4768, four_mark_moment_nonzero=216,
        single_L_double_D_zero_checks=40, single_L_pair_checks=10)
    assert record['maximum_absolute_integer'] == 81 < 2**53
    return dict(finite_counts=dict(totals), geometry_counts=dict(geometry),
                arithmetic_replay='Saved nilpotent exponent to four-mark coefficient; no CAR main rerun')


def joint(record):
    p = record['parameters']; rows = record['records']; stats = Counter()
    assert p['C'] == 1 and p['cyclic_representatives'] is True
    assert p['formal_points'] == [['1/3','1/4'], ['1/2','1/2'], ['3/4','1/8']]
    expected_counts = {tuple(counts) for counts in p['counts_J_W_H']}
    for r in rows:
        word = r['word']; assert set(word) <= set('JWH')
        assert tuple(word.count(letter) for letter in 'JWH') in expected_counts
        orbit = {word[i:]+word[:i] for i in range(len(word))}
        assert word == min(orbit) and len(orbit) == r['cyclic_orbit_size']
        assert r['denominator_power'] == word.count('J')
        assert r['exact'] is True and r['actual_numerator'] == r['expected_numerator']
        assert r['zero_moment'] == (not r['actual_numerator'])
        assert all(F(term['coefficient']) != 0 and term['lambda_degree'] >= 0 and term['q_degree'] >= 0
                   for term in r['actual_numerator'])
        stats.update({k:v for k,v in r['statistics'].items() if k != 'maximum_source_order'})
    counts = dict(cyclic_representatives=len(rows), words_covered=sum(r['cyclic_orbit_size'] for r in rows),
        zero_representatives=sum(r['zero_moment'] for r in rows), old_JH_interface_checks=sum('W' not in r['word'] for r in rows),
        maximum_source_order_observed=max(r['statistics']['maximum_source_order'] for r in rows),
        formal_Gram_diagnostics=len(record['diagnostics']), **stats)
    assert counts == record['counts'] == dict(cyclic_representatives=53, words_covered=233, zero_representatives=23,
        old_JH_interface_checks=4, maximum_source_order_observed=2, formal_Gram_diagnostics=3,
        accepted_diagrams=1763, negative_diagram_coefficients=110)
    for r in record['diagnostics']:
        q, lam = F(r['q']), F(r['lambda_value']); w = 4*lam+16*lam*lam*q/(1+q)
        v = {k:F(x) for k,x in r['moments'].items()}
        assert v == dict(HH=F(1), WW=lam, JJ=w, JW=F(0), WH=F(0), JJWW=w*lam,
                         JWJW=q*q*w*lam, WWWW=(2+q)*lam*lam, JJJJ=(2+q**4)*w*w,
                         WHWH=q*lam, WWHH=lam)
        assert F(r['commutator_squared_JW']) == 2*(v['JJWW']-v['JWJW']) > 0
        assert F(r['product_Gram_determinant_JW']) == v['JJWW']**2-v['JWJW']**2 > 0
        assert F(r['normalized_commutator_squared_JW']) == 2*(1-q*q)
    return dict(counts=counts, scope='Recorded formal common-source identities and exact Gram arithmetic')


def transfer(record):
    rows = record['cases']; assert [(r['N'],r['p']) for r in rows] == [(2,1),(3,2),(4,2)]
    for r in rows:
        n, p, d = r['N'], r['p'], r['dimension']; q = F(r['formal_q']); lam = F(r['formal_lambda'])
        assert d == 2**n and lam == F(p*p, n) and q == F(1,2) and F(r['target_C']) == 1
        assert F(r['kappa']) == 2*lam*(1-q)/(1+q)
        assert len(r['coupling']) == comb(n,p)
        assert all(len(row) == comb(n,p) for row in r['coupling'])
        assert all(r['coupling'][i][j] == r['coupling'][j][i] for i in range(comb(n,p)) for j in range(comb(n,p)))
        assert r['matrix_identity_count'] == 11 and r['scalar_entry_comparisons'] == 9*d*d+2
        assert F(r['commutator_norm_squared']) == 2*(F(r['J2_squared_W1_squared'])-F(r['J2_W1_J2_W1'])) >= 0
    counts = dict(cases=len(rows), matrix_identities=sum(r['matrix_identity_count'] for r in rows),
        scalar_entry_comparisons=sum(r['scalar_entry_comparisons'] for r in rows),
        nonzero_commutators=sum(F(r['commutator_norm_squared']) > 0 for r in rows))
    assert counts == dict(cases=3, matrix_identities=record['total_matrix_identities'],
        scalar_entry_comparisons=record['total_scalar_entry_comparisons'], nonzero_commutators=record['nonzero_commutators'])
    assert counts == dict(cases=3, matrix_identities=33, scalar_entry_comparisons=3030, nonzero_commutators=1)
    return dict(counts=counts, scope='Exact finite conserved-difference identities; no eigenvector approximation or DS fitting')


def mixed_otoc(record):
    p = record['parameters']; rows = record['records']; coefficients = record['coefficients']
    assert p['Js'] == 1 and p['physical_lambda'] == [.25,.1,.05]
    assert p['probe_lengths'] == [2,1] and p['swapped_probe_lengths'] == [1,2]
    assert p['times'] == sorted({i/20 for i in range(21)} | {.01})
    assert p['requested_truncation_tolerance'] == 1e-11 and p['refinement_increment'] == 16
    assert [r['order'] for r in coefficients] == p['Taylor_orders'] == [0,2,4,6]
    for r in coefficients:
        n = r['order']; pairings = factorial(n)//(2**(n//2)*factorial(n//2))
        assert r['species_swap_exact'] is True and r['arc_layouts'] == comb(n+3,3)
        assert r['colored_pairings'] == r['arc_layouts']*pairings
        assert all(F(v) != 0 for v in r['coefficients'].values())
    assert coefficients[0]['coefficients'] == {'2':'1'}
    points = [point for row in rows for point in row['points']]
    counts = dict(Taylor_coefficients=4, exact_species_swap_coefficients=4,
        arc_layouts_per_species_order=sum(r['arc_layouts'] for r in coefficients),
        colored_pairings_per_species_order=sum(r['colored_pairings'] for r in coefficients),
        physical_lambda_cases=len(rows), time_points=len(points), kernel_evaluations=3*len(points),
        independent_small_time_comparisons=sum('chord_series_absolute_difference' in x for x in points))
    assert counts == record['counts'] == dict(Taylor_coefficients=4, exact_species_swap_coefficients=4,
        arc_layouts_per_species_order=130, colored_pairings_per_species_order=1376,
        physical_lambda_cases=3, time_points=66, kernel_evaluations=198, independent_small_time_comparisons=6)
    for r in rows:
        lam = r['lambda_value']; q = r['q']; near(q, exp(-4*lam)); near(r['sigma'], 1/sqrt(lam))
        assert r['cutoff'] == max(x['L'] for x in r['source_cutoff_requests'])
        assert r['refined_cutoff'] == r['cutoff']+16
        assert [x['time'] for x in r['points']] == p['times']
        for b in r['interval_bounds']:
            d1,d2 = b['lengths']; a = q**d1; mutual = q**(d1*d2)
            near(b['Cb'],1/prod(1-q**j for j in range(1,2*d2)))
            near(b['K'],b['Cb']/(1-a)); near(b['source_tail'],q**(d2*b['inner']))
            value = 2*mutual*b['K']*(b['vacuum_bound']+b['propagator_bound']+2*b['source_tail'])
            near(value,b['otoc_absolute_truncation_bound']); near(2*value,b['commutator_absolute_truncation_bound'])
            assert 0 <= value < p['requested_truncation_tolerance']
        for x in r['points']:
            t, f, d = x['time'], x['otoc'], x['commutator_squared']; assert d == 2*(1-f)
            assert 0 <= f <= 1 and 0 <= d <= 2
            c = abs(x['coarse_otoc']-f); s = abs(x['coarse_otoc']-x['swapped_coarse_otoc'])
            near(c,x['cutoff_absolute_difference']); near(s,x['species_swap_absolute_difference'])
            assert max(c,s) < p['floating_comparison_tolerance']
            near(c/f,x['cutoff_relative_otoc_difference']); near(s/f,x['species_swap_relative_otoc_difference'])
            near(2*max(c,s),x['empirical_absolute_commutator_difference'])
            near(2*max(c,s)/d,x['empirical_relative_commutator_difference'])
            near((1+abs(f))/abs(1-f),x['commutator_subtraction_condition'])
            tail = x['otoc_absolute_analytic_truncation_bound']; assert 0 <= tail < p['requested_truncation_tolerance']
            near(tail/f,x['otoc_relative_analytic_truncation_bound'])
            near(2*tail,x['commutator_absolute_analytic_truncation_bound'])
            near(2*tail/d,x['commutator_relative_analytic_truncation_bound'])
            if t in p['Taylor_comparison_times']:
                terms = [sum(float(F(v))*q**int(k) for k,v in c['coefficients'].items())*(t/sqrt(lam))**c['order'] for c in coefficients]
                assert terms == x['chord_series_terms']; series = sum(terms)
                near(series,x['sixth_order_chord_series']); near(abs(f-series),x['chord_series_absolute_difference'])
                assert abs(f-series) < p['Taylor_comparison_tolerance']
            if t == 0: near(abs(f-q*q),x['equal_time_mutual_factor_error'])
    for name,key in [('maximum_cutoff_absolute_difference','cutoff_absolute_difference'),
                     ('maximum_species_swap_absolute_difference','species_swap_absolute_difference'),
                     ('maximum_commutator_subtraction_condition','commutator_subtraction_condition')]:
        assert record[name] == max(x[key] for x in points)
    figure = ROOT/record['figure']['path']; assert sha(figure) == record['figure']['sha256']
    assert record['figure']['dpi'] == p['figure_dpi'] == 180
    return dict(counts=counts, figure_sha256=sha(figure),
        maximum_cutoff_difference=record['maximum_cutoff_absolute_difference'],
        maximum_species_swap_difference=record['maximum_species_swap_absolute_difference'],
        scope='Existing unequal-probe kernel applied to normalized J2/W1 with mutual q^2. Analytic '
              'bounds cover Jacobi truncation; saved floating differences are not interval certification.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--visual-review-complete', action='store_true')
    args = parser.parse_args(); assert args.visual_review_complete, 'Finish final-page review before sealing'
    target = OUT/'local-evidence.json'
    if target.exists(): raise FileExistsError(target)
    previous = parents(); science = {name:load(OUT/name) for name in SCIENCE}
    verified = {name:inputs(record,name) for name,record in science.items()}
    checks = {name:check(science[name]) for name,check in zip(SCIENCE,[four_source,joint,transfer,mixed_otoc])}
    source = ROOT/'manuscript'/STAGE/'main.tex'; pdf = ROOT/'output/pdf'/f'de-sitter-{STAGE}.pdf'
    build = ROOT/'tmp/pdfs'/f'{STAGE}-build'; render = ROOT/'tmp/pdfs'/f'{STAGE}-render'
    log = read_text(OUT/'main.log'); count = page_count(log)
    assert count > 0 and not re.search(r'undefined|Rerun to get|Overfull|Fatal error',log)
    assert sha(pdf) == sha(build/'main.pdf') and sha(OUT/'main.log') == sha(build/'main.log')
    assert (build/'main.pdf').stat().st_mtime_ns >= source.stat().st_mtime_ns
    passes = [build/f'pass{i}.txt' for i in [1,2]]; pass_records = []
    for number,path in enumerate(passes,1):
        content = read_text(path); assert page_count(content) == count
        assert not re.search(r'Fatal error|! Emergency stop',content)
        if number == 2: assert not re.search(r'undefined|Rerun to get|Overfull',content)
        pass_records.append(dict(pass_number=number,path=path.relative_to(ROOT).as_posix(),sha256=sha(path)))
    assert passes[1].stat().st_mtime_ns >= passes[0].stat().st_mtime_ns >= source.stat().st_mtime_ns
    width = len(str(count)); pages = [render/f'page-{i:0{width}d}.png' for i in range(1,count+1)]
    assert set(render.glob('page-*.png')) == set(pages)
    assert all(path.stat().st_mtime_ns >= pdf.stat().st_mtime_ns for path in pages)
    metadata = {path.relative_to(ROOT).as_posix():png_metadata(path) for path in pages}
    names = set(DOCS) | {'code/record_bilocal_four_point.py','code/record_bilocal_reconstruction.py',
        'code/build_bilocal_four_point.ps1',f'manuscript/{STAGE}/main.tex',
        f'output/pdf/de-sitter-{STAGE}.pdf',f'results/{STAGE}/main.log',f'results/{STAGE}/normalized-commutator.png'}
    for mapping,_ in verified.values(): names.update(mapping)
    names.update(f'results/{STAGE}/{name}' for name in SCIENCE)
    diagnostics = sorted(path for path in (OUT/'diagnostics').rglob('*') if path.is_file())
    files = [ROOT/name for name in sorted(names)]+diagnostics
    assert len(files) == len(set(files)) and all(path.is_file() for path in files)
    result = dict(status='Local research milestone; long-term de Sitter research remains open',
        created_at_utc=datetime.now(timezone.utc).isoformat(),
        scientific_status='Three-species H/W1/J2 mixed-q moments and actual low-generator transfer '
            'support normalized J2-W1 real-time OTOC and commutator diagnostics.',
        proof_scope='Fixed low generators J2/W1 and their actual A0/B1 transfer. This stage does not '
            'establish the complete finite-relative-time bilocal four-point function, a vacuum CCR '
            'representation, a bulk mass operator or finite-G scattering. Fixed-word checks and '
            'recorded finite identities accompany separate arbitrary-word/high-moment transfer proofs.',
        parameters='Standard CAR, common complex Gaussian p-body couplings, normalized annealed '
            'beta=mu=0 full-Fock trace, hbar=1/string time; fixed positive lambda=p^2/N after DS. '
            'Numerical OTOC uses Js=1, lambda=.25/.1/.05 and fixed |T|<=1.',
        previous_evidence=previous,checks=checks,
        input_hashes_verified={name:value[0] for name,value in verified.items()},
        parameter_hashes_verified={name:value[1] for name,value in verified.items()},
        evidence_scopes={name:record['scope'] for name,record in science.items()},
        recorded_driver_environments={name:record.get('environment',dict(python=record.get('python'),arithmetic=record.get('arithmetic')))
                                      for name,record in science.items()},
        git_commit_created=False,git_push_performed=False,archive_created=False,external_reply_sent=False,
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_bilocal_four_point.ps1',
            pages=count,passes=2,pass_records=pass_records,compiler_header=log.splitlines()[0],
            render_command=f'pdftoppm -r 95 -png output/pdf/de-sitter-{STAGE}.pdf tmp/pdfs/{STAGE}-render/page',
            final_render_filenames=[path.name for path in pages],render_metadata=metadata,
            final_pages_visually_inspected=list(range(1,count+1)),
            visual_review_attestation='Caller supplied --visual-review-complete after final-page inspection. '
                'Recorder checks hashes, page counts, 95dpi PNG metadata and freshness; it does not see the pages.'),
        diagnostics=dict(preserved_files=[path.relative_to(ROOT).as_posix() for path in diagnostics]),
        environment=dict(python=platform.python_version(),platform=platform.platform(),dependencies='Python standard library only'),
        process_state='Saved science/build/render artifacts only; recorder executes no scientific main, compiler, renderer or Git.',
        hashes={path.relative_to(ROOT).as_posix():sha(path) for path in files},
        render_hashes={path.relative_to(ROOT).as_posix():sha(path) for path in pages},
        build_pass_hashes={path.relative_to(ROOT).as_posix():sha(path) for path in passes})
    with target.open('x',encoding='utf-8') as stream:
        json.dump(result,stream,ensure_ascii=False,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps(dict(status='passed',pages=count,scientific_outputs=4,current_files=len(files),
        final_renders=len(pages),final_build_passes=2,
        previous_scientific_files_verified=sum(len(r['unchanged_scientific_files']) for r in previous.values()))))


if __name__ == '__main__':
    main()
