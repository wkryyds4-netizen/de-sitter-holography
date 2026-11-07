"""Whole-field closed geometry versus the local singleton-normalized source.

This driver imports algorithms but never invokes any predecessor main.
G2 calibration distinguishes energy-only subtraction from charge subtraction.
"""
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import hashlib
import json
import platform
import time

import numpy as np

from memory_body_source_chords import exponent, G_exponent
from memory_oriented_marked_chords import exponent as pure_exponent

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT/'results/susskind-memory-oriented-transfer/body-geometry-checks.json'
PARAMETERS = {
    'layouts': [[2,2,2,0], [3,3,2,0], [4,2,2,0], [2,4,1,1], [2,2,2,2],
                [3,1,2,2], [4,4,0,0], [1,1,3,3], [2,2,0,0], [3,3,1,1]],
    'kinds': ['D', 'L'],
    'maximum_slots': 8,
    'restriction': 'D whole blocks allow internal pairs; L blocks forbid internal pairs',
    'G2_formal_rational_points': [['1/4','1/3'], ['1/2','1/2'], ['1','2/3']],
    'rational_point_scope': 'independent formal lambda,q for exact algebraic identities',
}
PARENTS = ['results/susskind-memory-oriented-transfer/body-source-checks.json',
           'results/susskind-memory-oriented-transfer/oriented-mark-checks.json',
           'results/susskind-memory-oriented-transfer/charge-background-checks.json']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def pairings(slots):
    if not slots:
        yield ()
    else:
        for i in range(1, len(slots)):
            for tail in pairings(slots[1:i]+slots[i+1:]):
                yield ((slots[0], slots[i]),)+tail


def cross(a, b):
    x, y = sorted(a)
    z, w = sorted(b)
    return (x < z < y < w) or (z < x < w < y)


def crossings(chords):
    return sum(cross(a, b) for a, b in combinations(chords, 2))


def block_geometry(pa, start, end):
    internal, opened = [], []
    for chord in pa:
        inside = sum(start <= endpoint < end for endpoint in chord)
        if inside == 2:
            internal.append(chord)
        elif inside == 1:
            opened.append(chord)
    assert len(opened) == end-start-2*len(internal)
    return tuple(internal), frozenset(opened)


def closed_geometry(length, pa, blocks):
    info = [block_geometry(pa, a, b) for a, b, _ in blocks]
    expected = {0: F(-4*crossings(pa))}
    for index, ((inside, opened), (_, _, kind)) in enumerate(zip(info, blocks)):
        d, i = len(opened), len(inside)
        assert kind == 'D' or i == 0
        expected[1 << index] = F(d*(d-1)-4*crossings(opened)-2*i if kind == 'D'
                                  else d*d-4*crossings(opened))
    if len(blocks) == 2:
        common = info[0][1] & info[1][1]
        b, cr = len(common), crossings(common)
        kinds = tuple(kind for _, _, kind in blocks)
        if kinds == ('D', 'D'):
            connected = 4*len(info[0][0])*len(info[1][0])+2*b*(b-1)-4*cr
        elif kinds == ('L', 'L'):
            connected = b*b-4*cr
        else:
            connected = b*(b-1)-4*cr
        expected[3] = F(connected)
    return expected, info


def serial(value):
    return {str(mask): str(v) for mask, v in sorted(value.items())}


def geometry_checks():
    rows, counts = [], Counter()
    for layout in PARAMETERS['layouts']:
        j, k, a, b = layout
        length = sum(layout)
        assert length <= PARAMETERS['maximum_slots']
        for pa in pairings(tuple(range(length))):
            for kinds in product(PARAMETERS['kinds'], repeat=2):
                blocks = ((0,j,kinds[0]), (j+a,j+a+k,kinds[1]))
                info = [block_geometry(pa, start, end) for start, end, _ in blocks]
                if any(kind == 'L' and internal for (internal, _), kind in zip(info, kinds)):
                    continue
                local = exponent(length, pa, blocks)
                expected, info = closed_geometry(length, pa, blocks)
                local_G = G_exponent(length, pa, blocks)
                expected_G = dict(expected)
                for index, ((_, opened), kind) in enumerate(zip(info, kinds)):
                    if kind == 'D':
                        d = len(opened)
                        expected_G[1 << index] -= d*(d-1)
                for mask in range(4):
                    assert local.get(mask, 0) == expected.get(mask, 0), (layout, pa, kinds, mask)
                    assert local_G.get(mask, 0) == expected_G.get(mask, 0)
                    counts['coefficient_comparisons'] += 2
                counts['geometry_cases'] += 1
                counts['cases_'+''.join(kinds)] += 1
                counts['cases_with_internal_pairs'] += any(internal for internal, _ in info)
                counts['mixed_D_L_with_D_internal'] += kinds[0] != kinds[1] and any(internal for internal, _ in info)
                counts['negative_linear_coefficients'] += sum(expected[bit] < 0 for bit in (1,2))
                counts['negative_connected_coefficients'] += expected[3] < 0
                counts['negative_G_linear_coefficients'] += sum(expected_G[bit] < 0 for bit in (1,2))
                common = info[0][1] & info[1][1]
                counts['no_bridge_positive_charge_connected'] += not common and expected[3] > 0
                pure_match = not any(internal for internal, _ in info)
                if pure_match:
                    pure = pure_exponent(length, pa, tuple((x,y) for x,y,_ in blocks), kinds)
                    assert all(pure.get(mask, 0) == local.get(mask, 0) for mask in range(4))
                    counts['pure_Wick_comparisons'] += 4
                    counts['pure_Wick_cases'] += 1
                rows.append({'layout': layout, 'pairing': pa, 'kinds': kinds,
                             'internal_pairs': [internal for internal, _ in info],
                             'open_chords': [sorted(opened) for _, opened in info],
                             'shared_open_chords': sorted(common),
                             'local_kappa': serial({mask: local.get(mask, F(0)) for mask in range(4)}),
                             'closed_kappa': serial(expected), 'G_kappa': serial(expected_G),
                             'closed_geometry_exact': True, 'G_shift_exact': True,
                             'pure_Wick_match_checked': pure_match})
    assert counts['cases_DD'] == 858
    assert counts['mixed_D_L_with_D_internal'] > 0
    assert counts['negative_connected_coefficients'] > 0
    assert counts['no_bridge_positive_charge_connected'] > 0
    return rows, counts


def p_add(*polys):
    out = Counter()
    for poly in polys:
        out.update(poly)
    return {key: value for key, value in out.items() if value}


def p_scale(poly, scalar):
    return {key: scalar*value for key, value in poly.items() if scalar*value}


def p_mul(a, b):
    out = Counter()
    for (i,j), x in a.items():
        for (k,l), y in b.items():
            out[i+k,j+l] += x*y
    return {key: value for key, value in out.items() if value}


def p_eval(poly, lam, q):
    return sum((value*lam**i*q**j for (i,j), value in poly.items()), F(0))


def p_serial(poly):
    return [{'lambda_power': i, 'q_power': j, 'coefficient': str(value)}
            for (i,j), value in sorted(poly.items())]


def G2_calibration():
    diagrams, norm, raw_energy = [], {}, {}
    blocks = ((0,2,'D'), (2,4,'D'))
    for pa in pairings(tuple(range(4))):
        local = G_exponent(4, pa, blocks)
        cr = crossings(pa)
        term = p_add({(1,cr): local.get(3,0)}, {(2,cr): local.get(1,0)*local.get(2,0)})
        norm = p_add(norm, term)
        one = G_exponent(4, pa, ((0,2,'D'),))
        raw_energy = p_add(raw_energy, {(1,cr): one.get(1,0)})
        diagrams.append({'pairing': pa, 'kappa_G': serial(local),
                         'norm_contribution_over_C_squared': p_serial(term)})
    mean_coefficient = G_exponent(2, ((0,1),), ((0,2,'D'),)).get(1,0)
    mean = {(1,0): mean_coefficient}
    h2_inner = p_add(raw_energy, p_scale(mean, -1))  # HH2=H^2-C
    assert norm == {(1,0): F(8), (2,0): F(4), (2,1): F(16)}
    assert mean == {(1,0): F(-2)} and h2_inner == {(1,1): F(-4)}
    one_plus_q = {(0,0): F(1), (0,1): F(1)}
    energy_removed_numerator = p_add(
        p_mul(norm, one_plus_q), p_scale(p_mul(p_mul(mean, mean), one_plus_q), -1),
        p_scale(p_mul(h2_inner, h2_inner), -1))
    charge_variance = {(1,0): F(4)}
    joint_removed_numerator = p_add(energy_removed_numerator,
                                   p_scale(p_mul(charge_variance, one_plus_q), -1))
    old_w_numerator = {(1,0): F(4), (1,1): F(4), (2,1): F(16)}
    assert joint_removed_numerator == old_w_numerator
    assert energy_removed_numerator == p_add(old_w_numerator, p_mul(charge_variance, one_plus_q))
    evaluations = []
    for lam_string, q_string in PARAMETERS['G2_formal_rational_points']:
        lam, q = F(lam_string), F(q_string)
        raw = p_eval(norm, lam, q)
        energy = p_eval(energy_removed_numerator, lam, q)/(1+q)
        joint = p_eval(joint_removed_numerator, lam, q)/(1+q)
        target = 4*(lam+4*lam*lam*q/(1+q))
        assert joint == target and energy-joint == 4*lam
        evaluations.append({'lambda': str(lam), 'q': str(q), 'raw_norm_over_C_squared': str(raw),
                            'energy_only_residual_over_C_squared': str(energy),
                            'charge_and_energy_residual_over_C_squared': str(joint),
                            'old_w_over_C_squared': str(target)})
    return {'diagrams': diagrams, 'raw_norm_over_C_squared': '8*lambda+4*lambda^2+16*lambda^2*q',
            'raw_norm_polynomial': p_serial(norm), 'mean_over_C': '-2*lambda',
            'HH2_inner_over_C_squared': '-4*lambda*q',
            'HH2_norm_over_C_squared': '1+q',
            'deterministic_energy_regression': '-2*lambda*C-4*lambda*q/(1+q)*HH2(H)',
            'energy_only_residual_over_C_squared': '8*lambda+16*lambda^2*q/(1+q)',
            'energy_only_residual_numerator_over_1_plus_q': p_serial(energy_removed_numerator),
            'centered_charge_static_variance_over_C_squared': '4*lambda',
            'charge_and_energy_residual_over_C_squared': '4*lambda+16*lambda^2*q/(1+q)',
            'joint_residual_numerator_over_1_plus_q': p_serial(joint_removed_numerator),
            'old_w_recovered_after_charge_subtraction': True,
            'energy_only_equals_old_w': False, 'formal_rational_evaluations': evaluations,
            'scope': 'Exact low-degree diagram algebra and subtraction using the separately identified charge-static variance; not a new proof of general joint regression or absence of other static mass.'}


def main():
    if OUT.exists():
        raise FileExistsError(OUT)
    started = time.perf_counter()
    parents = {}
    for name in PARENTS:
        data = json.loads((ROOT/name).read_text(encoding='utf-8'))
        assert data['status'] == 'passed'
        for input_name, expected in data['input_sha256'].items():
            assert sha(ROOT/input_name) == expected, (name, input_name)
        parents[name] = {'sha256': sha(ROOT/name), 'input_hashes_verified': True}
    rows, counts = geometry_checks()
    calibration = G2_calibration()
    counts['G2_diagrams'] = len(calibration['diagrams'])
    counts['G2_formal_rational_checks'] = len(calibration['formal_rational_evaluations'])
    names = ['code/check_memory_body_geometry.py', 'code/memory_body_source_chords.py',
             'code/memory_oriented_marked_chords.py', 'code/memory_size_marked_chords.py',
             'code/singlet_coefficient.py', 'code/singlet_multimark.py',
             'code/complex_syk_dictionary.py']+PARENTS
    result = {'status': 'passed', 'schema_version': 1, 'created_utc': datetime.now(timezone.utc).isoformat(),
              'scope': 'Exact finite-diagram closed geometry versus singleton-normalized local algorithm, including internal D pairs and mixed D/L. L remains coupling-Wick restricted. Signed terms are retained; this finite set does not replace the all-degree proof.',
              'parameters': PARAMETERS, 'parameter_sha256': digest(PARAMETERS),
              'environment': {'python': platform.python_version(), 'numpy': np.__version__, 'platform': platform.platform()},
              'input_sha256': {name: sha(ROOT/name) for name in names}, 'parents_verified': parents,
              'counts': dict(counts), 'geometry_rows': rows, 'G2_calibration': calibration,
              'elapsed_seconds': time.perf_counter()-started}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2)
        stream.write('\n')
    print(json.dumps({'status': result['status'], 'counts': result['counts'],
                      'parameter_sha256': result['parameter_sha256'], 'output_sha256': sha(OUT),
                      'elapsed_seconds': result['elapsed_seconds']}, sort_keys=True))


if __name__ == '__main__':
    main()
