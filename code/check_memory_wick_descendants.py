"""Fresh OU/Wick-boundary checks and signed descendant chord moments.

Finite comparisons use full-Fock half-Clifford-size before restriction, never
identify finite rank with size, and retain every real Gaussian coordinate.
The formal C=1 diagram calculations retain internal contractions, scalar chaos,
reversed adjoints, and the distinction between Wick and mark boundaries.
"""
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb, factorial
from pathlib import Path
import hashlib
import json
import platform
import sys
import time

import mpmath as mp
import numpy as np
import sympy as sp

from memory_charged_transfer import system
from memory_jacobi_bridge import step, raw_product
from memory_wick_rank import mu
from singlet_wick_moments import gue_basis
from singlet_coefficient import pairings
from memory_size_marked_chords import moment as prior_U_moment
import memory_wick_descendant_chords as dc


ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT/'results/susskind-memory-wick-descendants/descendant-checks.json'
ATOL, RTOL = 2e-10, 2e-10
LAMBDA, Q = sp.symbols('lambda q')
FINITE_LAYOUTS = [
    ((0, 1, 1), (0, 1, 1), 0, 0),
    ((1, 1, 0), (1, 1, 0), 0, 0),
    ((2, 1, 0), (0, 1, 0), 0, 0),
    ((1, 1, 1), (1, 1, 1), 0, 0),
    ((0, 2, 1), (0, 2, 1), 0, 0),
    ((0, 1, 2), (0, 2, 1), 0, 0),
    ((0, 1, 1), (0, 1, 1), 1, 1),
]


def complex_record(value):
    return dict(real=float(value.real), imag=float(value.imag))


def primitive_term(spec):
    prefix, v, w = spec
    return dc.Term(prefix+v+w, ((prefix, prefix+v+w),), ((prefix, prefix+v),))


def right_step(poly, basis):
    dual = {key: value.conj().T for key, value in poly.items()}
    return {key: value.conj().T for key, value in step(dual, basis).items()}


def finite_checks(N, p, m):
    cs, integer_basis, variances = gue_basis(N, p)
    Cbase = comb(m, p)*comb(N-m+p, p)
    basis = [matrix*np.sqrt(float(variance)/2/Cbase)
             for matrix, variance in zip(integer_basis, variances)]
    ids = [state for state in range(1 << N) if state.bit_count() == m]
    dimension = len(ids)
    comparison = system(N, p, m)
    assert len(basis) == len(comparison['up'])
    assert max(np.linalg.norm(matrix[np.ix_(ids, ids)]-other)
               for matrix, other in zip(basis, comparison['up'])) < 1e-13
    count = len(basis)
    zero = (0,)*count
    eye = np.eye(1 << N, dtype=complex)
    halfsize = lambda X: N*X/2-sum(c.conj().T@X@c+c@X@c.conj().T for c in cs)/2
    wick = [{zero: eye}]
    for _ in range(3):
        wick.append(step(wick[-1], basis, 'creation'))

    @lru_cache(None)
    def hermite(spec):
        prefix, v, w = spec
        marked = {key: v*p*value-halfsize(value) for key, value in wick[v].items()}
        result = raw_product(marked, wick[w])
        for _ in range(prefix):
            result = step(result, basis)
        return result

    def inner(left, right):
        total = 0j
        for key in left.keys() & right.keys():
            weight = 1
            for entry in key:
                weight *= factorial(entry)
            first = left[key][np.ix_(ids, ids)]
            second = right[key][np.ix_(ids, ids)]
            total += weight*np.vdot(first, second)/dimension
        return total

    filters = [
        ('identity', lambda degree: 1.0, 'identity', lambda degree: 1.0),
        ('OU_1/2', lambda degree: 0.5**degree, 'OU_3/4', lambda degree: 0.75**degree),
        ('finite_mu', lambda degree: float(mu(N, p, m, degree*p)),
         'OU_1/2', lambda degree: 0.5**degree),
    ]
    records = []
    for left_spec, right_spec, a, b in FINITE_LAYOUTS:
        left, right = primitive_term(left_spec), primitive_term(right_spec)
        _, length, wicks, marks, intervals = dc.layout(left, right, a, b)
        assert length <= 6
        starts = {start: end for start, end in marks}

        @lru_cache(None)
        def word_trace(labels):
            result = eye.copy()
            position = 0
            while position < length:
                if position in starts:
                    end = starts[position]
                    block = eye.copy()
                    for slot in range(position, end):
                        block = block@basis[labels[slot]]
                    block = (end-position)*p*block-halfsize(block)
                else:
                    end, block = position+1, basis[labels[position]]
                result = result@block
                position = end
            return np.trace(result[np.ix_(ids, ids)])/dimension

        grouped = Counter()
        allowed_pairings = coordinate_assignments = 0
        for pairing in pairings(tuple(range(length))):
            if not dc.allowed(pairing, wicks):
                continue
            degrees = tuple(dc.outgoing_degree(pairing, interval) for interval in intervals)
            total = 0j
            for colors in product(range(count), repeat=len(pairing)):
                labels = [None]*length
                for (x, y), color in zip(pairing, colors):
                    labels[x] = labels[y] = color
                total += word_trace(tuple(labels))
                coordinate_assignments += 1
            grouped[degrees] += total
            allowed_pairings += 1
        for left_name, left_weight, right_name, right_weight in filters:
            Fpoly = {key: left_weight(sum(key))*value for key, value in hermite(left_spec).items()}
            Gpoly = {key: right_weight(sum(key))*value for key, value in hermite(right_spec).items()}
            for _ in range(a):
                Gpoly = step(Gpoly, basis)
            for _ in range(b):
                Gpoly = right_step(Gpoly, basis)
            expected = inner(Fpoly, Gpoly)
            observed = sum(left_weight(dl)*right_weight(dr)*value
                           for (dl, dr), value in grouped.items())
            error = abs(observed-expected)
            tolerance = ATOL+RTOL*max(abs(observed), abs(expected))
            assert error <= tolerance, (N, p, m, left_spec, right_spec, left_name, error)
            records.append(dict(N=N, p=p, m=m, C_sector=1,
                                full_Fock_dimension=1 << N, sector_dimension=dimension,
                                Gaussian_coordinates=count, left_primitive=left_spec,
                                right_primitive=right_spec, left_energy=a, right_energy=b,
                                raw_slots=length, left_filter=left_name, right_filter=right_name,
                                allowed_pairings=allowed_pairings,
                                coordinate_assignments=coordinate_assignments,
                                distinct_CAR_words=word_trace.cache_info().currsize,
                                Hermite_value=complex_record(expected),
                                internal_pair_count_value=complex_record(observed),
                                absolute_error=float(error), tolerance=float(tolerance),
                                grouped_outgoing_chaoses=[dict(left_degree=dl, right_degree=dr,
                                                               value=complex_record(value))
                                                         for (dl, dr), value in sorted(grouped.items())]))
    print(f'finite full-Fock OU/Hermite N={N},p={p},m={m}: {len(records)} comparisons', flush=True)
    return records


def symbolic(poly):
    return sp.expand(sum(sp.Rational(value.numerator, value.denominator)*LAMBDA**a*Q**b
                         for (a, b), value in poly.items()))


def P(j):
    return sp.prod(sum(Q**r for r in range(k)) for k in range(1, j+1))


def check_low_orders(eta, V):
    expected3 = dc.fscale(dc.U(1), dc.mono(1, 1, -4), degree=3)
    assert eta[3].terms == expected3.terms
    anticommutator = dc.fadd(dc.hleft(dc.U(1)), dc.hright(dc.U(1)))
    expected4 = dc.fscale(anticommutator, dc.mono(1, 1, -4), degree=4)
    expected4 = dc.fadd(expected4, dc.fscale(dc.U(2), dc.mono(1, 2, -4), degree=4))
    expected4 = dc.fadd(expected4, dc.fscale(dc.identity(4), dc.mono(2, 1, 4)))
    difference = dc.fadd(eta[4], expected4, -1)
    assert dc.moment(difference, difference, max_slots=8)[0] == {}
    mean4 = dc.moment(dc.identity(), eta[4], max_slots=8)[0]
    assert mean4 == dc.mono(2, 1, -4)
    assert dc.moment(dc.identity(), dc.ou(eta[4]), max_slots=8)[0] == mean4
    assert dc.moment(dc.select(eta[4], 0), dc.select(eta[4], 0), max_slots=8)[0] == dc.mono(4, 2, 16)
    expected_norms = {
        1: {(3, 2): F(16), (4, 2): F(16)},
        2: {(3, 4): F(64), (4, 4): F(256)},
    }
    for j, expected in expected_norms.items():
        assert V[j].terms == dc.fscale(dc.U(j), dc.mono(1, j, 4)).terms
        assert dc.moment(V[j], V[j], max_slots=8)[0] == expected
    calibrations = []
    for j in [1, 2]:
        for k in [1, 2]:
            for a, b in [(0, 0), (1, 1), (2, 0)]:
                prior, _ = prior_U_moment(j, k, a, b, 3)
                expected = dc.pmul(prior, dc.mono(2, j+k, 16))
                actual, _ = dc.moment(V[j], V[k], a, b, max_slots=8)
                assert actual == expected
                calibrations.append(dict(j=j, k=k, left_energy=a, right_energy=b,
                                         polynomial=dc.serialize_polynomial(actual)))
    return dict(eta3_identity='eta3=-4 lambda q U1 (C=1)',
                eta4_identity='eta4=-4 lambda [q {H,U1}+q^2 U2-lambda q I] (C=1)',
                eta4_scalar=dc.serialize_polynomial(mean4),
                eta4_scalar_squared_norm=dc.serialize_polynomial(dc.mono(4, 2, 16)),
                eta4_explicit_difference_norm_squared=[],
                V1_V2_norms={str(j): dc.serialize_polynomial(poly) for j, poly in expected_norms.items()},
                prior_V1_V2_moment_calibrations=calibrations)


def regression(field):
    result = {}
    for degree in range(field.degree+1):
        numerator = symbolic(dc.moment(dc.Z(degree), field, max_slots=8)[0])
        if numerator != 0:
            result[degree] = sp.factor(numerator/P(degree))
    return result


def no_bridge_checks(eta):
    fields = {'eta3': eta[3], 'eta4': eta[4],
              'OU_eta4': dc.ou(eta[4]), 'OU_eta5': dc.ou(eta[5])}
    regressions = {name: regression(field) for name, field in fields.items()}
    assert regressions['eta3'] == {1: -4*LAMBDA**2*Q}
    assert regressions['eta4'][0] == -4*LAMBDA**2*Q
    for degree, value in regressions['eta4'].items():
        assert sp.cancel(regressions['OU_eta4'][degree]-Q**degree*value) == 0
    cases = [('eta3', 'eta3', 1, 1), ('eta3', 'eta4', 1, 2),
             ('eta4', 'eta4', 2, 2), ('OU_eta4', 'eta4', 2, 2),
             ('OU_eta5', 'eta3', 2, 2), ('OU_eta5', 'OU_eta5', 1, 1)]
    rows = []
    for first, second, a, b in cases:
        left, right = fields[first], fields[second]
        calculate = lambda Ffield, Gfield: symbolic(dc.moment(
            Ffield, Gfield, a, b, max_slots=8, no_field_bridge=True)[0])
        raw = calculate(left, right)
        left_sub = sum(value*calculate(dc.Z(degree), right)
                       for degree, value in regressions[first].items())
        right_sub = sum(value*calculate(left, dc.Z(degree))
                        for degree, value in regressions[second].items())
        added = sum(value*other*calculate(dc.Z(degree), dc.Z(other_degree))
                    for degree, value in regressions[first].items()
                    for other_degree, other in regressions[second].items())
        residual = sp.cancel(raw-left_sub-right_sub+added)
        assert residual == 0, (first, second, a, b, residual)
        assert raw != 0, (first, second, 'vacuous zero-bridge fixture')
        rows.append(dict(left=first, right=second, left_energy=a, right_energy=b,
                         uncentered_no_bridge=str(sp.factor(raw)),
                         subtract_left_regression=str(sp.factor(left_sub)),
                         subtract_right_regression=str(sp.factor(right_sub)),
                         add_double_regression=str(sp.factor(added)), centered_residual=str(residual)))
    return dict(regressions={name: {str(degree): str(value) for degree, value in values.items()}
                             for name, values in regressions.items()}, rows=rows,
                scope='Exact finite-degree diagram cancellations after the displayed deterministic '
                      'regressions. These finite checks alone do not establish general spectral '
                      'absolute continuity, convergence of exact energy pinching, or a time kernel.')


def source_regression_checks():
    rows = []
    for n in range(1, 5):
        for v in range(1, n+1):
            w = n-v
            actual = symbolic(dc.moment(dc.Z(n), dc.composite(v, w), max_slots=8)[0])
            pv = sp.expand(P(v))
            predicted = sp.cancel(LAMBDA*(v*v*pv-4*Q*sp.diff(pv, Q))*P(n)/pv)
            assert sp.cancel(actual-predicted) == 0
            rows.append(dict(v=v, w=w, total_degree=n, numerator=str(sp.factor(actual)),
                             regression_ratio=str(sp.factor(actual/P(n))),
                             predicted_ratio=str(sp.factor(predicted/P(n)))))
    return dict(rows=rows, scope='Low-degree test of the a_v regression prediction for the '
                                 'composite source, not an all-degree proof by enumeration.')


def formal_moments(eta, V):
    rows = []
    energies = [(0, 0), (1, 0), (1, 1), (2, 0)]
    for family, fields, indices in [('eta', eta, [3, 4, 5]), ('V', V, [1, 2, 3, 4])]:
        for j in indices:
            for k in indices:
                if j > k:
                    continue
                for a, b in energies:
                    largest = max((tl.length+tr.length+a+b for tl in fields[j].terms
                                   for tr in fields[k].terms), default=0)
                    if largest > 8:
                        continue
                    polynomial, info = dc.moment(fields[j], fields[k], a, b, max_slots=8)
                    reverse, _ = dc.moment(fields[k], fields[j], a, b, max_slots=8)
                    # Real coefficients plus cyclic trace give equality under swapping fields
                    # together with the corresponding Hermitian adjoint in these limits.
                    if a == b:
                        assert polynomial == reverse
                    rows.append(dict(family=family, left_degree=j, right_degree=k,
                                     left_energy=a, right_energy=b, **info,
                                     C_restoration_exponent=str(F(j+k+a+b, 2)),
                                     polynomial=dc.serialize_polynomial(polynomial),
                                     negative_coefficients=sum(value < 0 for value in polynomial.values())))
    chaos = []
    for j in [3, 4, 5]:
        total, _ = dc.moment(eta[j], eta[j], max_slots=8)
        assembled = {}
        for degree in range(j-1):
            value, _ = dc.moment(dc.select(eta[j], degree), dc.select(eta[j], degree), max_slots=8)
            assembled = dc.padd(assembled, value)
            chaos.append(dict(eta_degree=j, coupling_chaos=degree,
                              squared_norm=dc.serialize_polynomial(value)))
        assert assembled == total
    return dict(rows=rows, eta_chaos_norms=chaos)


def source_hashes():
    paths = {Path(__file__).resolve()}
    for module in tuple(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename:
            path = Path(filename).resolve()
            if path.parent == ROOT/'code' and path.suffix == '.py':
                paths.add(path)
    return {path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(paths)}


def main():
    if OUT.exists():
        raise FileExistsError(OUT)
    started = time.time()
    finite = []
    for N, p, m in [(2, 1, 1), (3, 1, 2)]:
        finite.extend(finite_checks(N, p, m))
    eta, V, sources = dc.descendants(5)
    low = check_low_orders(eta, V)
    regression_checks = source_regression_checks()
    centered = no_bridge_checks(eta)
    moments = formal_moments(eta, V)
    parameters = dict(finite_systems=[(2, 1, 1), (3, 1, 2)], finite_layouts=FINITE_LAYOUTS,
                      finite_maximum_slots=6, formal_maximum_slots=8, maximum_descendant_degree=5,
                      absolute_tolerance=ATOL, relative_tolerance=RTOL, C=1)
    encoded = json.dumps(parameters, sort_keys=True, separators=(',', ':')).encode()
    counts = dict(finite_OU_Hermite_comparisons=len(finite),
                  prior_V1_V2_calibrations=len(low['prior_V1_V2_moment_calibrations']),
                  composite_source_regressions=len(regression_checks['rows']),
                  nonvacuous_centered_no_bridge_cases=len(centered['rows']),
                  formal_moment_polynomials=len(moments['rows']),
                  eta_chaos_norms=len(moments['eta_chaos_norms']))
    result = dict(status='passed', created_at_utc=datetime.now(timezone.utc).isoformat(),
                  parameters=parameters, parameter_sha256=hashlib.sha256(encoded).hexdigest(),
                  counts=counts, finite_comparisons=finite,
                  maximum_finite_absolute_error=max(row['absolute_error'] for row in finite),
                  fields=dict(eta=[dc.serialize_field(eta[j]) for j in [3, 4, 5]],
                              V=[dc.serialize_field(V[j]) for j in [1, 2, 3, 4]],
                              sources=[dc.serialize_field(sources[j]) for j in sorted(sources)]),
                  low_order_checks=low, composite_source_regression=regression_checks,
                  centered_zero_bridge=centered, moment_checks=moments,
                  cache_statistics=dict(term_moments=dc.term_moment.cache_info()._asdict(),
                                        local_exponent=dc.exponent.cache_info()._asdict()),
                  conventions='Moment order is <F^dagger H^a G H^b>. Dagger reverses marks '
                              'and Wick groups. Composite sources have a whole v+w Wick block '
                              'but a size mark only on its v subblock. The OU degree is raw '
                              'field slots minus twice the number of pairs internal to that '
                              'whole field; external energy slots are excluded. Every field '
                              'retains its physical energy degree, so C restores as '
                              'C^((j+k+a+b)/2). No samplewise pinching is executed.',
                  environment=dict(python=platform.python_version(), numpy=np.__version__,
                                   sympy=sp.__version__, mpmath=mp.__version__),
                  input_sha256=source_hashes(), elapsed_seconds=time.time()-started,
                  scope='Full finite Gaussian-Hermite versus coordinate-pairing comparisons '
                        'use full-Fock half-size and floating normalized CAR matrices. The '
                        'formal signed lambda/q polynomials and centered no-bridge cancellations '
                        'are exact fixed-degree chord computations under the established '
                        'pair-Poisson rule. No finite-N rank=size substitution, arbitrary-degree '
                        'proof by finite enumeration, interval arithmetic, or full-time '
                        'pinched/charged memory kernel is claimed.',
                  reproducibility='Only this new main is executed. Historical sources and '
                                  'outputs are unchanged; output is exclusive-create.')
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps(dict(status='passed', counts=counts,
                          maximum_finite_absolute_error=result['maximum_finite_absolute_error'],
                          elapsed_seconds=result['elapsed_seconds'],
                          output_sha256=hashlib.sha256(OUT.read_bytes()).hexdigest())), flush=True)


if __name__ == '__main__':
    main()
