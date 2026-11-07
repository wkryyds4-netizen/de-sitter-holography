"""Exact complete-CAR entry checks for the oriented transfer defect.

Uses Gaussian-integer matrices and generalized Hermite polynomials with
coordinate variances 2 (diagonal) and 1 (off diagonal). No sample, rank cut,
scientific project import, historical main or spectrum computation is used.
"""
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, factorial, prod
from pathlib import Path
import hashlib
import json
import platform
import time

import numpy as np


ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT/'results/susskind-memory-oriented-transfer/finite-transfer-checks.json'
PARAMETERS = dict(
    systems=[[2, 1, 4], [3, 1, 3], [3, 2, 3], [4, 2, 2]],
    coordinate_variances=dict(diagonal=2, real_off_diagonal=1, imaginary_off_diagonal=1),
    finite_energy_normalization='C_m=2 binom(m,p) binom(N-m+p,p); no physical Js rescaling',
    alpha='max(1-p/m,0) for m>0; m=0 defect defined as zero',
    hermite='Monic multivariate Hermite with recorded coordinate variances; number=sum(multiindex)',
    zero_power_convention='0^0=1, including alpha=0 at m=p',
    taylor_small_m=list(range(1, 25)), taylor_small_k=list(range(41)),
    taylor_large_p=[1, 2, 3, 5, 13, 34, 89],
    taylor_large_m='unique sorted {p,p+1,2p,2p+1,7p+3,p^2+1}',
    taylor_large_k=[0, 1, 2, 3, 4, 7, 16, 41, 64, 127, 256],
    arithmetic='Exact numpy int64 real/imaginary arrays with conservative pre-operation overflow guards; Fraction for scalar Taylor tests',
)
LIMIT = 2**60
MAXIMUM_INTEGER = 0
COUNTS = Counter()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def bound(matrix):
    global MAXIMUM_INTEGER
    value = int(np.max(np.abs(matrix), initial=0))
    assert value < LIMIT
    MAXIMUM_INTEGER = max(MAXIMUM_INTEGER, value)
    return value


def gi(real, imag=None):
    real = np.asarray(real, dtype=np.int64)
    return np.stack([real, np.zeros_like(real) if imag is None else np.asarray(imag, dtype=np.int64)])


def add(left, right, factor=1):
    assert bound(left)+abs(factor)*bound(right) < LIMIT
    answer = left+factor*right
    bound(answer)
    return answer


def scale(matrix, factor):
    factor = int(factor)
    assert abs(factor)*bound(matrix) < LIMIT
    answer = factor*matrix
    bound(answer)
    return answer


def mul(left, right):
    assert 2*left.shape[-1]*bound(left)*bound(right) < LIMIT
    answer = np.stack([left[0]@right[0]-left[1]@right[1],
                       left[0]@right[1]+left[1]@right[0]])
    bound(answer)
    return answer


def dagger(matrix):
    return np.stack([matrix[0].T, -matrix[1].T])


def restrict(matrix, ids):
    return matrix[:, ids][:, :, ids]


def same(left, right):
    assert np.array_equal(left, right)
    COUNTS['matrix_identities'] += 1
    COUNTS['integer_component_equalities'] += int(left.size)


def zero(matrix):
    same(matrix, np.zeros_like(matrix))


def car_system(N, p):
    dimension = 1 << N
    eye = gi(np.eye(dimension, dtype=np.int64))
    cs = []
    for index in range(N):
        c = np.zeros((dimension, dimension), dtype=np.int64)
        for state in range(dimension):
            if state & (1 << index):
                c[state ^ (1 << index), state] = (-1)**((state & ((1 << index)-1)).bit_count())
        cs.append(gi(c))
    for i in range(N):
        for j in range(N):
            same(add(mul(cs[i], dagger(cs[j])), mul(dagger(cs[j]), cs[i])),
                 eye if i == j else np.zeros_like(eye))
    charges = np.array([state.bit_count() for state in range(dimension)], dtype=np.int64)
    charge = gi(np.diag(charges))
    creators = []
    for indices in combinations(range(N), p):
        operator = eye
        for index in indices:
            operator = mul(operator, dagger(cs[index]))
        creators.append(operator)
    basis, variances = [], []
    for i, first in enumerate(creators):
        basis.append(mul(first, dagger(first)))
        variances.append(2)
        for second in creators[:i]:
            term = mul(first, dagger(second))
            basis.append(add(term, dagger(term)))
            antisymmetric = add(term, dagger(term), -1)
            basis.append(np.stack([-antisymmetric[1], antisymmetric[0]]))
            variances.extend([1, 1])
    assert len(basis) == comb(N, p)**2

    def D(matrix):
        value = mul(charge, matrix)
        for c in cs:
            value = add(value, mul(mul(dagger(c), matrix), c), -1)
        return value

    zero(D(eye))
    for matrix in basis:
        same(D(matrix), scale(matrix, p))
        same(matrix, dagger(matrix))
        same(mul(charge, matrix), mul(matrix, charge))
    variance = np.zeros_like(eye)
    for matrix, value in zip(basis, variances):
        variance = add(variance, mul(matrix, matrix), value)
    sectors = [np.flatnonzero(charges == m) for m in range(N+1)]
    for m, ids in enumerate(sectors):
        C = 2*comb(m, p)*comb(N-m+p, p)
        same(restrict(variance, ids), scale(restrict(eye, ids), C))
    return eye, cs, basis, variances, sectors, D


def poly_add(out, key, matrix, factor=1):
    out[key] = add(out[key], matrix, factor) if key in out else scale(matrix, factor)


def step(poly, basis, variances, hermite):
    out = {}
    for key, matrix in poly.items():
        for i, B in enumerate(basis):
            value = mul(B, matrix)
            raised = list(key)
            raised[i] += 1
            poly_add(out, tuple(raised), value)
            if hermite and key[i]:
                lowered = list(key)
                lowered[i] -= 1
                poly_add(out, tuple(lowered), value, key[i]*variances[i])
    return out


def reductions(key, variances, hermite_to_monomial, alpha_numerator=1, noise_numerator=0):
    """Independent scalar product formulas, without matrix or chaos recursion."""
    for removed in product(*(range(entry//2+1) for entry in key)):
        output, factor = [], 1
        for degree, pairs, variance in zip(key, removed, variances):
            unpaired = degree-2*pairs
            coefficient = factorial(degree)//(2**pairs*factorial(pairs)*factorial(unpaired))
            output.append(unpaired)
            if hermite_to_monomial:
                factor *= coefficient*(-variance)**pairs
            else:
                factor *= coefficient*variance**pairs*noise_numerator**pairs*alpha_numerator**unpaired
        yield tuple(output), factor


def polynomial_same(left, right, shape):
    blank = np.zeros(shape, dtype=np.int64)
    for key in left.keys() | right.keys():
        same(left.get(key, blank), right.get(key, blank))


def finite_check(N, p, maximum_degree):
    eye, cs, basis, variances, sectors, D = car_system(N, p)
    before = COUNTS.copy()
    unit_checks = 0
    # Every matrix unit in the complete neutral algebra, supported in one charge.
    for source_m, source_ids in enumerate(sectors):
        for row in source_ids:
            for column in source_ids:
                matrix = np.zeros_like(eye)
                matrix[0, row, column] = 1
                defect = D(matrix)
                for m in range(1, N+1):
                    ids, low = sectors[m], sectors[m-1]
                    removal = np.zeros((2, len(ids), len(ids)), dtype=np.int64)
                    for c in cs:
                        rectangular = c[:, low][:, :, ids]
                        removal = add(removal, mul(mul(dagger(rectangular), restrict(matrix, low)), rectangular))
                    same(restrict(defect, ids), add(scale(restrict(matrix, ids), m), removal, -1))
                    unit_checks += 1
                zero(restrict(defect, sectors[0]))
    key0 = (0,)*len(basis)
    hermite, raw = {key0: eye}, {key0: eye}
    rows = []
    for degree in range(maximum_degree+1):
        defect = {key: D(value) for key, value in hermite.items()}
        for m, ids in enumerate(sectors):
            if m == 0:
                for key, matrix in hermite.items():
                    zero(restrict(defect[key], ids))
                    if degree:
                        zero(restrict(matrix, ids))
                rows.append(dict(degree=degree, m=m, alpha=None, Hermite_coefficients=len(hermite),
                                 zero_charge_convention=True, E0_zero=True))
                continue
            numerator = max(m-p, 0)
            low = sectors[m-1]
            weighted_hermite, conditional_raw = {}, {}
            nonzero_E = 0
            for key, matrix in hermite.items():
                chaos = sum(key)
                upper, lower = restrict(matrix, ids), restrict(matrix, low)
                removal = np.zeros_like(upper)
                for c in cs:
                    rectangular = c[:, low][:, :, ids]
                    removal = add(removal, mul(mul(dagger(rectangular), lower), rectangular))
                same(restrict(defect[key], ids), add(scale(upper, m), removal, -1))
                # Clear m^chaos; no division and 0^0 remains 1 at m=p.
                left = add(scale(removal, m**chaos), scale(upper, m*numerator**chaos), -1)
                right = add(scale(upper, m**(chaos+1)-m*numerator**chaos),
                            scale(restrict(defect[key], ids), m**chaos), -1)
                same(left, right)
                if degree <= 1 or m < p:
                    zero(left)
                nonzero_E += int(np.any(left))
                COUNTS['Hermite_transfer_coefficient_cases'] += 1
                for output_key, coefficient in reductions(key, variances, True):
                    poly_add(weighted_hermite, output_key, upper,
                             coefficient*m**(degree-chaos)*numerator**chaos)
            # Direct conditional Gaussian expectation of the raw polynomial,
            # using an independent coupling J' and integrating all its coordinates.
            for key, matrix in raw.items():
                assert sum(key) == degree
                for output_key, coefficient in reductions(key, variances, False,
                                                           numerator, m*m-numerator*numerator):
                    poly_add(conditional_raw, output_key, restrict(matrix, ids), coefficient)
            polynomial_same(weighted_hermite, conditional_raw, (2, len(ids), len(ids)))
            COUNTS['complete_OU_polynomial_cases'] += 1
            rows.append(dict(degree=degree, m=m, alpha=str(F(numerator, m)),
                             Hermite_coefficients=len(hermite), raw_coefficients=len(raw),
                             OU_output_monomials=len(conditional_raw), nonzero_E_coefficients=nonzero_E,
                             transfer_and_number_identity=True, conditional_OU_identity=True,
                             initial_E_zero=(degree <= 1), low_charge_E_zero=(m < p)))
        if degree < maximum_degree:
            hermite = step(hermite, basis, variances, True)
            raw = step(raw, basis, variances, False)
    assert sum(row.get('nonzero_E_coefficients', 0) for row in rows) > 0
    return dict(N=N, p=p, maximum_degree=maximum_degree, full_Fock_dimension=len(eye[0]),
                real_Gaussian_coordinates=len(basis), all_sectors=list(range(N+1)),
                neutral_matrix_units=sum(len(ids)**2 for ids in sectors),
                matrix_unit_sector_checks=unit_checks,
                full_basis_DH_equals_pH=len(basis), variance_check_all_sectors=True,
                rows=rows, counts=dict(COUNTS-before))


def taylor_checks():
    cases = {(m, p, k) for m in PARAMETERS['taylor_small_m'] for p in range(1, m+1)
             for k in PARAMETERS['taylor_small_k']}
    for p in PARAMETERS['taylor_large_p']:
        for m in {p, p+1, 2*p, 2*p+1, 7*p+3, p*p+1}:
            cases.update((m, p, k) for k in PARAMETERS['taylor_large_k'])
    rows = []
    for m, p, k in sorted(cases):
        remainder = m*(1-F(m-p, m)**k)-k*p+F(comb(k, 2)*p*p, m)
        upper = F(comb(k, 3)*p**3, m*m)
        assert 0 <= remainder <= upper, (m, p, k, remainder, upper)
        if k <= 2:
            assert remainder == upper == 0
        if k == 3:
            assert remainder == upper
        rows.append(dict(m=m, p=p, k=k, remainder=str(remainder), upper_bound=str(upper)))
    return dict(rows=rows, counts=dict(cases=len(rows),
                endpoint_m_equals_p=sum(row['m'] == row['p'] for row in rows),
                k_zero_to_three=sum(row['k'] <= 3 for row in rows),
                zero_remainder=sum(row['remainder'] == '0' for row in rows),
                upper_equality=sum(row['remainder'] == row['upper_bound'] for row in rows)),
                rows_sha256=canonical_hash(rows),
                scope='Exact scalar finite tests of the all-k Taylor inequality; the general proof is separate.')


def main():
    if OUT.exists():
        raise FileExistsError(OUT)
    started = time.perf_counter()
    systems = []
    for args in PARAMETERS['systems']:
        result = finite_check(*args)
        systems.append(result)
        print(json.dumps(dict(system=args, status='passed', counts=result['counts'])), flush=True)
    taylor = taylor_checks()
    data = dict(schema_version=1, status='passed',
                created_utc=datetime.now(timezone.utc).isoformat(),
                parameters=PARAMETERS, parameter_sha256=canonical_hash(PARAMETERS),
                systems=systems, counts=dict(COUNTS), taylor=taylor,
                maximum_absolute_integer=MAXIMUM_INTEGER, pre_operation_integer_limit=LIMIT,
                environment=dict(python=platform.python_version(), numpy=np.__version__,
                                 platform=platform.platform()),
                input_sha256={'code/check_memory_oriented_transfer.py': sha(Path(__file__))},
                elapsed_seconds=time.perf_counter()-started,
                scope='Complete finite neutral matrix units and every Gaussian-Hermite coefficient '
                      'of the recorded ordinary powers, exact common-J transfer/D/OU identities, '
                      'DH=pH, E(0)=E-prime(0)=0, m=0 convention, m<p vanishing and m=p boundary. '
                      'Direct conditional Gaussian integration independently checks the number filter. '
                      'Finite tests do not prove a double-scaling rate, a full-time weighted estimate '
                      'or an AA/AB/BA spectral kernel; physical Js units require separate rescaling.')
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(data, indent=2)+'\n')
    print(json.dumps(dict(status='passed', counts=data['counts'], taylor=taylor['counts'],
                         maximum_absolute_integer=MAXIMUM_INTEGER,
                         elapsed_seconds=data['elapsed_seconds'], output_sha256=sha(OUT))), flush=True)


if __name__ == '__main__':
    main()
