"""Fresh exact finite-CAR checks for the fixed-chaos rank/Clifford-size bridge.

The full-Fock D/Db compositions are evaluated before charge restriction.
Integer fixtures and every highest-Gaussian-chaos coefficient of one small
model are checked separately. No asymptotic bound is inferred from samples.
"""
from datetime import datetime, timezone
from itertools import combinations
from math import factorial
from pathlib import Path
import hashlib
import json
import platform
import sys
import time

import numpy as np
import sympy as sp

from complex_syk_dictionary import fermions, creators
from memory_rank_blocks import casimir_matrix, projectors
from singlet_wick_moments import gue_basis


ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'results/susskind-memory-size-bridge/finite-bridge.json'
SEED = 11213764
MAX_PRODUCT_BOUND = 0


def checked(matrix):
    matrix = np.asarray(matrix, dtype=np.complex128)
    assert np.array_equal(matrix.real, np.rint(matrix.real))
    assert np.array_equal(matrix.imag, np.rint(matrix.imag))
    assert max(float(np.max(abs(matrix.real))), float(np.max(abs(matrix.imag)))) < 2**48
    return matrix


def entry_bound(matrix):
    checked(matrix)
    return max(int(np.max(abs(matrix.real))), int(np.max(abs(matrix.imag))))


def mul(left, right):
    global MAX_PRODUCT_BOUND
    left, right = checked(left), checked(right)
    bound = 2*left.shape[1]*entry_bound(left)*entry_bound(right)
    assert bound < 2**48, ('integer-product guard', bound)
    MAX_PRODUCT_BOUND = max(MAX_PRODUCT_BOUND, bound)
    return checked(left @ right)


def add(*matrices):
    assert matrices
    bound = sum(entry_bound(matrix) for matrix in matrices)
    assert bound < 2**48
    return checked(sum(matrices))


def scale(number, matrix):
    assert int(number) == number
    assert abs(int(number))*entry_bound(matrix) < 2**48
    return checked(int(number)*matrix)


def equal(left, right, label):
    checked(left)
    checked(right)
    assert np.array_equal(left, right), label


def squared_integer_norm(matrix):
    matrix = checked(matrix)
    return sum(int(r)**2+int(i)**2 for r, i in zip(matrix.real.flat, matrix.imag.flat))


class CAR:
    def __init__(self, N):
        self.N = N
        self.cs = [checked(c) for c in fermions(N)]
        self.eye = checked(np.eye(1 << N))
        self.Q = checked(np.diag([state.bit_count() for state in range(1 << N)]))
        self.NminusQ = add(scale(N, self.eye), scale(-1, self.Q))
        self.sectors = {}

    def D(self, X):
        return add(mul(self.Q, X),
                   *[scale(-1, mul(mul(c.conj().T, X), c)) for c in self.cs])

    def Db(self, X):
        return add(mul(self.NminusQ, X),
                   *[scale(-1, mul(mul(c, X), c.conj().T)) for c in self.cs])

    def sector(self, m):
        if m not in self.sectors:
            ids = [i for i in range(1 << self.N) if i.bit_count() == m]
            C2 = casimir_matrix(self.N, m)
            Ps = projectors(C2, self.N, min(m, self.N-m))
            denominator = int(sp.ilcm(*[den for _, den in Ps]))
            rank_numerator = sum(rank*num*(denominator//den)
                                 for rank, (num, den) in enumerate(Ps))
            self.sectors[m] = (ids, checked(C2), denominator, checked(rank_numerator))
        return self.sectors[m]

    def check_word(self, X, degree, m):
        """All identities are linear in X, including at a finite resonance."""
        N, d = self.N, degree
        equal(mul(self.Q, X), mul(X, self.Q), 'neutral full-Fock word')
        DX, DbX = self.D(X), self.Db(X)
        DDbX, DbDX = self.D(DbX), self.Db(DX)
        equal(add(DDbX, scale(-1, DbDX)), add(DX, scale(-1, DbX)),
              '[D,Db]=D-Db, full-Fock order')
        AX, BX = add(scale(d, X), scale(-1, DX)), add(scale(d, X), scale(-1, DbX))
        ABX = add(scale(d, BX), scale(-1, self.D(BX)))
        ids, C2, den, rank_num = self.sector(m)
        restrict = lambda field: checked(field[np.ix_(ids, ids)])
        vector = lambda field: restrict(field).reshape(-1, 1)
        x = vector(X)
        casimir_rhs = add(scale(N-m+1, DX), scale(m, DbX), scale(-1, DDbX))
        equal(mul(C2, x), scale(2, vector(casimir_rhs)), 'general-charge Casimir')
        identity = checked(np.eye(len(x)))
        U_num = add(scale(d*den, identity), scale(-1, rank_num))
        Ux = mul(U_num, x)
        U2x = mul(U_num, Ux)
        # 2 den^2 E, E=(d-R)X - (A+B)X/2.
        error_num = add(scale(2*den, Ux), scale(-den*den, vector(add(AX, BX))))
        left = scale(N+1-2*d, error_num)
        right = add(scale(2*den*den, vector(ABX)), scale(-2, U2x),
                    scale(-(N-2*m+1)*den*den, vector(add(BX, scale(-1, AX)))))
        equal(left, right, 'cleared rank/size error identity')
        # Independent substitution of half-Clifford-size=(D+Db)/2.
        size_error_num = add(scale(-2*den, mul(rank_num, x)),
                             scale(den*den, vector(add(DX, DbX))))
        equal(error_num, size_error_num, 'error equals half-size minus rank')
        return dict(N=N, m=m, degree=degree, rank_denominator=den,
                    resonance_factor=N+1-2*d,
                    nonzero_sector_word=bool(np.any(x)),
                    nonzero_finite_bridge_error=bool(np.any(error_num)),
                    error_numerator_squared_norm=str(squared_integer_norm(error_num)),
                    error_denominator=2*den*den,
                    max_cleared_identity_entry=entry_bound(left), identities=4)


def sparse_hamiltonian(car, p, rng):
    basis = [checked(b) for b in creators(car.cs, p)]
    count = len(basis)
    J = np.zeros((count, count), dtype=np.complex128)
    for index in rng.choice(count, size=min(3, count), replace=False):
        J[index, index] += int(rng.choice([-1, 1]))
    pairs = list(combinations(range(count), 2))
    for which in rng.choice(len(pairs), size=min(5, len(pairs)), replace=False):
        a, b = pairs[which]
        value = int(rng.choice([-1, 1])) + 1j*int(rng.choice([-1, 1]))
        J[a, b] += value
        J[b, a] += value.conjugate()
    terms = []
    for a in range(count):
        for b in range(count):
            if J[a, b] != 0:
                # Complex Gaussian-integer coefficients remain exactly represented.
                term = mul(basis[a], basis[b].conj().T)
                terms.append(checked(J[a, b]*term))
    H = add(*terms)
    equal(H, H.conj().T, 'Hermitian sparse integer fixture')
    equal(car.D(H), scale(p, H), 'D H=p H')
    coupling_record = [[a, b, int(J[a, b].real), int(J[a, b].imag)]
                       for a in range(count) for b in range(count) if J[a, b] != 0]
    return H, coupling_record


def fixture_checks():
    rng = np.random.default_rng(SEED)
    rows, fixtures = [], []
    for N, p, sectors in [(5, 1, [2, 3]), (5, 2, [2, 3]),
                          (6, 1, [2, 4]), (6, 2, [2, 4])]:
        car = CAR(N)
        generated = [sparse_hamiltonian(car, p, rng) for _ in range(4)]
        fields = [pair[0] for pair in generated]
        fixtures.append(dict(N=N, p=p, sectors=sectors,
                             Gaussian_integer_couplings=[pair[1] for pair in generated]))
        for case in range(2):
            word = car.eye
            for j in range(1, 5):
                word = mul(word, fields[case])
                for m in sectors:
                    row = car.check_word(word, j*p, m)
                    row.update(j=j, case=case, word_kind='single_fixture_power')
                    rows.append(row)
        word = car.eye
        for j in range(1, 5):
            word = mul(word, fields[j-1])
            for m in sectors:
                row = car.check_word(word, j*p, m)
                row.update(j=j, word_kind='ordered_distinct_fixture_product')
                rows.append(row)
        print(f'general-j CAR fixtures N={N}, p={p} passed', flush=True)
    assert any(row['nonzero_finite_bridge_error'] for row in rows)
    assert any(row['resonance_factor'] == 0 for row in rows)
    return dict(rows=rows, fixtures=fixtures,
                scope='Sparse deterministic integer fixtures and ordered products, not Gaussian '
                      'ensemble sampling. Resonant finite cases check the cleared identity only.')


def highest_chaos_checks():
    N, p = 3, 1
    car = CAR(N)
    _, basis, variance_numerators = gue_basis(N, p)
    basis = [checked(b) for b in basis]
    count = len(basis)
    zero = (0,)*count
    coefficients = {zero: car.eye}
    summaries = []
    for j in range(1, 5):
        next_coefficients = {}
        for multiindex, value in coefficients.items():
            for coordinate, matrix in enumerate(basis):
                index = list(multiindex)
                index[coordinate] += 1
                index = tuple(index)
                term = mul(matrix, value)
                next_coefficients[index] = (add(next_coefficients[index], term)
                                            if index in next_coefficients else term)
        coefficients = next_coefficients
        expected = factorial(count+j-1)//(factorial(j)*factorial(count-1))
        assert len(coefficients) == expected
        for m in [1, 2]:
            nonzero_words = nonzero_errors = 0
            greatest_entry = 0
            for value in coefficients.values():
                record = car.check_word(value, j*p, m)
                nonzero_words += int(record['nonzero_sector_word'])
                nonzero_errors += int(record['nonzero_finite_bridge_error'])
                greatest_entry = max(greatest_entry, record['max_cleared_identity_entry'])
            summaries.append(dict(N=N, p=p, m=m, j=j, coordinates=count,
                                  homogeneous_coefficients=len(coefficients),
                                  nonzero_sector_coefficients=nonzero_words,
                                  nonzero_bridge_error_coefficients=nonzero_errors,
                                  identities=4*len(coefficients),
                                  max_cleared_identity_entry=greatest_entry,
                                  resonance_factor=N+1-2*j*p))
        print(f'complete highest-chaos coefficients degree {j} passed', flush=True)
    assert any(row['nonzero_bridge_error_coefficients'] for row in summaries)
    return dict(rows=summaries, Gaussian_coordinate_variances_over_common_w=
                [str(sp.Rational(int(r), 2)) for r in variance_numerators],
                scope='All nine real Gaussian coordinates of the complete N=3,p=1 model. '
                      'For each homogeneous coefficient g^alpha of H^j, replacing g^alpha '
                      'by its product monic Hermite polynomial gives W_j(H^j). The linear '
                      'operator identities are checked on every coefficient for j=1..4 '
                      'and m=1,2. This is not a test of every lower Hermite-chaos component.')


def clifford_checks():
    rows = []
    for N in [2, 3, 4]:
        car = CAR(N)
        gammas = []
        for c in car.cs:
            gammas.extend([add(c, c.conj().T), checked(1j*add(c.conj().T, scale(-1, c)))])
        for a, ga in enumerate(gammas):
            for b, gb in enumerate(gammas):
                equal(add(mul(ga, gb), mul(gb, ga)),
                      scale(2 if a == b else 0, car.eye), 'Majorana CAR')
        lengths = {}
        for length in range(0, 2*N+1, 2):
            count = 0
            for subset in combinations(range(2*N), length):
                word = car.eye
                for index in subset:
                    word = mul(word, gammas[index])
                size = add(car.D(word), car.Db(word))
                equal(size, scale(length, word), 'even Clifford degree equals D+Db')
                conjugated = add(*[mul(mul(gamma, word), gamma) for gamma in gammas])
                equal(scale(2, size), add(scale(2*N, word), scale(-1, conjugated)),
                      'size equals ungraded even-Clifford conjugation generator')
                count += 1
            lengths[str(length)] = count
        rows.append(dict(N=N, even_strings_by_length=lengths,
                         total_even_strings=sum(lengths.values())))
    return dict(rows=rows,
                convention='gamma_(2i)=c_i+c_i^dagger, gamma_(2i+1)=i(c_i^dagger-c_i); '
                           '{gamma_a,gamma_b}=2 delta_ab. On even Clifford strings of '
                           'length ell, (D+Db)X=ell X and (D+Db)X/2=(ell/2)X. '
                           'The ungraded formula is not asserted for odd strings.')


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
    start = time.time()
    fixtures = fixture_checks()
    chaos = highest_chaos_checks()
    clifford = clifford_checks()
    counts = dict(fixture_sector_word_cases=len(fixtures['rows']),
                  fixture_identities=sum(row['identities'] for row in fixtures['rows']),
                  highest_chaos_sector_coefficient_cases=sum(row['homogeneous_coefficients']
                                                            for row in chaos['rows']),
                  highest_chaos_identities=sum(row['identities'] for row in chaos['rows']),
                  even_clifford_strings=sum(row['total_even_strings'] for row in clifford['rows']))
    data = dict(status='passed', created_at_utc=datetime.now(timezone.utc).isoformat(),
                fixture_checks=fixtures, highest_chaos_checks=chaos, Clifford_size=clifford,
                counts=counts, seed=SEED,
                exact_error_identity='(N+1-2d)[(d-R)X-(A+B)X/2] = '
                                     'ABX-(d-R)^2X-(N-2m+1)(B-A)X/2, '
                                     'A=d-D, B=d-Db, d=jp; compositions precede restriction.',
                arithmetic=dict(kind='Gaussian integers stored exactly in complex128',
                                every_matrix_product_absolute_sum_bound_less_than=2**48,
                                largest_product_bound=MAX_PRODUCT_BOUND,
                                rank_projectors='Exact Casimir polynomial projectors; '
                                                'all rank denominators cleared before comparison.'),
                environment=dict(python=platform.python_version(), numpy=np.__version__,
                                 sympy=sp.__version__), input_sha256=source_hashes(),
                elapsed_seconds=time.time()-start,
                scope='Finite algebra and normalization checks only. The central strong bridge '
                      'requires the separate fixed-j overlap-moment and rank-tail proof; no '
                      'large-N convergence, growing-j estimate, or finite-time kernel is tested.',
                reproducibility='Only this new driver main was run. Output is exclusive-create; '
                                'all old sources and scientific outputs remain untouched.')
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open('x', encoding='utf-8') as handle:
        json.dump(data, handle, indent=2)
        handle.write('\n')
    print(json.dumps(dict(status='passed', counts=counts,
                          output_sha256=hashlib.sha256(OUT.read_bytes()).hexdigest(),
                          elapsed_seconds=data['elapsed_seconds'])), flush=True)


if __name__ == '__main__':
    main()
