"""Exact finite CAR checks for differentiated marked-Wick algebra.

All matrix-unit checks use int64 after multiplying size maps by two.
These are algebraic fixtures, not samples of the Gaussian ensemble.
No successful scientific output is overwritten.
"""
from pathlib import Path
import hashlib
import itertools
import json
import platform

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/susskind-sector-gradient-tail/graded-checks.json'


def setup(n):
    dim = 1 << n
    q = np.array([v.bit_count() for v in range(dim)], dtype=np.int64)
    c = []
    for i in range(n):
        a = np.zeros((dim, dim), dtype=np.int64)
        for v in range(dim):
            if v & (1 << i):
                a[v ^ (1 << i), v] = (-1) ** ((v & ((1 << i) - 1)).bit_count())
        c.append(a)
    return q, c


def partial(a, x, parity):
    return a @ x - (-1) ** parity * x @ a


def size2(x, parity, q, c, which):
    d = (q[:, None] + q[None, :]) * x - 2 * (-1) ** parity * sum(a.T @ x @ a for a in c)
    if which == 'L':
        d = d - sum(partial(a.T, partial(a, x, parity), 1-parity) for a in c)
    return d


def main():
    if OUT.exists():
        raise SystemExit('Saved scientific output exists; do not rerun.')
    counts = {'matrix_units': 0, 'graded_derivative_commutators': 0,
              'repeated_shift_identities': 0, 'normal_monomial_euler': 0,
              'normal_monomial_deletions': 0, 'nonzero_deletions': 0}
    max_integer = 0
    rows = []

    def eq(a, b, context):
        nonlocal max_integer
        max_integer = max(max_integer, int(np.max(np.abs(a))), int(np.max(np.abs(b))))
        assert np.array_equal(a, b), context

    for n in range(1, 5):
        before = dict(counts)
        q, c = setup(n)
        dim = 1 << n
        derivs = [(i, False, a) for i, a in enumerate(c)] + [(i, True, a.T) for i, a in enumerate(c)]
        for u, v in itertools.product(range(dim), repeat=2):
            x = np.zeros((dim, dim), dtype=np.int64)
            x[u, v] = 1
            parity = (int(q[u]) + int(q[v])) % 2
            counts['matrix_units'] += 1
            for which in ('D', 'L'):
                dx = size2(x, parity, q, c, which)
                for i, raising, a in derivs:
                    y = partial(a, x, parity)
                    eq(partial(a, dx, parity) - size2(y, 1-parity, q, c, which), y,
                       ('comm', n, u, v, which, i, raising))
                    counts['graded_derivative_commutators'] += 1

        # All normal-ordered monomials, including I intersect J and odd degree.
        # The list order, rather than a reordered sign convention, fixes signs.
        for ibits, jbits in itertools.product(range(dim), repeat=2):
            letters = [(i, True, c[i].T) for i in range(n) if ibits & (1 << i)]
            letters += [(j, False, c[j]) for j in range(n) if jbits & (1 << j)]
            degree = len(letters)
            parity = degree % 2

            def product(omitted=None):
                a = np.eye(dim, dtype=np.int64)
                for t, (_, _, op) in enumerate(letters):
                    if t != omitted:
                        a = a @ op
                return a

            x = product()
            eq(size2(x, parity, q, c, 'D'), degree * x, ('Euler', n, ibits, jbits))
            counts['normal_monomial_euler'] += 1
            for i, raising, a in derivs:
                expected = np.zeros_like(x)
                for t, (site, creation, _) in enumerate(letters):
                    if site == i and creation != raising:
                        expected += (-1) ** t * product(t)
                actual = partial(a, x, parity)
                eq(actual, expected, ('deletion', n, ibits, jbits, i, raising))
                counts['normal_monomial_deletions'] += 1
                counts['nonzero_deletions'] += int(np.any(actual))

        # Non-selfadjoint neutral fixtures with two/three repeated deficits.
        rng = np.random.default_rng(20261108+n)
        for trial in range(3):
            x = rng.integers(-2, 3, (dim, dim), dtype=np.int64)
            x[q[:, None] != q[None, :]] = 0
            b = 1 + trial
            for k in (1, 2, 3):
                for word in itertools.product(('D', 'L'), repeat=k):
                    left = x.copy()
                    for which in reversed(word):
                        left = 2*b*left - size2(left, 0, q, c, which)
                    for i, raising, a in derivs:
                        right = partial(a, x, 0)
                        for which in reversed(word):
                            right = (2*b-1)*right - size2(right, 1, q, c, which)
                        eq(partial(a, left, 0), right, ('repeat', n, trial, word, i, raising))
                        counts['repeated_shift_identities'] += 1
        rows.append({'N': n, 'counts': {k: counts[k]-before[k] for k in counts}})

    assert max_integer < 10**12  # all intermediates have much smaller conservative bounds
    result = {'status': 'passed', 'counts': counts, 'rows': rows,
              'parameters': {'N': [1, 2, 3, 4], 'repeated_deficits_max': 3,
                             'all_matrix_units': True, 'all_normal_monomials': True},
              'environment': {'python': platform.python_version(), 'numpy': np.__version__,
                              'arithmetic': 'exact int64; doubled D_g/L_g maps'},
              'max_observed_integer': max_integer,
              'input_sha256': {'code/check_sector_graded_gradient.py': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
              'scope': 'Exact graded CAR derivation commutators, half-unit deficit shift and membership-preserving monomial deletion. Does not numerically establish degree-uniform Gaussian estimates.'}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open('x', encoding='utf-8') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    print(json.dumps({'status': 'passed', 'counts': counts,
                      'sha256': hashlib.sha256(OUT.read_bytes()).hexdigest()}))


if __name__ == '__main__':
    main()
