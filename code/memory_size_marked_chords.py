"""Local Clifford-size deficit marks in common-J Wick chord diagrams.

Fractions and nilpotent block variables retain signs and same-site marks.
The finite CAR checks are independent of the pair-Poisson limiting rule.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import datetime
import hashlib
import json
import math
import platform
import numpy as np

from singlet_coefficient import pairings
from singlet_multimark import crossing, nmul
from complex_syk_dictionary import fermions, creators

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'results/susskind-memory-size-bridge'
I2 = (F(1), F(0), F(0), F(1))
C2 = (F(0), F(1), F(0), F(0))
CD2 = (F(0), F(0), F(1), F(0))
ZERO = (F(0),) * 4


def mm(a, b):
    return (a[0]*b[0]+a[1]*b[2], a[0]*b[1]+a[1]*b[3],
            a[2]*b[0]+a[3]*b[2], a[2]*b[1]+a[3]*b[3])


def deficit(a, incidences):
    halfsize = ((a[0]-a[3])/2, a[1]/2, a[2]/2, (a[3]-a[0])/2)
    return tuple(F(incidences, 2)*x-y for x, y in zip(a, halfsize))


def multiply_polynomials(a, b):
    result = {}
    for i, x in a.items():
        for j, y in b.items():
            if i & j:
                continue
            term = mm(x, y)
            result[i | j] = tuple(u+v for u, v in zip(result.get(i | j, ZERO), term))
    return result


def local_polynomial(length, pi, blocks, occupied_labels):
    occupied = set(occupied_labels)
    starts = {a: (i, b) for i, (a, b) in enumerate(blocks)}
    result = {0: I2}
    pos = 0
    while pos < length:
        if pos in starts:
            mark, end = starts[pos]
        else:
            mark, end = None, pos+1
        matrix = I2
        incidences = 0
        for slot in range(pos, end):
            if slot in occupied:
                matrix = mm(matrix, CD2)
                incidences += 1
            if pi[slot] in occupied:
                matrix = mm(matrix, C2)
                incidences += 1
        factor = {0: matrix}
        if mark is not None:
            factor[1 << mark] = deficit(matrix, incidences)
        result = multiply_polynomials(result, factor)
        pos = end
    return {mask: (mat[0]+mat[3])/2 for mask, mat in result.items()
            if mat[0]+mat[3]}


def sign_of_grouping(sequence):
    return (-1) ** sum(sequence[i] > sequence[j]
                       for i in range(len(sequence)) for j in range(i+1, len(sequence)))


def validate_blocks(length, pairing, blocks):
    assert all(0 <= a < b <= length for a, b in blocks)
    assert all(blocks[i][1] <= blocks[i+1][0] for i in range(len(blocks)-1))
    assert not any(a <= x < b and a <= y < b
                   for a, b in blocks for x, y in pairing)


@lru_cache(None)
def exponent(length, pairing, blocks):
    """Kappa in exp(lambda*Kappa); constant is exactly -4*crossings."""
    validate_blocks(length, pairing, blocks)
    pi = dict(pairing) | {b: a for a, b in pairing}
    answer = Counter({0: F(-math.comb(length, 2)-length//2)})
    for i, j in combinations(range(length), 2):
        labels = list(range(length))
        labels[j] = labels[i]
        sequence = [labels[x] for slot in range(length) for x in (slot, pi[slot])]
        sign = sign_of_grouping(sequence)
        local = local_polynomial(length, pi, blocks, (i, j))
        ratio = {mask: 4*sign*value for mask, value in local.items()}
        # Independent determinant minor is valid also when the local trace is zero.
        ent = lambda a, b: 1 if a <= pi[b] else -1
        minor = (ent(i, i)*ent(j, j)-ent(i, j)*ent(j, i))*ent(i, i)*ent(j, j)
        assert ratio.get(0, F(0)) == minor
        for mask, value in ratio.items():
            answer[mask] += value
    cr = sum(crossing(a, b) for a, b in combinations(pairing, 2))
    assert answer[0] == -4*cr
    return dict(answer)


def coefficient_polynomial(kappa, number_marks):
    """Coefficients by lambda degree, with q^crossings stripped off."""
    nonconstant = {k: v for k, v in kappa.items() if k and v}
    term = {0: F(1)}
    full = (1 << number_marks)-1
    result = {}
    for power in range(number_marks+1):
        value = term.get(full, F(0))/math.factorial(power)
        if value:
            result[power] = value
        term = nmul(term, nonconstant)
    return result


def moment(j, k, left_energy=0, right_energy=0, mark_mask=3):
    length = j+left_energy+k+right_energy
    if length % 2:
        return {}, 0
    blocks_all = ((0, j), (j+left_energy, j+left_energy+k))
    marked = tuple(block for a, block in enumerate(blocks_all) if mark_mask >> a & 1)
    result = Counter()
    count = 0
    for pairing in pairings(tuple(range(length))):
        if any(a <= x < b and a <= y < b for a, b in blocks_all for x, y in pairing):
            continue
        kappa = exponent(length, pairing, marked)
        cr = int(-kappa[0]/4)
        for power, value in coefficient_polynomial(kappa, len(marked)).items():
            result[power, cr] += value
        count += 1
    return {key: value for key, value in result.items() if value}, count


def finite_local(N, p, sets, pairing, blocks):
    length = len(sets)
    validate_blocks(length, pairing, blocks)
    pi = dict(pairing) | {b: a for a, b in pairing}
    sequence = [site for slot in range(length)
                for site in (*sets[slot], *reversed(sets[pi[slot]]))]
    sign = sign_of_grouping(sequence)
    result = {0: F(sign)}
    for site in range(N):
        occupied = tuple(slot for slot, values in enumerate(sets) if site in values)
        result = nmul(result, local_polynomial(length, pi, blocks, occupied))
    return [result.get(mask, F(0)) for mask in range(1 << len(blocks))]


def finite_checks():
    rng = np.random.default_rng(110709)
    rows = []
    for N, p, j, k, gap in [(3, 1, 1, 1, 0), (4, 2, 2, 2, 0),
                            (4, 1, 3, 3, 0), (5, 2, 2, 2, 2),
                            (5, 2, 1, 3, 2)]:
        length = j+k+gap
        blocks = ((0, j), (j+gap, length))
        allowed = [pa for pa in pairings(tuple(range(length))) if not any(
            a <= x < b and a <= y < b for a, b in blocks for x, y in pa)]
        cs = fermions(N)
        subsets = list(combinations(range(N), p))
        creation = dict(zip(subsets, creators(cs, p)))
        eye = np.eye(1 << N)
        half_size = lambda X: N*X/2-sum(c.conj().T@X@c+c@X@c.conj().T for c in cs)/2
        for case in range(32):
            pairing = allowed[case % len(allowed)]
            pi = dict(pairing) | {b: a for a, b in pairing}
            sets = ([subsets[case % len(subsets)]]*length if case < 8 else
                    [subsets[int(rng.integers(len(subsets)))] for _ in range(length)])
            ops = [creation[sets[i]]@creation[sets[pi[i]]].conj().T for i in range(length)]
            starts = {a: (ii, b) for ii, (a, b) in enumerate(blocks)}
            actual = []
            for mask in range(4):
                prod = eye.copy()
                pos = 0
                while pos < length:
                    if pos in starts:
                        mark, end = starts[pos]
                    else:
                        mark, end = None, pos+1
                    block = eye.copy()
                    for slot in range(pos, end):
                        block = block@ops[slot]
                    if mark is not None and mask >> mark & 1:
                        block = (end-pos)*p*block-half_size(block)
                    prod = prod@block
                    pos = end
                trace = np.trace(prod)
                scaled = trace*4
                assert scaled.imag == 0 and scaled.real == round(scaled.real)
                actual.append(F(int(scaled.real), 4*(1 << N)))
            expected = finite_local(N, p, sets, pairing, blocks)
            assert actual == expected, (N, p, blocks, sets, pairing, actual, expected)
            rows.append(dict(N=N, p=p, blocks=blocks, sets=sets, pairing=pairing,
                             all_mark_subsets=list(map(str, actual))))
    return rows


def serialize_polynomial(poly):
    return [dict(lambda_power=i, q_power=j, coefficient=str(value))
            for (i, j), value in sorted(poly.items())]


def main():
    path = OUT / 'marked-chord-checks.json'
    if path.exists():
        raise FileExistsError(path)
    finite = finite_checks()
    print('finite local/global CAR comparisons passed', flush=True)
    moments = []
    for j, k in [(1, 1), (2, 2), (3, 3), (4, 4), (1, 3), (2, 4)]:
        for a, b in [(0, 0), (1, 1), (2, 0)]:
            if j+k+a+b > 10:
                continue
            for mask in [0, 1, 3]:
                poly, count = moment(j, k, a, b, mask)
                moments.append(dict(j=j, k=k, left_energy=a, right_energy=b,
                                    mark_mask=mask, allowed_pairings=count,
                                    polynomial=serialize_polynomial(poly)))
    assert moment(1, 1, mark_mask=0)[0] == {(0, 0): F(1)}
    assert moment(1, 1, mark_mask=1)[0] == {(1, 0): F(1)}
    assert moment(1, 1, mark_mask=3)[0] == {(1, 0): F(1), (2, 0): F(1)}
    assert moment(2, 2, mark_mask=1)[0] == {(1, 0): F(4)}
    names = ['code/memory_size_marked_chords.py', 'code/singlet_coefficient.py',
             'code/singlet_multimark.py', 'code/complex_syk_dictionary.py']
    data = dict(status='passed', created_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                seed=110709, finite_words=len(finite), finite_comparisons=4*len(finite),
                finite_records=finite, moment_rows=moments,
                scope='Exact finite-word CAR versus local graded tensor algebra for all subsets of two size marks. Limiting chord polynomials use the stated pair-Poisson theorem; they are not finite-size equalities or a closed finite-time spectral kernel. Rank marks require the separate strong bridge.',
                normalization='Moments are divided by C^((j+k+a+b)/2), lambda=p^2/N and q=exp(-4lambda). Marked blocks are coupling-Wick ordered; internal chords are excluded for both marked and unmarked Wick blocks.',
                environment=dict(python=platform.python_version(), numpy=np.__version__),
                input_sha256={name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names})
    OUT.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8') as handle:
        handle.write(json.dumps(data, indent=2)+'\n')
    print(json.dumps(dict(status='passed', finite_words=len(finite),
                          comparisons=4*len(finite), moment_cases=len(moments),
                          cached_exponents=exponent.cache_info().currsize)), flush=True)


if __name__ == '__main__':
    main()
