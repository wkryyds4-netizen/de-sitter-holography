"""Whole-field nilpotent loss sources, with singleton normalization.

Ordinary words may have Gaussian pairs internal to each marked field.
This module is separate from the sealed no-internal-pair size module.
All arithmetic is Fraction; the callers impose Gaussian Wick exclusions.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import comb

from memory_size_marked_chords import mm, sign_of_grouping, I2, C2, CD2, ZERO


def mul(a, b):
    out = Counter()
    for i, x in a.items():
        for j, y in b.items():
            if not i & j:
                out[i | j] += x*y
    return {i: x for i, x in out.items() if x}


def inverse(a, marks):
    assert a.get(0, 0)
    constant = a[0]
    rest = {i: -x/constant for i, x in a.items() if i}
    out, power = Counter({0: F(1)}), {0: F(1)}
    for _ in range(marks):
        power = mul(power, rest)
        out.update(power)
    return {i: x/constant for i, x in out.items() if x}


def loss(x, kind):
    a, b, c, d = x
    if kind == 'D':
        return (F(0), b/2, c/2, d-a)
    if kind == 'L':
        return ((a-d)/2, b/2, c/2, (d-a)/2)
    raise ValueError(kind)


def local(length, pi, blocks, occupied):
    starts = {a: (index, b, kind) for index, (a, b, kind) in enumerate(blocks)}
    result, pos = {0: I2}, 0
    while pos < length:
        mark, end, kind = starts.get(pos, (None, pos+1, None))
        value = I2
        for slot in range(pos, end):
            if slot in occupied:
                value = mm(value, CD2)
            if pi[slot] in occupied:
                value = mm(value, C2)
        factor = {0: value}
        if mark is not None:
            factor[1 << mark] = tuple(-x for x in loss(value, kind))
        out = {}
        for i, x in result.items():
            for j, y in factor.items():
                if not i & j:
                    term = mm(x, y)
                    out[i | j] = tuple(a+b for a, b in zip(out.get(i | j, ZERO), term))
        result, pos = out, end
    return {mask: (x[0]+x[3])/2 for mask, x in result.items() if x[0]+x[3]}


@lru_cache(None)
def exponent(length, pairing, blocks):
    """Return K in exp(lambda K); p*chaos singleton terms cancel first.

    D blocks permit internal pairs. L blocks here are limited to Wick
    source blocks; their mark is the established incidence/size deficit.
    The constant coefficient is -4 times the chord crossing number.
    """
    assert length % 2 == 0
    assert sorted(x for pair in pairing for x in pair) == list(range(length))
    assert all(0 <= a < b <= length for a, b, _ in blocks)
    assert all(blocks[i][1] <= blocks[i+1][0] for i in range(len(blocks)-1))
    assert all(kind in ['D', 'L'] for _, _, kind in blocks)
    assert not any(kind == 'L' and a <= x < y < b
                   for a, b, kind in blocks for x, y in pairing)
    pi = dict(pairing) | {b: a for a, b in pairing}
    singletons = [local(length, pi, blocks, (i,)) for i in range(length)]
    product = {0: F(1)}
    for single in singletons:
        product = mul(product, {mask: 2*value for mask, value in single.items()})
    degrees = []
    for index, (a, b, _) in enumerate(blocks):
        degree = b-a-2*sum(a <= x < y < b for x, y in pairing)
        degrees.append(degree)
        product = mul(product, {0: F(1), 1 << index: F(degree)})
    assert product == {0: F(1)}, ('singleton cancellation', pairing, blocks, product)
    answer = Counter({0: F(-comb(length, 2)-length//2)})
    for i, j in combinations(range(length), 2):
        labels = list(range(length))
        labels[j] = i
        sequence = [labels[x] for slot in range(length) for x in (slot, pi[slot])]
        sign = sign_of_grouping(sequence)
        numerator = {mask: sign*value for mask, value in local(length, pi, blocks, (i, j)).items()}
        denominator = mul(singletons[i], singletons[j])
        answer.update(mul(numerator, inverse(denominator, len(blocks))))
    crossing = sum(a < c < b < d or c < a < d < b
                   for (a, b), (c, d) in combinations(pairing, 2))
    assert answer[0] == -4*crossing
    return {mask: value for mask, value in answer.items() if value}


def G_exponent(length, pairing, blocks):
    out = Counter(exponent(length, pairing, blocks))
    for index, (a, b, kind) in enumerate(blocks):
        if kind == 'D':
            degree = b-a-2*sum(a <= x < y < b for x, y in pairing)
            out[1 << index] -= degree*(degree-1)
    return {mask: value for mask, value in out.items() if value}
