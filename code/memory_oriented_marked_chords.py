"""Exact local oriented deficits on neutral, coupling-Wick source blocks.

D is Q X - sum c_i^* X c_i; Db is its opposite orientation, and
L=(D+Db)/2. Blocks must have no internal Gaussian pairing. This is NOT
the whole-field Gaussian-degree mark on primitives with internal pairs.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import comb

I2 = (F(1), F(0), F(0), F(1))
C2 = (F(0), F(1), F(0), F(0))
CD2 = (F(0), F(0), F(1), F(0))
ZERO = (F(0),) * 4
KINDS = ("D", "L", "Db")


def mm(a, b):
    return (a[0]*b[0]+a[1]*b[2], a[0]*b[1]+a[1]*b[3],
            a[2]*b[0]+a[3]*b[2], a[2]*b[1]+a[3]*b[3])


def loss(x, kind, form="asymmetric"):
    """On globally even source words; local factors themselves may be odd.

    Jordan-Wigner grouping signs remain outside the local product. The
    graded sandwich has sign (-1)^local_parity, but vanishes for odd local
    factors. Global odd words require a different treatment.
    """
    a, b, c, d = x
    if form == "asymmetric":
        normal, reverse = (F(0), F(0), c, d-a), (a-d, b, F(0), F(0))
    elif form == "symmetric":
        normal, reverse = (F(0), b/2, c/2, d-a), (a-d, b/2, c/2, F(0))
    else:
        raise ValueError(form)
    if kind == "D":
        return normal
    if kind == "Db":
        return reverse
    if kind == "L":
        return tuple((u+v)/2 for u, v in zip(normal, reverse))
    raise ValueError(kind)


def deficit(x, creation_count, annihilation_count, kind, form="asymmetric"):
    if form == "symmetric" or kind == "L":
        scalar = F(creation_count+annihilation_count, 2)
    else:
        scalar = creation_count if kind == "D" else annihilation_count
    return tuple(scalar*u-v for u, v in zip(x, loss(x, kind, form)))


def poly_multiply(a, b):
    out = {}
    for mask1, x in a.items():
        for mask2, y in b.items():
            if mask1 & mask2:
                continue
            mask = mask1 | mask2
            out[mask] = tuple(u+v for u, v in zip(out.get(mask, ZERO), mm(x, y)))
    return out


def local_polynomial(length, pi, blocks, kinds, occupied, form="asymmetric"):
    starts = {a: (mark, b) for mark, (a, b) in enumerate(blocks)}
    result, pos = {0: I2}, 0
    while pos < length:
        mark, end = starts.get(pos, (None, pos+1))
        matrix, nc, na = I2, 0, 0
        for slot in range(pos, end):
            if slot in occupied:
                matrix, nc = mm(matrix, CD2), nc+1
            if pi[slot] in occupied:
                matrix, na = mm(matrix, C2), na+1
        factor = {0: matrix}
        if mark is not None:
            factor[1 << mark] = deficit(matrix, nc, na, kinds[mark], form)
        result = poly_multiply(result, factor)
        pos = end
    return {mask: (mat[0]+mat[3])/2 for mask, mat in result.items() if mat[0]+mat[3]}


def grouping_sign(sequence):
    return (-1) ** sum(sequence[i] > sequence[j]
                       for i in range(len(sequence)) for j in range(i+1, len(sequence)))


def validate(length, pairing, blocks, kinds):
    assert len(blocks) == len(kinds) and all(kind in KINDS for kind in kinds)
    assert all(0 <= a < b <= length for a, b in blocks)
    assert all(blocks[i][1] <= blocks[i+1][0] for i in range(len(blocks)-1))
    assert sorted(x for chord in pairing for x in chord) == list(range(length))
    assert not any(a <= x < b and a <= y < b
                   for a, b in blocks for x, y in pairing)


def pair_site_terms(length, pairing, blocks, kinds, form="asymmetric"):
    validate(length, pairing, blocks, kinds)
    pi = dict(pairing) | {b: a for a, b in pairing}
    rows = []
    for i, j in combinations(range(length), 2):
        labels = list(range(length))
        labels[j] = labels[i]
        sequence = [labels[x] for slot in range(length) for x in (slot, pi[slot])]
        sign = grouping_sign(sequence)
        local = local_polynomial(length, pi, blocks, kinds, {i, j}, form)
        rows.append(((i, j), {mask: 4*sign*value for mask, value in local.items()}))
    return rows


@lru_cache(None)
def exponent(length, pairing, blocks, kinds, form="asymmetric"):
    result = Counter({0: F(-comb(length, 2)-length//2)})
    for _, local in pair_site_terms(length, pairing, blocks, kinds, form):
        result.update(local)
    return dict(result)


def crosses(a, b):
    x, y = sorted(a)
    z, w = sorted(b)
    return (x < z < y < w) or (z < x < w < y)


def crossing_count(chords):
    return sum(crosses(a, b) for a, b in combinations(chords, 2))


def geometric_exponent(length, pairing, blocks, kinds):
    validate(length, pairing, blocks, kinds)
    incident = [frozenset(chord for chord in pairing if any(a <= x < b for x in chord))
                for a, b in blocks]
    result = {0: -4*crossing_count(pairing)}
    t = {"D": F(0), "L": F(1, 2), "Db": F(1)}
    for i, chords in enumerate(incident):
        j = len(chords)
        result[1 << i] = j*(j-1+2*t[kinds[i]])-4*crossing_count(chords)
    for i, j in combinations(range(len(blocks)), 2):
        common = incident[i] & incident[j]
        m, cr = len(common), crossing_count(common)
        u, v = t[kinds[i]], t[kinds[j]]
        same_orientation = (1-u)*(1-v)+u*v
        result[(1 << i) | (1 << j)] = 4*same_orientation*comb(m, 2)-4*cr+4*u*v*m
    return result


def finite_local(length, sets, pairing, blocks, kinds, number_sites, form="asymmetric"):
    validate(length, pairing, blocks, kinds)
    pi = dict(pairing) | {b: a for a, b in pairing}
    sequence = [site for slot in range(length)
                for site in (*sets[slot], *reversed(sets[pi[slot]]))]
    result = {0: F(grouping_sign(sequence))}
    for site in range(number_sites):
        occupied = {slot for slot, sites in enumerate(sets) if site in sites}
        local = local_polynomial(length, pi, blocks, kinds, occupied, form)
        updated = Counter()
        for a, x in result.items():
            for b, y in local.items():
                if not a & b:
                    updated[a | b] += x*y
        result = dict(updated)
    return [result.get(mask, F(0)) for mask in range(1 << len(blocks))]
