"""Ordered Gaussian-Wick descendants with separate Wick, mark and OU domains.

Formal moments use C=1.  A field carries its physical energy degree j;
<F_j^dagger H^a G_k H^b> restores C^((j+k+a+b)/2).  A term's implicit
C power is (physical_degree-raw_slot_count)/2.  No pinching is performed.
"""
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations

from memory_size_marked_chords import exponent, coefficient_polynomial
from singlet_coefficient import pairings


def clean(poly):
    return {key: F(value) for key, value in poly.items() if value}


def padd(left, right, factor=1):
    out = Counter(left)
    for key, value in right.items():
        out[key] += factor*value
    return clean(out)


def pmul(left, right):
    out = Counter()
    for (a, b), x in left.items():
        for (c, d), y in right.items():
            out[a+c, b+d] += x*y
    return clean(out)


def mono(lam=0, q=0, coefficient=1):
    return {(lam, q): F(coefficient)} if coefficient else {}


def qinteger(j):
    return {(0, r): F(1) for r in range(j)}


@dataclass(frozen=True)
class Term:
    length: int
    wick_blocks: tuple = ()
    mark_blocks: tuple = ()
    ou_power: int = 0
    selected_chaos: int = -1

    def dagger(self):
        reflect = lambda blocks: tuple(sorted((self.length-b, self.length-a) for a, b in blocks))
        return Term(self.length, reflect(self.wick_blocks), reflect(self.mark_blocks),
                    self.ou_power, self.selected_chaos)

    def left_h(self):
        assert self.ou_power == 0 and self.selected_chaos == -1
        shift = lambda blocks: tuple((a+1, b+1) for a, b in blocks)
        return Term(self.length+1, shift(self.wick_blocks), shift(self.mark_blocks))

    def right_h(self):
        assert self.ou_power == 0 and self.selected_chaos == -1
        return Term(self.length+1, self.wick_blocks, self.mark_blocks)


@dataclass
class Field:
    degree: int
    terms: dict
    name: str = ''

    def validate(self):
        for term in self.terms:
            assert 0 <= term.length <= self.degree and (self.degree-term.length) % 2 == 0
            assert all(0 <= a < b <= term.length for a, b in term.wick_blocks+term.mark_blocks)
            assert all(any(a <= x < y <= b for a, b in term.wick_blocks)
                       for x, y in term.mark_blocks)
        return self


def fadd(left, right, factor=1, name=''):
    assert left.degree == right.degree
    terms = {term: dict(poly) for term, poly in left.terms.items()}
    for term, poly in right.terms.items():
        value = padd(terms.get(term, {}), poly, factor)
        if value:
            terms[term] = value
        else:
            terms.pop(term, None)
    return Field(left.degree, terms, name).validate()


def fscale(field, polynomial, degree=None, name=''):
    degree = field.degree if degree is None else degree
    return Field(degree, {term: pmul(poly, polynomial) for term, poly in field.terms.items()
                         if pmul(poly, polynomial)}, name).validate()


def hleft(field):
    return Field(field.degree+1, {term.left_h(): poly for term, poly in field.terms.items()}).validate()


def hright(field):
    return Field(field.degree+1, {term.right_h(): poly for term, poly in field.terms.items()}).validate()


def ou(field):
    return Field(field.degree,
                 {Term(term.length, term.wick_blocks, term.mark_blocks, term.ou_power+1,
                       term.selected_chaos): poly for term, poly in field.terms.items()}).validate()


def select(field, chaos):
    return Field(field.degree,
                 {Term(term.length, term.wick_blocks, term.mark_blocks, term.ou_power, chaos): poly
                  for term, poly in field.terms.items()}).validate()


def identity(degree=0):
    return Field(degree, {Term(0): mono()}, 'I').validate()


def Z(j):
    return identity() if j == 0 else Field(j, {Term(j, ((0, j),)): mono()}, f'Z{j}').validate()


def composite(v, w):
    assert v >= 1 and w >= 0
    n = v+w
    return Field(n, {Term(n, ((0, n),), ((0, v),)): mono()}, f'U{v}_diamond_Z{w}').validate()


def U(j):
    return composite(j, 0) if j else Field(0, {}, 'U0')


def descendants(maximum):
    eta = [Field(j, {}, f'eta{j}') for j in range(maximum+1)]
    sources = {}
    for j in range(2, maximum):
        source = Field(j+1, {}, f'rhat{j}')
        for v in range(1, j):
            source = fadd(source, fscale(composite(v, j-1-v), mono(1, v, 4), degree=j+1))
        sources[j] = source
        value = hleft(eta[j])
        value = fadd(value, fscale(eta[j-1], qinteger(j), degree=j+1), -1)
        eta[j+1] = fadd(value, source, -1, f'eta{j+1}')
    V = [Field(0, {}, 'V0')]
    for j in range(1, maximum+1):
        value = fscale(U(j), mono(1, j, 4))
        value = fadd(value, fscale(eta[j], mono(0, j)))
        V.append(fadd(value, ou(eta[j]), -1, f'V{j}'))
    return eta, V, sources


def allowed(pairing, wick_blocks):
    return not any(a <= x < b and a <= y < b
                   for a, b in wick_blocks for x, y in pairing)


def outgoing_degree(pairing, interval):
    a, b = interval
    return b-a-2*sum(a <= x < b and a <= y < b for x, y in pairing)


def layout(left, right, energy_a, energy_b):
    left = left.dagger()
    offset = left.length+energy_a
    length = offset+right.length+energy_b
    shift = lambda blocks: tuple((a+offset, b+offset) for a, b in blocks)
    wick_blocks = left.wick_blocks+shift(right.wick_blocks)
    mark_blocks = left.mark_blocks+shift(right.mark_blocks)
    return left, length, wick_blocks, mark_blocks, ((0, left.length), (offset, offset+right.length))


def crosses(a, b):
    x, y = sorted(a)
    u, v = sorted(b)
    return (x < u < y < v) or (u < x < v < y)


def geometric_kappa(pairing, marks):
    count = lambda chords: sum(crosses(a, b) for a, b in combinations(chords, 2))
    incident = [tuple(edge for edge in pairing if any(a <= endpoint < b for endpoint in edge))
                for a, b in marks]
    out = {0: -4*count(pairing)}
    for index, (a, b) in enumerate(marks):
        assert len(incident[index]) == b-a
        out[1 << index] = (b-a)**2-4*count(incident[index])
    for a, b in combinations(range(len(marks)), 2):
        bridges = tuple(edge for edge in incident[a] if edge in incident[b])
        out[(1 << a) | (1 << b)] = len(bridges)**2-4*count(bridges)
    return out


@lru_cache(None)
def term_moment(left, right, energy_a=0, energy_b=0, no_field_bridge=False):
    dual, length, wicks, marks, intervals = layout(left, right, energy_a, energy_b)
    if length % 2:
        return (), 0, ()
    result = Counter()
    histogram = Counter()
    count = 0
    for pairing in pairings(tuple(range(length))):
        if not allowed(pairing, wicks):
            continue
        if no_field_bridge:
            (la, lb), (ra, rb) = intervals
            if any((la <= x < lb and ra <= y < rb) or
                   (la <= y < lb and ra <= x < rb) for x, y in pairing):
                continue
        dleft, dright = [outgoing_degree(pairing, interval) for interval in intervals]
        if ((left.selected_chaos >= 0 and dleft != left.selected_chaos) or
                (right.selected_chaos >= 0 and dright != right.selected_chaos)):
            continue
        kappa = exponent(length, pairing, marks)
        geometry = geometric_kappa(pairing, marks)
        assert all(kappa.get(mask, 0) == geometry.get(mask, 0) for mask in range(1 << len(marks)))
        crossings = -kappa[0]/4
        assert crossings.denominator == 1
        power_q = int(crossings)+left.ou_power*dleft+right.ou_power*dright
        for power_lambda, value in coefficient_polynomial(kappa, len(marks)).items():
            result[power_lambda, power_q] += value
        histogram[dleft, dright] += 1
        count += 1
    return tuple(sorted(clean(result).items())), count, tuple(sorted(histogram.items()))


def moment(left, right, energy_a=0, energy_b=0, max_slots=None, no_field_bridge=False):
    result = {}
    actual_slots = 0
    diagram_calls = 0
    for term_left, poly_left in left.terms.items():
        for term_right, poly_right in right.terms.items():
            slots = term_left.length+term_right.length+energy_a+energy_b
            if max_slots is not None:
                assert slots <= max_slots
            actual_slots = max(actual_slots, slots)
            polynomial, count, _ = term_moment(term_left, term_right, energy_a, energy_b, no_field_bridge)
            result = padd(result, pmul(pmul(poly_left, poly_right), dict(polynomial)))
            diagram_calls += count
    return clean(result), dict(maximum_raw_slots=actual_slots, allowed_term_pairings=diagram_calls,
                              energy_degree=left.degree+right.degree+energy_a+energy_b)


def serialize_polynomial(poly):
    return [dict(lambda_power=a, q_power=b, coefficient=str(value))
            for (a, b), value in sorted(poly.items())]


def serialize_field(field):
    return dict(name=field.name, physical_degree=field.degree,
                terms=[dict(raw_slots=term.length, implicit_C_power=(field.degree-term.length)//2,
                            wick_blocks=term.wick_blocks, mark_blocks=term.mark_blocks,
                            OU_q_power=term.ou_power, selected_chaos=term.selected_chaos,
                            coefficient=serialize_polynomial(poly))
                       for term, poly in sorted(field.terms.items(), key=lambda item: repr(item[0]))])
