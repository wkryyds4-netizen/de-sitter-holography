"""Pure exact sparse CAR and scalar helpers for covariance path tangents."""
from collections import defaultdict
from fractions import Fraction as F
from math import comb, prod


def clean(x):
    return {k: v for k, v in x.items() if v}


def add(x, y, scale=F(1)):
    z = defaultdict(F, x)
    for key, value in y.items():
        z[key] += scale * value
    return clean(z)


def scale(x, value):
    return clean({k: value*v for k, v in x.items()})


def action(state, site, create):
    occupied = bool(state & (1 << site))
    if occupied == create:
        return None
    return state ^ (1 << site), (-1) ** ((state & ((1 << site)-1)).bit_count())


def sandwich(x, n, which):
    # T: c^dagger X c raises both row and column states; S lowers both.
    create = which == 'T'
    out = defaultdict(F)
    for (row, col), value in x.items():
        for site in range(n):
            a, b = action(row, site, create), action(col, site, create)
            if a is not None and b is not None:
                out[a[0], b[0]] += value * a[1] * b[1]
    return clean(out)


def channel(x, n, which):
    out = sandwich(x, n, which)
    return {k: v / (k[0].bit_count() if which == 'T' else n-k[0].bit_count())
            for k, v in out.items()}


def deficit(x, n, b, which):
    diag = {key: (b-(key[0].bit_count() if which == 'T' else n-key[0].bit_count()))*v
            for key, v in x.items()}
    return add(diag, sandwich(x, n, which))


def restrict(x, m):
    return {k: v for k, v in x.items() if k[0].bit_count() == m}


def path_apply(x, n, word):
    for which in reversed(word):
        x = channel(x, n, which)
    return x


def path_data(n, p, m, b):
    assert 1 <= p <= m <= n
    word = ['T']*p + ['S']*p
    denominators = list(range(m, m-p, -1)) + list(range(n-m+p, n-m, -1))
    factors = [1-F(b, d) for d in denominators]
    return word, denominators, factors


def tangent_expansion(x, n, p, m, b):
    """Exact expansion: path = theta X + first order + two-deficit remainder."""
    word, den, a = path_data(n, p, m, b)
    theta = prod(a)
    first = {}
    remainder = {}
    coefficients = {'T': F(0), 'S': F(0)}
    for j, which in enumerate(word):
        weight = prod((a[h] for h in range(2*p) if h != j), start=F(1)) / den[j]
        coefficients[which] += weight
        first = add(first, restrict(deficit(x, n, b, which), m), weight)
        for i in range(j):
            weight2 = prod((a[h] for h in range(i+1, 2*p) if h != j), start=F(1)) / (den[i]*den[j])
            inner = deficit(deficit(x, n, b, which), n, b, word[i])
            term = restrict(path_apply(inner, n, word[:i]), m)
            remainder = add(remainder, term, weight2)
    rhs = add(add(scale(restrict(x, m), theta), first), remainder)
    return rhs, first, remainder, coefficients


def norm2_squared(x, n, m):
    return sum((v*v for (row, col), v in x.items() if row.bit_count() == m), start=F(0)) / comb(n, m)


def scalar_tangent(n, p, m, v):
    word, den, a = path_data(n, p, m, v*p)
    assert min(den) > v*p
    theta = prod(a)
    zt = p*sum(F(1, d-v*p) for d in den[:p])
    zs = p*sum(F(1, d-v*p) for d in den[p:])
    deviation = max(abs(F(n, 2)-(d-v*p)) for d in den)
    assert deviation < F(n, 2)
    coefficient_error_bound = F(p*p)*deviation/(F(n, 2)*(F(n, 2)-deviation))
    remainder_bound = F(p*comb(2*p, 2), min(den)**2)
    return {'theta': theta, 'zeta_T': zt, 'zeta_S': zs,
            'target_each': F(2*p*p, n),
            'coefficient_error_bound': coefficient_error_bound,
            'remainder_bound_without_C_or_second_deficit': remainder_bound}
