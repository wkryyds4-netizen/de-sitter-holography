"""Exact adjacent-charge coherent Racah transform; no positivity assumption."""
from fractions import Fraction as F
from math import comb
from functools import lru_cache
from cubic_linear_dynamics import crossing, weights


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def charged(N, m, a, t):
    term = F(1)
    result = term
    for k in range(min(a, t)):
        term *= F((-a+k)*(a-N-1+k)*(-t+k)*(t-N+k),
                  (1-m+k)*(m-N+k)*(-N-1+k)*(k+1))
        result += term
    return result


def exchanged(N, a, b):
    m = N // 2
    term = F(1)
    result = term
    for k in range(min(a, b)):
        term *= F((-a+k)*(a-N-1+k)*(-b+k)*(b-N-1+k),
                  (1-m+k)*(-m+k)*(-N-1+k)*(k+1))
        result += term
    return result


def rates(N, a):
    m = N // 2
    return (F((a+1-m)*(a-N-1)**2, 2*(2*a-N-1)),
            F(a*a*(a-m-2), 2*(2*a-N-1)))


@lru_cache(None)
def coherent_table(N, m):
    R = min(m-1, N-m)
    dp, dm = comb(N, m), comb(N, m-1)
    ds = [choose(N, r)**2-choose(N, r-1)**2 for r in range(R+1)]
    es = [choose(N, t)*(choose(N, t+1)-choose(N, t-1)) for t in range(R+1)]
    ell2 = [F(m*(m-r)*(N-m+1-r), N-m+1) for r in range(R+1)]
    T = [[charged(N, m, a, t) for t in range(R+1)] for a in range(R+1)]
    mu = {}
    for a in range(R+1):
        for b in range(R+1):
            mu[a,b] = [F(ds[r], m*m*dp*dm)*ell2[r]*sum(
                es[t]*T[a][t]*T[b][t]*T[r][t] for t in range(R+1))
                for r in range(R+1)]
    return ds, es, ell2, T, mu


def mixed_moments(N, a, b):
    """Return unsigned/signed moments of X=x_(a+b)-x_r under mu."""
    m = N // 2
    x = lambda r: r*(N+1-r)
    ct = F(m*(m-1)*(N+1), N-1)
    xi = exchanged(N, a, b)
    delta = (F(x(b), N+1-2*b)*(charged(N,m,a,b-1)-charged(N,m,a,b))) if b else F(0)
    return [F(1), -2*a*b+F(x(a)*x(b), ct)], [xi, -2*a*b*xi-delta]


def plus_moments(N, a, b):
    m = N // 2
    x = lambda r: r*(N+1-r)
    cf = F(m*m*(N+1), N)
    chi = crossing(N, m, a, b)
    return [F(1), -2*a*b+F(x(a)*x(b), cf)], [chi, -2*a*b*chi]


def cross_BX_A(N, p, a, b):
    """Per ordered input-rank pair, divided by C_a C_b; X has /N."""
    m = N//2
    assert N % 2 == 0 and p < m and a <= p and b <= p
    mu, smu = mixed_moments(N,a,b)
    pi, spi = plus_moments(N,a,b)
    gain = F((m-p)**2,m)
    # Combine exact fractions before evaluating; O(m) cancellation is essential.
    return [gain*(mu[k]+smu[k])+(2*p-m)*(pi[k]+spi[k]) for k in range(2)][0:1] + [
        (gain*(mu[1]+smu[1])+(2*p-m)*(pi[1]+spi[1]))/N]
