"""Energy pinches and the Gamma-regularized pure-recoil topology graph.

No physical insertion prescription or finite-G correlator is assumed.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, platform, sys
import mpmath as mp
import sympy as sp
from recoil_internal import parameters, residue, vertical, VERT, QUAD, sb_scaled
from recoil_vertex import triples, logT

ROOT = Path(__file__).resolve().parents[1]

def setup(M, nu, omega, u=0, v=0):
    ts, m = triples(M, nu, omega, u, v)
    edges = [m[0], nu, m[1], m[2], nu, m[3]]
    a, r = parameters(edges)
    return a, r, edges, ts

def ps(a, r):
    return [[mp.mpf('.5')+1j*(x-y) for y in r] for x in a]

def gamma_product(a, r):
    return mp.fprod(mp.gamma(x) for row in ps(a, r) for x in row)

def regular_left(a, r):
    """Four left residue families, divided by all 12 pinch Gamma factors."""
    p = ps(a, r); terms = []
    for j, x in enumerate(a):
        others = [k for k in range(4) if k != j]
        coeff = (mp.exp(-mp.pi*x)*mp.rgamma(mp.mpf('.5')-1j*x)
                 *mp.fprod(mp.gamma(1j*(a[k]-x)) for k in others)
                 *mp.fprod(mp.rgamma(p[k][l]) for k in others for l in range(3)))
        if coeff == 0:
            terms.append(mp.mpc(0)); continue
        h = mp.hyper([mp.mpf('.5')+1j*x]+p[j],
                     [1+1j*(x-a[k]) for k in others], 1)
        terms.append(coeff*h)
    value = mp.fsum(terms)
    return value, mp.fsum(abs(z) for z in terms)/abs(value), terms

def regular_right(a, r):
    """Three right residue families; works at coincident left a parameters."""
    p = ps(a, r); terms = []
    for l, y in enumerate(r):
        others = [h for h in range(3) if h != l]
        coeff = (-1j*mp.exp(-mp.pi*y)*mp.rgamma(1-1j*y)
                 *mp.fprod(mp.gamma(1j*(y-r[h])) for h in others)
                 *mp.fprod(mp.rgamma(p[j][h]) for j in range(4) for h in others))
        if coeff == 0:
            terms.append(mp.mpc(0)); continue
        h = mp.hyper([p[j][l] for j in range(4)],
                     [1-1j*y]+[1+1j*(r[h]-y) for h in others], 1)
        terms.append(coeff*h)
    value = mp.fsum(terms)
    return value, mp.fsum(abs(z) for z in terms)/abs(value), terms

def leg_product(edges):
    return mp.fprod(x/mp.sinh(mp.pi*x) if x else 1/mp.pi for x in edges)

def regular_graph(M, nu, omega, u=0, v=0, method='right'):
    a, r, edges, _ = setup(M, nu, omega, u, v)
    fun = regular_right if method == 'right' else regular_left
    h, c1, _ = fun(a, r)
    ac, rc, _, _ = setup(M, nu, mp.conj(omega), u, v)
    hs0, c2, _ = fun(ac, rc); hs = mp.conj(hs0)
    return mp.pi**6*leg_product(edges)*h*hs, h, hs, max(c1, c2)

def simple_residue(a, r, j, l, N, sign):
    """Coefficient of (omega-omega0)^(-1); dp/domega=-i sign/2."""
    p = ps(a, r); x = a[j]
    others = [k for k in range(4) if k != j]
    pp = p[j].copy(); pp[l] = -N
    poly = mp.hyper([mp.mpf('.5')+1j*x]+pp,
                    [1+1j*(x-a[k]) for k in others], 1)
    return (2j*sign*(-1)**N/mp.factorial(N)*mp.exp(-mp.pi*x)
            *mp.rgamma(mp.mpf('.5')-1j*x)
            *mp.fprod(mp.gamma(1j*(a[k]-x)) for k in others)
            *mp.fprod(mp.gamma(p[j][h]) for h in range(3) if h != l)*poly)

def direct_c(ms):
    """Mass-dependent chiral C0 recoil factor, up to a common constant."""
    S = sum(ms)
    return (mp.exp(mp.pi*S/2)*mp.gamma((1-1j*S)/2)
            *mp.fprod(mp.gamma((1-1j*(S-2*x))/2) for x in ms)
            *mp.fprod(mp.rgamma(1-1j*x) for x in ms))

def direct_c_pair(ms):
    """Schwarz continuation of |C0|^2 in the recoil limit."""
    S = sum(ms)
    return (mp.pi*mp.exp(mp.pi*S)*mp.fprod(mp.sinh(mp.pi*x)/x for x in ms)
            /mp.cosh(mp.pi*S/2)
            /mp.fprod(mp.cosh(mp.pi*(S-2*x)/2) for x in ms))

def v2_limit(ms):
    S = sum(ms)
    return mp.exp(mp.pi*S/2)*mp.gamma((1-1j*S)/2)*mp.fprod(
        mp.rgamma((1-1j*(S-2*x))/2) for x in ms)

def v2_finite(ms, kappa):
    """Exact double-sine V^2, using inverse-b shifts before evaluation."""
    S = sum(ms); L = -1j*kappa
    zz = [(1-1j*S)/2]+[(1+1j*(S-2*x))/2 for x in ms]
    return (mp.fprod(sb_scaled(z, kappa) for z in zz)
            *mp.fprod(2*mp.sin(mp.pi*z) for z in zz[1:])
            /(2*mp.sin(mp.pi*(zz[0]-L))))

if __name__ == '__main__':
    raise SystemExit('Use the separate validation driver; this module does not overwrite frozen outputs.')
