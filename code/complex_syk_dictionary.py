"""Explicit complex-SYK conventions and finite-dimensional operator construction."""
from itertools import combinations
from math import comb, factorial
from fractions import Fraction
import numpy as np


def fermions(n):
    d = 1 << n
    result = []
    for i in range(n):
        c = np.zeros((d, d), dtype=complex)
        for state in range(d):
            if state & (1 << i):
                sign = (-1) ** ((state & ((1 << i) - 1)).bit_count())
                c[state ^ (1 << i), state] = sign
        result.append(c)
    return result


def creators(cs, p):
    out = []
    for inds in combinations(range(len(cs)), p):
        a = np.eye(len(cs[0]), dtype=complex)
        for i in inds:
            a = a @ cs[i].conj().T
        out.append(a)
    return out


def covariance(n, p, js=1.0):
    """Variance in CAR=2 units, for E[J_IJ J_KL]=v delta_IL delta_JK."""
    return js**2 * float(Fraction(n, p*p) * Fraction(factorial(p), n**p)**2)


def exact_moments(n, p):
    """Rational ensemble H^2 and A2 means in standard CAR, Js=1."""
    counts = [(k, comb(p,k)*comb(n-p,p-k)*2**k)
              for k in range(max(0,2*p-n), p+1)]
    v = Fraction(n, p*p) * Fraction(factorial(p), n**p)**2
    w2 = v * comb(n,p) * sum(w for _,w in counts)
    a2 = v * comb(n,p) * sum((2*p-k)*w for k,w in counts)
    return w2, a2


def sample_hamiltonian(n, p, seed, js=1.0):
    cs = fermions(n)
    basis = creators(cs,p)
    rng = np.random.default_rng(seed)
    m = len(basis)
    raw = rng.normal(size=(m,m)) + 1j*rng.normal(size=(m,m))
    j = np.sqrt(covariance(n,p,js)) * (raw + raw.conj().T)/2
    h = np.zeros_like(cs[0])
    for a in range(m):
        for b in range(m):
            h += (2**p * j[a,b]) * (basis[a] @ basis[b].conj().T)
    return cs, j, h


def comm(a,b):
    return a@b-b@a


def a2_operator(h, cs):
    return sum(comm(h,c).conj().T @ comm(h,c) for c in cs)


def relative(a,b):
    return float(np.linalg.norm(a-b)/max(1.0,np.linalg.norm(a),np.linalg.norm(b)))


def charge_spectrum(h, op, n):
    """Resolve charge exactly before diagonalization; never mix degenerate blocks."""
    d = 1 << n
    charges = np.array([s.bit_count() for s in range(d)])
    omega, weight, charge = [], [], []
    min_gap = float('inf')
    for k in range(n+1):
        inds = np.flatnonzero(charges == k)
        e,u = np.linalg.eigh(h[np.ix_(inds,inds)])
        o = u.conj().T @ op[np.ix_(inds,inds)] @ u
        gaps = e[None,:]-e[:,None]
        off = np.abs(gaps[np.triu_indices(len(e),1)])
        nonzero = off[off > 1e-10]
        if len(nonzero):
            min_gap = min(min_gap,float(nonzero.min()))
        omega.extend(gaps.ravel())
        weight.extend((np.abs(o)**2/d).ravel())
        charge.extend([k]*(len(e)**2))
    return np.array(omega),np.array(weight),np.array(charge),min_gap
