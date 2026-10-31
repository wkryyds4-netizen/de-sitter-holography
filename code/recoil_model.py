"""Pure-recoil heavy-light dS3 baseline from arXiv:2607.14042v1 (2.78).

This module reproduces a specified semiclassical correlator. It does not
implement finite-clock preparation or select a physical contour past a cut.
"""
from __future__ import annotations

import mpmath as mp

SOURCE_URL = 'https://arxiv.org/html/2607.14042v1'
SOURCE_EQUATIONS = ['2.76', '2.77', '2.78', '2.82', '2.83']


def _nu(value):
    return mp.mpc(value)


def _real_q(q):
    q = mp.mpf(q)
    if not 0 <= q < 1:
        raise ValueError('This representation requires real 0 <= q < 1.')
    return q


def recoil_hyper(nu, q):
    """Analytic continuation of 3F2, including an explicitly complex q."""
    nu = _nu(nu)
    return mp.hyper([1, 1 - 1j * nu, 1 + 1j * nu],
                    [mp.mpf('1.5'), mp.mpf('1.5')], q * q)


def recoil_series(nu, q, tol=None, max_terms=100000):
    """Positive-series sum, with a certified geometric remainder bound.

    Supported nu is real or purely imaginary with |Im nu|<1. Once
    n >= ceil(|nu|^2)+2, a_(n+1)/a_n <= q^2. `tail_bound` bounds the
    unsummed infinite tail, apart from numerical roundoff.
    """
    nu = _nu(nu)
    q = _real_q(q)
    if nu.imag and (nu.real or abs(nu.imag) >= 1):
        raise ValueError('Tail proof supports real nu or nu=i*mu, |mu|<1.')
    tol = mp.mpf(tol) if tol is not None else mp.power(10, -mp.mp.dps + 12)
    z = q * q
    total = term = mp.mpf(1)
    if not z:
        return {'value': total, 'terms': 1, 'tail_bound': mp.mpf(0)}
    safe_n = int(mp.ceil(abs(nu) ** 2)) + 2
    nu2 = mp.re(nu * nu)
    for n in range(max_terms):
        next_term = term * z * ((n + 1) ** 2 + nu2) / (mp.mpf(n) + mp.mpf('1.5')) ** 2
        if n >= safe_n:
            bound = next_term / (1 - z)
            if bound <= tol:
                return {'value': total, 'terms': n + 1, 'tail_bound': bound}
        total += next_term
        term = next_term
    raise RuntimeError('Series did not meet its requested tail bound.')


def geometry_integrand(nu, x):
    """Normalized propagator evaluated geometrically, without hypergeometric calls."""
    nu = _nu(nu)
    if not x:
        return mp.mpf(1)
    root = mp.sqrt(1 - x*x)
    if not nu:
        return mp.asin(x) / (x * root)
    return mp.sinh(2 * nu * mp.asin(x)) / (2 * nu * x * root)


def recoil_angle(nu, q):
    """Independent ordinary integral over u=sin(theta), for q<1."""
    q = _real_q(q)
    return mp.quad(lambda u: geometry_integrand(nu, q * mp.sqrt(1 - u*u)),
                   [0, mp.mpf('0.5'), 1])


def recoil_conformal(q):
    """nu=i/2 closed form; complex q fixes the side of the logarithm."""
    return mp.atanh(q) / q if q else mp.mpf(1)


def conformal_boundary(q, side):
    """Boundary value at real q>1, with side=+1/-1 in the q plane."""
    q = mp.mpf(q)
    if q <= 1 or side not in (-1, 1):
        raise ValueError('Require q>1 and side in {-1,1}.')
    return mp.log((q + 1)/(q - 1))/(2*q) + side*1j*mp.pi/(2*q)


def q_from_time(m, nu_heavy, t):
    return mp.mpf(nu_heavy) * mp.cosh(t) / (2 * mp.mpf(m))


def time_branch_q(q_real, ratio, time_sign, contour, epsilon):
    """Return t and q(t+i*contour*epsilon).

    contour=+1 labels the F12 continuation in (2.82); contour=-1 labels
    F14 in (2.83). At negative real t, the side in the q/z plane flips.
    No physical prescription is selected by this helper.
    """
    if time_sign not in (-1, 1) or contour not in (-1, 1):
        raise ValueError('Signs must be +1 or -1.')
    q_real, ratio, epsilon = map(mp.mpf, (q_real, ratio, epsilon))
    if ratio <= 0 or 2*q_real/ratio < 1 or epsilon <= 0:
        raise ValueError('Invalid branch parameters.')
    t = time_sign * mp.acosh(2*q_real/ratio)
    return t, ratio * mp.cosh(t + contour*1j*epsilon)/2
