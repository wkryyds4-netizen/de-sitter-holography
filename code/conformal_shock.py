"""Two conformal probes in the finite-G dS3 eikonal formula.

Source: arXiv:2607.14042v1, G.2, 2.25-2.26, 2.33-2.35.
Neither contour is designated as a quantum-state expectation value.
"""
from __future__ import annotations
import mpmath as mp

SOURCE = 'https://arxiv.org/html/2607.14042v1'


def alpha(deficit):
    deficit=mp.mpf(deficit)
    if not 0 <= deficit < 1:
        raise ValueError('Require 0 <= 4 G m < 1.')
    return 2*mp.pi*(1-deficit)


def shock_profile(phi, m, deficit):
    m,deficit,phi=map(mp.mpf,(m,deficit,phi))
    if m <= 0:
        raise ValueError('Observer mass must be positive.')
    alpha(deficit)
    sinc=mp.sin(mp.pi*deficit)/(mp.pi*deficit) if deficit else mp.mpf(1)
    return -mp.cos((1-deficit)*(mp.pi-abs(phi)))/(2*m*sinc)


def angular_points(deficit):
    k=1-mp.mpf(deficit)
    if k > mp.mpf('0.5'):
        return [mp.mpf(0),mp.pi-mp.pi/(2*k),mp.pi]
    return [mp.mpf(0),mp.pi]


def angular_average(func, deficit):
    return mp.quad(func,angular_points(deficit))/mp.pi


def shock_moments(m, deficit):
    m,deficit=map(mp.mpf,(m,deficit))
    a=alpha(deficit)
    G=deficit/(4*m)
    sinc=mp.sin(mp.pi*deficit)/(mp.pi*deficit) if deficit else mp.mpf(1)
    return -4*mp.pi*G/a, (1+mp.sin(a)/a)/(8*m*m*sinc*sinc)


def j_boundary(c, side):
    """J(c+i*side*0); explicit Ei formula on the negative real cut."""
    c=mp.mpf(c)
    if side not in (-1,1):
        raise ValueError('side must be +1 or -1.')
    if not c:
        return mp.mpf(1)
    if c > 0:
        a=1/c
        return a*mp.exp(a)*mp.e1(a)
    a=-1/c
    return a*mp.exp(-a)*mp.ei(a)-side*1j*mp.pi*a*mp.exp(-a)


def j_complex(c):
    c=mp.mpc(c)
    if not c:
        return mp.mpf(1)
    if not c.imag and c.real < 0:
        raise ValueError('Select a side explicitly with j_boundary for negative real c.')
    return mp.exp(1/c)*mp.e1(1/c)/c


def j_integral(c):
    """Independent convergent rational integral away from the negative-real cut."""
    c=mp.mpc(c)
    if not c.imag and c.real < 0:
        raise ValueError('Use a principal-value calculation on the cut.')
    candidate=mp.re(-1/c) if c else 0
    points=[0,1,mp.inf]
    if candidate > 0:
        points=sorted(set([mp.mpf(0),mp.mpf(1),candidate/2,candidate,2*candidate]))+[mp.inf]
    return mp.quad(lambda x:mp.exp(-x)/(1+c*x),points)


def j_pv_real(c):
    """Real principal value for c<0, by symmetric pole subtraction."""
    c=mp.mpf(c)
    if c >= 0:
        raise ValueError('PV formula requires c<0.')
    a=-1/c
    def subtracted(x):
        y=x-a
        return -mp.exp(-a) if not y else mp.exp(-a)*mp.expm1(-y)/y
    core=mp.quad(subtracted,[0,a/2,a,3*a/2,2*a])
    tail=mp.quad(lambda x:mp.exp(-x)/(x-a),[2*a,2*a+1,mp.inf])
    return -a*(core+tail)


def disconnected(m, deficit):
    """Same-geometry disconnected/TOC baseline at retained eikonal order."""
    return 1/(16*alpha(deficit)**2)


def conformal_ratio(m, deficit, t, contour=1, epsilon=None):
    """R12 (contour=+1) or R14 (-1) at the symmetric configuration.

    epsilon=None evaluates the exact lateral boundary within the eikonal
    model. epsilon>0 evaluates its explicitly regulated Euclidean approach.
    The total is the two original sectors minus the duplicated disconnected 1.
    """
    m,deficit,t=map(mp.mpf,(m,deficit,t))
    if contour not in (-1,1):
        raise ValueError('contour must be +1 (12) or -1 (14).')
    et,emt=mp.exp(t),mp.exp(-t)
    if epsilon is None:
        def f(phi):
            h=shock_profile(phi,m,deficit)
            return j_boundary(h*et/4,contour)+j_boundary(h*emt/4,-contour)-1
    else:
        epsilon=mp.mpf(epsilon)
        if not 0 < epsilon < mp.pi/2:
            raise ValueError('Require 0 < epsilon < pi/2.')
        ep=mp.exp(-1j*contour*epsilon)
        def f(phi):
            h=shock_profile(phi,m,deficit)
            return j_complex(h*et*ep/4)+j_complex(h*emt/ep/4)-1
    return angular_average(f,deficit)


def imaginary_ratio_integral(m, deficit, t):
    """Independent real angular integral for Im R12, from pole residues."""
    m,deficit,t=map(mp.mpf,(m,deficit,t))
    et,emt=mp.exp(t),mp.exp(-t)
    def f(phi):
        h=shock_profile(phi,m,deficit)
        if h >= 0:
            return mp.mpf(0)
        a=-h
        return 4*mp.pi/a*(et*mp.exp(-4*et/a)-emt*mp.exp(-4*emt/a))
    return angular_average(f,deficit)


def perturbative_ratio(m, deficit, t):
    h1,h2=shock_moments(m,deficit)
    return 1-h1*(mp.exp(t)+mp.exp(-t))/4+h2*(mp.exp(2*t)+mp.exp(-2*t))/8
