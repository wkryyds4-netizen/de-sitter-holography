"""A regulated Rindler midpoint realization of the small-lambda singlet kernel.

All boost variables dimensionless, hbar=c=1, free massive scalar in 1+1D.
KL inverse/Plancherel are standard; no finite-G identification is assumed.
"""
import mpmath as mp
import numpy as np
from scipy.special import roots_legendre


def spectrum(nu):
    nu = mp.mpf(nu)
    if not nu:
        return mp.mpf(8)/3
    return 8*mp.pi/3*nu*(1+nu**2)/mp.sinh(mp.pi*nu)


def target(eta):
    return mp.sech(eta/2)**4


def amplitude(nu, epsilon):
    return mp.sqrt(mp.pi*spectrum(nu))*mp.exp(-epsilon*nu**2/2)


def midpoint(eta, epsilon):
    return mp.quad(lambda n: spectrum(n)*mp.exp(-epsilon*n*n)*mp.cos(n*eta)/mp.pi,
                   [0, 1, 3, 6, 12, 24, mp.inf])


def radial_norm(epsilon):
    e = mp.mpf(epsilon)
    return 4*mp.sqrt(mp.pi)/(3*e**mp.mpf('1.5'))+2*mp.sqrt(mp.pi)/e**mp.mpf('2.5')


def radial_profiles(log_x, epsilons, order=128, cutoff=24, dps=25):
    """Evaluate KL inverse for x=m*rho, sharing K evaluations across regulators.

    Gauss-Legendre nodes and final dot products are float64. Bessel values are
    computed with mpmath at dps; this is not an interval quadrature certificate.
    """
    nodes, weights = roots_legendre(order)
    nodes = (nodes+1)*cutoff/2
    weights = weights*cutoff/2
    with mp.workdps(dps):
        multipliers = np.array([[float(2/mp.pi**2*n*mp.sinh(mp.pi*n)*amplitude(n,e))
                                 for n in nodes] for e in epsilons])
        values=[]
        conditions=[]
        for s in log_x:
            x=mp.exp(float(s))
            kvals=np.array([float(mp.re(mp.besselk(1j*mp.mpf(float(n)),x))) for n in nodes])
            terms=multipliers*(weights*kvals)
            val=terms.sum(axis=1)
            values.append(val)
            conditions.append(np.sum(np.abs(terms),axis=1)/np.maximum(np.abs(val),1e-300))
    return np.array(values),np.array(conditions)


def radial_tail_bound(x, epsilon, cutoff):
    """Conservative exact-arithmetic pointwise frequency tail using |K_iν|<=K_0.

    For ν>=L>=1, sqrt(ν(1+ν²))<=sqrt(2)ν^(3/2)<=sqrt(2)ν².
    sqrt(sinh(pi ν))<=exp(pi ν/2)/sqrt(2). Complete square,
    bound exp(pi ν/2-eν²/2)<=exp(-eν²/4) if L>=2pi/e.
    This deliberately conservative optional bound is not used for main plots.
    """
    e,L=mp.mpf(epsilon),mp.mpf(cutoff)
    if L < max(1, 2*mp.pi/e):
        raise ValueError('cutoff too small for this simple tail bound')
    # inverse integrand <= 4 sqrt(2/3)/pi * K0(x) * nu^3 exp(-e nu²/4)
    return 4*mp.sqrt(mp.mpf(2)/3)/mp.pi*mp.besselk(0,x)*mp.gammainc(2,e*L*L/4,mp.inf)/(2*(e/4)**2)
