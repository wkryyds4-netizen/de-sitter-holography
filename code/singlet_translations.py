"""Kinematic positive-energy extensions of the singlet boost kernel.

These are candidate one-particle representations, not a microscopic bulk map.
Light-cone convention Pplus=P0+P1, Pminus=P0-P1; no sqrt(2).
"""
import mpmath as mp

def kernel(m,mu,a,b,eta,d=2,eps=0):
    """<psi| exp[-(eps+ia)P0+ibP1] exp[-i eta B] |psi>."""
    m,mu,a,b,eta,d,eps=map(mp.mpf,[m,mu,a,b,eta,d,eps])
    alpha=1+mp.exp(eta)+mu*(eps+1j*(a-b))/2
    beta=m*m*(eps+1j*(a+b))/(2*mu)
    pref=2**(2*d)*mp.exp(d*eta)/mp.gamma(2*d)
    if beta==0:return pref*mp.gamma(2*d)/alpha**(2*d)
    return pref*2*(beta/alpha)**d*mp.besselk(2*d,2*mp.sqrt(alpha*beta))

def integral(m,mu,a,b,eta,d=2,eps=0,lo='0.0005',hi=64):
    """Real-p quadrature with explicit absolute bounds on the omitted tails."""
    m,mu,a,b,eta,d,eps,lo,hi=map(mp.mpf,[m,mu,a,b,eta,d,eps,lo,hi])
    alpha=1+mp.exp(eta)+mu*(eps+1j*(a-b))/2
    beta=m*m*(eps+1j*(a+b))/(2*mu)
    pref=2**(2*d)*mp.exp(d*eta)/mp.gamma(2*d)
    edges=[lo];x=lo
    while x<hi:
        x=min(x*mp.mpf('1.5'),hi);edges.append(x)
    v=pref*mp.quad(lambda p:p**(2*d-1)*mp.exp(-alpha*p-beta/p),edges)
    lower=pref*lo**(2*d)/(2*d)
    upper=pref*mp.gammainc(2*d,alpha.real*hi,mp.inf)/alpha.real**(2*d)
    absolute_integral_bound=pref*mp.gamma(2*d)/alpha.real**(2*d)
    return v,lower+upper,absolute_integral_bound

def two_mass_density(s,d=2):
    s=mp.mpf(s);k=2*mp.mpf(d)
    return 2/(mp.beta(k,k)*s**(k+1)*mp.sqrt(1-4/s)) if s>4 else mp.mpf(0)
