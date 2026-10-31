"""Specified CSTY lateral eikonal model, real principal probes, radius one.

Rapidity integrals keep signed cosine weights. No Hilbert positivity claim.
Sobol errors are empirical across independent scrambles, not certified bounds.
"""
import numpy as np
import mpmath as mp
from scipy.special import exp1, expi, ndtri
from scipy.stats import qmc


def kernel(c, side=1):
    """-J'(c); side=+1 has negative J imaginary part on c<0."""
    c=np.asarray(c,dtype=float)
    out=np.empty(c.shape,dtype=complex)
    small=np.abs(c)<1/128
    term=np.ones(np.count_nonzero(small)); value=term.copy()
    for n in range(1,25):
        term *= -c[small]*(n+1)**2/n
        value += term
    out[small]=value
    pos=c>=1/128
    a=1/c[pos]; j=a*np.exp(a)*exp1(a)
    out[pos]=((1+c[pos])*j-1)/c[pos]**2
    neg=c<=-1/128
    a=-1/c[neg]; j=a*np.exp(-a)*expi(a)
    out[neg]=((1+c[neg])*j-1)/c[neg]**2
    neg=c<0; a=-1/c[neg]
    # Polynomial times exp: retain the residue even in the small-c branch.
    im=np.zeros_like(a); keep=a<740
    im[keep]=np.pi*a[keep]**2*(a[keep]-1)*np.exp(-a[keep])
    out[neg]-=side*1j*im
    return out


def spectral_moment(n,nu,sigma,dps=50):
    with mp.workdps(dps):
        nu=mp.mpf(str(nu)); sigma=mp.mpf(str(sigma)); b=mp.mpf(n+2)/4
        def density(w):
            M=2**(mp.mpf(n)/2-1-1j*w)*mp.gamma(b+1j*(nu-w)/2)*mp.gamma(b-1j*(nu+w)/2)
            return mp.exp(-sigma*sigma*w*w)*abs(M)**2/mp.pi
        return mp.quad(density,[0,1,4,12,mp.inf])


def sobol_integral(nu,nup,sigma,t,m,gamma,power,seed,base,rate=1.):
    """Full 8D expectation, plus an exactly centered baseline control variate.

    Fixed denominator is separately integrated B0(nu)*B0(nup).
    Both raw and control-variate estimates are retained.
    """
    if m<=0 or not 0<=gamma<1 or sigma<0:raise ValueError('invalid parameters')
    sampler=qmc.Sobol(d=8,scramble=True,seed=seed)
    sums=np.zeros(7,dtype=complex)
    # Blocks avoid holding millions of intermediate complex values in memory.
    blockpower=min(power,15)
    for block in range(2**(power-blockpower)):
        u=sampler.random(2**blockpower)
        u=np.clip(u,np.finfo(float).eps,1-np.finfo(float).eps)
        X,Y,Z=(sigma*ndtri(u[:,:3])).T
        X/=np.sqrt(2);Y/=np.sqrt(2)
        r=-np.log(u[:,3:7])/rate
        A=np.exp(Y)*np.cosh(r[:,0])+np.exp(-Y)*np.cosh(r[:,1])
        B=np.exp(X)*np.cosh(r[:,2])+np.exp(-X)*np.cosh(r[:,3])
        weight=np.exp(rate*r.sum(axis=1))/(rate**4*A*A*B*B)
        weight*=np.cos(nu*r[:,0])*np.cos(nu*r[:,1])*np.cos(nup*r[:,2])*np.cos(nup*r[:,3])
        phi=np.pi*u[:,7]
        h=-np.cos((1-gamma)*(np.pi-phi))/(2*m*np.sinc(gamma))
        T=t+Z
        delta=kernel(h*np.exp(T)/(A*B),1)+kernel(h*np.exp(-T)/(A*B),-1)-2
        correction=weight*delta
        sums+=np.array([weight.sum(),correction.sum(),np.abs(weight).sum(),
                        np.abs(correction.real).sum(),np.abs(correction.imag).sum(),
                        (weight+correction).sum(),np.abs(weight+correction).sum()])
    means=sums/2**power
    cv=1+means[1]/base; raw=means[5]/base
    return dict(power=power,seed=seed,rate=rate,value=[float(cv.real),float(cv.imag)],
                raw_value=[float(raw.real),float(raw.imag)],baseline_ratio=float(means[0].real/base),
                baseline_absolute_ratio=float(means[2].real/base),
                absolute_correction=[float(means[3].real/base),float(means[4].real/base)],
                absolute_raw=float(means[6].real/base))
