"""Four independent real Gaussian smearings of the specified CSTY eikonal contours.

This is a smeared approximation, not a nonperturbative gravity expectation.
SciPy doubles; small-c asymptotics are checked against mpmath separately.
"""
import numpy as np
from scipy.special import exp1,expi
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss

def j(c,side=1):
    c=np.asarray(c,dtype=float);out=np.empty(c.shape,dtype=complex)
    small=np.abs(c)<1/128
    term=np.ones(np.count_nonzero(small));value=term.copy()
    for n in range(1,25):
        term=term*(-c[small])*n;value+=term
    out[small]=value
    pos=(c>=1/128);a=1/c[pos];out[pos]=a*np.exp(a)*exp1(a)
    neg=(c<=-1/128);a=-1/c[neg];out[neg]=a*np.exp(-a)*expi(a)
    neg=c<0;a=-1/c[neg]
    out[neg]-=side*1j*np.pi*a*np.exp(-a)
    return out

def normal(n):
    x,w=hermgauss(n);return np.sqrt(2)*x,w/np.sqrt(np.pi)

def angular(n,gamma):
    if not 0<=gamma<1:raise ValueError('0 <= 4Gm < 1')
    x,w=leggauss(n);edges=[0,np.pi]
    if gamma<.5:edges.insert(1,np.pi-np.pi/(2*(1-gamma)))
    phi=np.concatenate([l+(r-l)*(x+1)/2 for l,r in zip(edges,edges[1:])])
    weight=np.concatenate([(r-l)*w/(2*np.pi) for l,r in zip(edges,edges[1:])])
    H=-np.cos((1-gamma)*(np.pi-phi))/(2*np.sinc(gamma))
    return H,weight

def moment(r,sigma,n=64):
    x,w=normal(n);return float(w@(1/np.cosh(sigma*x/np.sqrt(2)))**r)

def enhancement(n,sigma,order=64):
    return np.exp(n*n*sigma*sigma/2)*(moment(n+1,sigma,order)/moment(1,sigma,order))**2

def smeared(m,gamma,t,sigma,n=24,na=32,boosted=False,center_only=False):
    if m<=0 or sigma<0:raise ValueError('positive m and nonnegative sigma required')
    z,w=normal(n);x=sigma*z/np.sqrt(2)
    C=(np.cosh(x)[:,None]*np.cosh(x)[None,:]).reshape(-1,1)
    ww=(w[:,None]*w[None,:]).reshape(-1,1)*w[None,:]
    if center_only:C=np.ones_like(C)
    T=t+sigma*z[None,:];H,hw=angular(na,gamma)
    total=0j;residue=0.;zero=0j
    for h,v in zip(H,hw):
        cp=h*np.exp(T)/(4*C) if boosted else h*np.exp(T)/(4*m*C)
        cm=np.zeros_like(cp) if boosted else h*np.exp(-T)/(4*m*C)
        val=j(cp,1)+j(cm,-1)-1
        total+=v*np.sum(ww*val/C)
        if h<0:
            ap=-1/cp
            im=-np.pi*ap*np.exp(-ap)
            if not boosted:
                am=-1/cm;im+=np.pi*am*np.exp(-am)
            residue+=v*np.sum(ww*im/C)
    base=float(np.sum(ww/C))
    a=2*np.pi*(1-gamma)
    return dict(value=[float(total.real/base),float(total.imag/base)],
                residue_imag=float(residue/base),baseline=base/(16*a*a),
                A1_squared=base,n=n,angular_order=na,nodes=int(n**3*len(H)),
                boosted=boosted,center_only=center_only)
