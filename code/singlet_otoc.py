"""BINT (3.9),(4.9),(4.10), with symmetric Jacobi normalization.

The selected crossed channel excludes its mutual-probe weight q^4.
All units Js=1. This module does no filesystem writes on import.
"""
from math import exp,log,lgamma,sqrt
import numpy as np


def taylor_tail(a,order):
    if not a:return 0.
    if order<=a:return float('inf')
    return exp(order*log(a)-lgamma(order+1)-log(1-a/(order+1)))


def cutoff(lam,tmax,tol=1e-10):
    q=exp(-4*lam);r=q*q;Q=q**4
    Cq=1/((1-q)*(1-q*q)*(1-q**3));K=Cq/(1-r)
    a=2*tmax/sqrt(lam*(1-q))
    inner=max(1,int(log(tol/(100*(1+K)))/log(r))+1)
    L=inner+max(12,int(a)+2)
    while True:
        epsvac=2*taylor_tail(a,L)
        epsprop=2*taylor_tail(a,L-inner)
        delta=epsvac+epsprop+2*r**inner
        bound=(4+2*Q*K)*delta
        if bound<tol:break
        L+=1
    return dict(L=L,inner=inner,K=K,RH=2/sqrt(lam*(1-q)),
        vacuum_vector_bound=epsvac,low_source_propagator_bound=epsprop,
        source_tail_bound=r**inner,total_otoc_truncation_bound=bound,
        crossed_reduced_truncation_bound=2*K*delta,uncrossed_truncation_bound=2*delta)


class Kernel:
    def __init__(self,lam,L):
        self.lam=lam;self.L=L;self.q=exp(-4*lam)
        self.r=self.q**2;self.Q=self.q**4
        n=np.arange(1,L)
        off=np.sqrt(-np.expm1(-4*lam*n)/(-np.expm1(-4*lam)*lam))
        self.e,self.u=np.linalg.eigh(np.diag(off,1)+np.diag(off,-1))
        self.rn=np.exp(-8*lam*np.arange(L))
        # log (q;q)_n avoids factorial similarity factors.
        self.logf=np.r_[0.,np.cumsum(np.log(-np.expm1(-4*lam*n)))]
        i=np.arange(L)
        logci=sum(np.log(-np.expm1(-4*lam*(i+j))) for j in [1,2,3])-sum(log(-np.expm1(-4*lam*j)) for j in [1,2,3])-8*lam*i
        l=i[:,None];ii=i[None,:];index=l+ii;mask=index<L
        safe=np.minimum(index,L-1)
        self.indices=safe
        self.weight=np.where(mask,np.exp(.5*logci[None,:]-8*lam*l+.5*(self.logf[safe]-self.logf[l])),0.)

    def at(self,t):
        minus=(self.u*np.exp(-1j*t*self.e)[None,:])@self.u.T
        plus=minus.conj();vac=minus[:,0]
        v=plus@(self.rn*vac)
        uncrossed=np.dot(v,v)
        # Positive sum for the selected crossed channel, including every i.
        V=self.weight*vac[self.indices]
        evolved=plus@V
        crossed=float(np.sum(self.rn[:,None]*abs(evolved)**2))
        g2=float(np.dot(self.rn,abs(vac)**2))
        g4=float(np.dot(self.rn**2,abs(vac)**2))
        otoc=2*uncrossed.real+self.Q*crossed
        f22=1+(1+self.Q)*g4
        return dict(time=t,uncrossed=float(uncrossed.real),uncrossed_imag=float(uncrossed.imag),
            crossed_reduced=crossed,crossed_weighted=self.Q*crossed,
            g2=g2,g4=g4,otoc=otoc,f22=f22,commutator_squared=2*(f22-otoc),
            connected_otoc=otoc-1-2*g2*g2,
            vacuum_norm_residual=float(abs(np.vdot(vac,vac).real-1)))
