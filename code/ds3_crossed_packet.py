"""Actual regulated two-creation overlap for a compact initial observer packet.

No internal spectral filter is part of the model. L is a numerical cutoff;
the omitted physical state is controlled by positive Hilbert-norm tails.
"""
import numpy as np
from scipy.special import roots_legendre


def probability(c,c0=.4,c1=.8):
    x=(c-c0)/(c1-c0)
    return np.where((x>=0)&(x<=1),630*x**4*(1-x)**4/(c1-c0),0)


def g(x):
    return x*np.exp(np.pi/(np.sqrt(1+x*x)+x))/(1+np.exp(-2*np.pi*x))


def balanced(a,b,p):
    A=np.cosh(np.pi*(a+b));B=np.cosh(np.pi*(a-b));C=np.cosh(np.pi*p)
    return (A-B)/(8*np.pi*(A+C)*(B+C))


def mesh(order,s=.2,L=8):
    breaks=[s]+[x for x in [1,2,4,6,8,10,12,16,20] if s<x<L]+[L]
    xs,ws=roots_legendre(order)
    return (np.concatenate([l+(r-l)*(xs+1)/2 for l,r in zip(breaks,breaks[1:])]),
            np.concatenate([ws*(r-l)/2 for l,r in zip(breaks,breaks[1:])]))


def initial_mesh(order,c0=.4,c1=.8):
    x,w=roots_legendre(order);c=c0+(c1-c0)*(x+1)/2
    return c,w*(c1-c0)/2*probability(c,c0,c1)


def packet_weights(a,b,d,c,p=.6,q=.9):
    return g(a)/g(c)*np.sqrt(balanced(a,b,q)*balanced(b,c,p)*
                            balanced(a,d,p)*balanced(d,c,q))


def norm_tail_bound(L,s=.2,c1=.8):
    """Analytic positive norm tail when a>L or internal label>L.

    Uniform over the normalized initial probability on [c0,c1].
    The full double-difference plane overcounts physical a,b>=s.
    """
    R=L-c1
    if R<=0:raise ValueError('L must exceed c1')
    pi=np.pi;B=2*np.exp(pi/(np.sqrt(1+s*s)+s))
    return B*np.exp(-pi*R)/(4*pi*pi)*(
        R/pi+1/pi**2+(R*R/pi+2*R/pi**2+2/pi**3)/s)


def positive_norms(order=24,nc=12,L=16,s=.2,p=.6,q=.9):
    x,w=mesh(order,s,L);cs,cw=initial_mesh(nc)
    a=x[:,None];b=x[None,:];ww=w[:,None]*w[None,:]
    ns=nt=0.
    for c,wc in zip(cs,cw):
        ns+=wc*np.sum(ww*g(a)/g(c)*balanced(a,b,q)*balanced(b,c,p))
        nt+=wc*np.sum(ww*g(a)/g(c)*balanced(a,b,p)*balanced(b,c,q))
    return dict(ns=float(ns),nt=float(nt),far_norm_bound=float(norm_tail_bound(L,s)),
                order=order,nc=nc,L=L,initial_norm=float(sum(cw)))
