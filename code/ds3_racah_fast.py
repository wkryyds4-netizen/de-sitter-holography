"""Batched complex128 LPRS scalar 6j with a balanced-4F3 zeta tail.

Empirical accelerator, not interval arithmetic. The known Gamma-ratio
asymptotics (DLMF 5.11.8) give integer-power tails because excess=1.
Return cancellation diagnostics; compare independent precision/refinements.
"""
import math
import numpy as np
from scipy.special import loggamma, zeta, bernoulli


def gamma_ratio(top, bottom):
    return np.exp(sum(loggamma(x) for x in top)-sum(loggamma(x) for x in bottom))


def balanced_4f3(top,bottom,N=128,K=12):
    """Parameters are broadcast arrays; all 4F3 have exact excess one."""
    top=np.broadcast_arrays(*[np.asarray(x,dtype=complex) for x in top])
    bottom=np.broadcast_arrays(*[np.asarray(x,dtype=complex) for x in bottom],top[0])[:3]
    if np.max(abs(sum(bottom)-sum(top)-1))>1e-10:
        raise ValueError('Requires excess exactly one')
    n=np.arange(N-1,dtype=float)
    ratio=np.ones(top[0].shape+(N-1,),complex)
    for a in top:ratio*=a[...,None]+n
    for b in bottom:ratio/=b[...,None]+n
    ratio/=n+1
    terms=np.concatenate([np.ones(top[0].shape+(1,),complex),np.cumprod(ratio,axis=-1)],axis=-1)
    # Terms n=0,...,N-1; tail begins at N.
    B=bernoulli(K+1)
    def bp(j,x):
        value=np.zeros_like(x)
        for r in range(j+1):value+=math.comb(j,r)*B[r]*x**(j-r)
        return value
    logs=[np.zeros_like(top[0])]
    for k in range(1,K+1):
        diff=sum(bp(k+1,a) for a in top)-sum(bp(k+1,b) for b in bottom)-bp(k+1,np.ones_like(top[0]))
        logs.append((-1)**(k+1)*diff/(k*(k+1)))
    coeff=[np.ones_like(top[0])]
    for k in range(1,K+1):
        coeff.append(sum(j*logs[j]*coeff[k-j] for j in range(1,k+1))/k)
    constant=gamma_ratio(bottom,top)
    tail=constant*sum(c*zeta(2+k,float(N)) for k,c in enumerate(coeff))
    value=terms.sum(axis=-1)+tail
    condition=(np.abs(terms).sum(axis=-1)+abs(tail))/abs(value)
    return value,condition


def omega(h,hp,hs,N=128,K=12):
    h1,h2,h3,h4=hs;p=h2+h3
    h12,h34,h23,h14,h13,h42=h1-h2,h3-h4,h2-h3,h1-h4,h1-h3,h4-h2
    va,ca=balanced_4f3([hp+h23,hp-h14,hp-h12+h34-p+1,hp-p+1],
                       [2*hp,hp-h12+h-p+1,hp-h12-h-p+2],N,K)
    a=gamma_ratio([2*h,hp-p+1,hp-h12+h34-p+1,-hp+h12+h+p-1],
                  [h12+h,h34+h,hp-h12+h-p+1])*va
    vb,cb=balanced_4f3([h13+h+p-1,h42+h+p-1,h34+h,h12+h],
                       [hp+h12+h+p-1,2*h,-hp+h12+h+p],N,K)
    b=gamma_ratio([2*hp,hp-h12-h-p+1,h13+h+p-1,h42+h+p-1],
                  [hp+h23,hp-h14,hp+h12+h+p-1])*vb
    value=a+b
    return value,np.maximum(ca,cb)*(abs(a)+abs(b))/abs(value)


def shadow(h,hi,hj):
    d=hi-hj
    return np.pi*gamma_ratio([2*h-1,1-h-d,1-h+d],[2-2*h,h-d,h+d])


def raw_6j(b,d,external,N=128,K=12):
    b,d,*e=np.broadcast_arrays(np.asarray(b),np.asarray(d),*[np.asarray(x) for x in external])
    hs=[(1+1j*x)/2 for x in e];h,hp=(1+1j*b)/2,(1+1j*d)/2
    h1,h2,h3,h4=hs
    pref=gamma_ratio([h+h1-h2,h-h1+h2,1-h+h4-h3,1-h+h3-h4],
                     [2*h,2-2*h])/4
    terms=[];conditions=[]
    for v,co in [(hp,shadow(1-hp,h1,h4)),(1-hp,shadow(hp,h3,h2))]:
        w,c1=omega(h,v,hs,N,K);w2,c2=omega(1-h,v,hs,N,K)
        terms.append(co*pref*np.sin(np.pi*(v-h1-h4))*np.sin(np.pi*(v-h2-h3))*w*w2)
        conditions.append(np.maximum(c1,c2))
    value=sum(terms)
    condition=np.maximum(*conditions)*sum(abs(x) for x in terms)/abs(value)
    return value,condition


def period_phase(x,y,z):
    S=x+y+z;alphas=[-.5-.5j*(S-2*n) for n in [x,y,z]]
    A=sum(alphas)
    log=A*np.log(4)+loggamma(A+2)+sum(loggamma(v+1) for v in alphas)-sum(loggamma(A-v+2) for v in alphas)
    return np.exp(1j*log.imag)


def physical(b,d,external,N=128,K=12):
    a,q,c,p=external
    phase=period_phase(a,q,b)*period_phase(c,p,-b)/(period_phase(c,q,d)*period_phase(a,p,-d))
    raw,condition=raw_6j(b,d,external,N,K)
    return b*d/np.pi**4*phase*raw,condition
