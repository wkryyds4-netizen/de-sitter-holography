"""Mehler bound for a Gaussian/finite-Hermite overlap, complex outer argument.

This supplies a uniform enclosure, not a point evaluation. The exact inner
integral is holomorphic on any ball accepted by outer_guard.
"""
from residual_overlap_rectangle import outer_guard
from flint import arb,acb

def gaussian_outer_bound(y,n,normw,terms):
    if not outer_guard(y):return acb('nan')
    pi=arb.pi();u=pi.sqrt()*y.real;v=pi.sqrt()*y.imag
    cosine=abs((u/2).cos())
    if cosine.contains(0):return acb('nan')
    phase=(2*terms*abs(v)).exp()/cosine**terms
    for j in range(terms):
        a=(u-(2*j+1)*pi).exp();den=abs(1-a)
        if den.contains(0):return acb('nan')
        phase*=(1+a)/den
    candidates=[]
    for k in range(1,9):
        r=1-arb(1)/(2**k)
        bound=normw*r**(-(n-1))/(pi*(1-r)).sqrt()
        bound*=(-y.real*y.real/(1+r)+y.imag*y.imag/(1-r)).exp()*phase
        candidates.append(bound.upper())
    best=min(candidates,key=float)
    return acb(arb(0,best),arb(0,best))
