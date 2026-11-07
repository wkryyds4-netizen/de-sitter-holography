"""Open-chord transfer, oldest endpoint first. No file mutation on import.

Basis strings use H (Hamiltonian opening) and P (selected creation).
A,a are Hamiltonian creation/annihilation; P,m are selected creation/annihilation.
H=A+a. Coefficients support exact Fractions or floating point.
"""
from collections import defaultdict
from functools import lru_cache

def clean(v):return {w:c for w,c in v.items() if c}
def add(*vs):
    out=defaultdict(int)
    for v in vs:
        for w,c in v.items():out[w]+=c
    return clean(out)
def scale(c,v):return clean({w:c*x for w,x in v.items()})
def step(v,op,gamma):
    q=gamma**4;out=defaultdict(int)
    for w,c in v.items():
        n=len(w)
        if op in ('H','A'):out[w+'H']+=c
        if op=='P':out[w+'P']+=c
        if op in ('H','a','m') and n:
            for j,t in enumerate(w):
                if op=='m':weight=q**(n-1)*(gamma if t=='H' else 1)
                else:weight=q**(n-1-j) if t=='H' else gamma*q**(n-1)
                out[w[:j]+w[j+1:]]+=c*weight
    return clean(out)
def moment(word,gamma):
    v={'':1}
    for op in reversed(word):v=step(v,{'H':'H','+':'P','-':'m'}[op],gamma)
    return v.get('',0)
@lru_cache(None)
def basis_inner(u,v,gamma):
    if len(u)!=len(v):return 0
    state={v:1}
    for letter in reversed(u):state=step(state,'a' if letter=='H' else 'm',gamma)
    return state.get('',0)
def inner(u,v,gamma):return sum(c*d*basis_inner(a,b,gamma) for a,c in u.items() for b,d in v.items())
def qnum(n,q):return sum(q**j for j in range(n))
def qfac(n,q):
    result=1
    for j in range(1,n+1):result*=qnum(j,q)
    return result
def primitive(v,gamma):
    """Finite-grade q-binomial projection to ker a."""
    if not v:return {}
    degree=max(map(len,v));q=gamma**4;out={};lower=v
    for k in range(degree+1):
        raised=lower
        for _ in range(k):raised=step(raised,'A',gamma)
        factor=(-1)**k*q**(k*(k-1)//2)/qfac(k,q)
        out=add(out,scale(factor,raised))
        lower=step(lower,'a',gamma)
    return out
