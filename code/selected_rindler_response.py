"""Pure selected-energy response helpers; q=e^-h, K=H/(2 J_s)."""
from functools import lru_cache
from math import comb, factorial, sqrt, tanh, cosh
from fractions import Fraction
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply

def add(a,b):
    c=dict(a)
    for k,v in b.items():c[k]=c.get(k,0)+v
    return {k:v for k,v in c.items() if v}

def mul(a,b):
    c={}
    for i,v in a.items():
        for j,w in b.items():c[i+j]=c.get(i+j,0)+v*w
    return {k:v for k,v in c.items() if v}

def scale(a,c):return {k:v*c for k,v in a.items() if v*c}

@lru_cache(None)
def power(n):
    if not n:return ({0:1},)
    prev=power(n-1);out=[{} for _ in range(n+1)]
    for j,c in enumerate(prev):
        out[j+1]=add(out[j+1],c)
        if j:out[j-1]=add(out[j-1],mul(c,{k:1 for k in range(j)}))
    return tuple(out)

@lru_cache(None)
def qfactorial(n):
    return {0:1} if n==0 else mul(qfactorial(n-1),{k:1 for k in range(n)})

@lru_cache(None)
def energy_pair(a,b):
    out={}
    for j in range(min(a,b)+1):
        term=mul(mul(power(a)[j],power(b)[j]),qfactorial(j))
        out=add(out,{k+2*j:v for k,v in term.items()})
    return out

def gap_difference(n):
    out={}
    for a in range(n+1):
        term=add(energy_pair(a+2,n-a),scale(energy_pair(a,n-a),-1))
        out=add(out,scale(term,(-1)**(n-a)*comb(n,a)))
    return out

def exponential_coefficient(poly,k,shift=Fraction(1,2)):
    return sum((Fraction(v)*(-Fraction(j)-shift)**k for j,v in poly.items()),Fraction(0))/factorial(k)

def response(u):
    return (u*tanh(u)-5*tanh(u)**2)/cosh(u)**4

def packets(h,L,times):
    q=np.exp(-h)
    off=np.sqrt(-np.expm1(-h*np.arange(1,L))/(-np.expm1(-h)*h))
    K=diags([off,off],[-1,1],shape=(L,L),format='csr')
    seeds=np.zeros((L,2),dtype=complex);seeds[0,0]=1;seeds[2,1]=1
    times=np.asarray(times)
    assert times[0]==0 and np.allclose(times,np.linspace(times[0],times[-1],len(times)))
    wave=expm_multiply(1j*K,seeds,start=0,stop=float(times[-1]),num=len(times),traceA=0)
    weights=np.exp(-2*h*np.arange(L))
    a,b=wave[:,:,0],wave[:,:,1]
    c0=np.sum(abs(a)**2*weights,axis=1)
    # Compute the difference directly, avoiding subtracting two nearly equal correlators.
    diff=sqrt(q*(1+q))*np.sum(np.conj(b)*a*weights,axis=1)
    cond=np.sum(abs(np.conj(b)*a)*weights,axis=1)*sqrt(q*(1+q))/np.maximum(abs(diff),1e-300)
    return c0,diff,cond,float(np.max(abs(np.sum(abs(a)**2,axis=1)-1)))

def packet_tail_log(h,L,u):
    """Common exact-arithmetic norm tail; -infinity at zero time."""
    from math import log,lgamma
    if not u:return float('-inf')
    q=np.exp(-h);z=2*abs(u)/sqrt(h*(1-q))
    return log(4*(1+sqrt(q*(1+q))))+z+(L-2)*log(z)-lgamma(L-1)
