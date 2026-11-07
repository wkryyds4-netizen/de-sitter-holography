"""Pure two-time neutral selected-state kernels; U has variance one."""
from fractions import Fraction as F
from functools import lru_cache
from math import cosh,tanh,sqrt,pi,log
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.special import roots_legendre

def padd(a,b,c=F(1)):
    out=dict(a)
    for k,v in b.items():out[k]=out.get(k,F(0))+c*v
    return {k:v for k,v in out.items() if v}

@lru_cache(None)
def hermite(n,q):
    if n==0:return {0:F(1)}
    if n==1:return {1:F(1)}
    return padd({k+1:v for k,v in hermite(n-1,q).items()},hermite(n-2,q),-sum(q**j for j in range(n-1)))

@lru_cache(None)
def power_coeff(n,q):
    if n==0:return {0:F(1)}
    out={}
    for j,c in power_coeff(n-1,q).items():
        out[j+1]=out.get(j+1,F(0))+c
        if j:out[j-1]=out.get(j-1,F(0))+c*sum(q**k for k in range(j))
    return out

def transition_power(n,q):
    out={}
    for j,c in power_coeff(n,q).items():out=padd(out,hermite(j,q),c*q**j)
    return out

def chain_moment(a,b,c,q):
    return sum((x*y*power_coeff(i+b+j,q).get(0,F(0))
                for i,x in transition_power(a,q).items()
                for j,y in transition_power(c,q).items()),F(0))

def pair_moment(a,b,q,r=None):
    if r is None:r=q
    result=F(0)
    aa,bb=power_coeff(a,q),power_coeff(b,q)
    for j,x in aa.items():
        fac=F(1)
        for k in range(1,j+1):fac*=sum(q**i for i in range(k))
        result+=x*bb.get(j,F(0))*r**j*fac
    return result

def predicted(u,v):
    delta=u-v;g=cosh(delta)**-2
    R=g*(delta*tanh(delta)-3*tanh(delta)**2)/2
    B=cosh(u)**-2*cosh(v)**-2
    return dict(trace_leading=g,response=R+B/2,selected_first=R/2+B/2,mixture_first=B/4)

class SpectralKernel:
    def __init__(self,eps,L):
        self.eps=eps;self.q=np.exp(-eps);self.L=L
        self.E,self.V=eigh_tridiagonal(np.zeros(L),np.sqrt(-np.expm1(-eps*np.arange(1,L))/-np.expm1(-eps)))
        self.V*=np.where(self.V[0,:]<0,-1.,1.)[None,:]
        self.root=self.V[0,:];self.D=np.exp(-eps*np.arange(L))
    def conditional(self,t):
        return self.V.T@(self.D*(self.V@(self.root*np.exp(1j*t*self.E/sqrt(self.eps)))))
    def pair(self,t,d=1,tilt=False):
        ph=np.exp(1j*t*self.E/sqrt(self.eps));a=self.V@(self.root*ph)
        if tilt:
            weight=1-sqrt(self.q)+sqrt(self.q)*self.E**2
            b=self.V@(self.root*weight*ph)
        else:b=a
        return np.vdot(a,self.D**d*b)
    def chain(self,u,v):
        terms=self.conditional(u)*self.conditional(-v)*np.exp(1j*(v-u)*self.E/sqrt(self.eps))
        value=np.sum(terms)
        return value,float(np.sum(abs(terms))/max(abs(value),1e-300))
    def selected(self,u,v):
        b,condition=self.chain(u,v)
        base=self.pair(u-v,tilt=True)
        return base+sqrt(self.q)*self.eps/4*(1+self.q)*b,b,condition

def angular_chain(q,points,terms,u,v):
    x,w=roots_legendre(points);theta=(x+1)*pi/2;w=w*pi/2
    logmu=np.log(4*np.sin(theta)**2/(2*pi))
    for n in range(1,terms+1):
        z=q**n;logmu+=np.log1p(-z)+np.log(1-2*z*np.cos(2*theta)+z*z)
    mu=w*np.exp(logmu);E=2*np.cos(theta)/sqrt(1-q);eps=-log(q)
    cp=np.cos(theta[:,None]+theta[None,:]);cm=np.cos(theta[:,None]-theta[None,:])
    logK=np.zeros_like(cp)
    for n in range(terms):
        z=q**(n+1)
        logK+=np.log1p(-q**(n+2))-np.log(1-2*z*cp+z*z)-np.log(1-2*z*cm+z*z)
    K=np.exp(logK)
    a=K@(mu*np.exp(1j*u*E/sqrt(eps)));c=K@(mu*np.exp(-1j*v*E/sqrt(eps)))
    value=np.sum(mu*a*c*np.exp(1j*(v-u)*E/sqrt(eps)))
    return value,float(mu.sum()),float(np.max(abs(K@mu-1)))
