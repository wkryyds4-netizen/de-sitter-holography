"""Pure fixed-lambda projected-bilocal preparation Gram; no scientific main."""
from math import sqrt,log,factorial,exp,expm1
import numpy as np
from scipy.linalg import eigh_tridiagonal
from numpy.polynomial.legendre import leggauss
from neutral_tower_projection import tower,qi

def projected_tail(q,lam,t,J):
    """L2 tail of P B, exact-arithmetic bound inherited from mean stage."""
    r=1/sqrt(1-q);a=q*r*abs(t);b=3*r*abs(t);force=4*lam*q/(1-q)
    if a>=J+1 or b>=J+2:return float('inf')
    return t*t*sqrt(lam)*(4*lam*a**(J+1)/factorial(J)/(1-a/(J+1))+
                            force*b**(J+1)/factorial(J+1)/(1-b/(J+2)))

def beta_norm_bound(q,lam,t):
    r=1/sqrt(1-q);a=q*r*abs(t);b=3*r*abs(t);force=4*lam*q/(1-q)
    return t*t*(4*lam*a*exp(a)+force*expm1(b))

class PreparationGram:
    def __init__(self,q,J=32,L=64,taylor=80):
        self.q=q;self.lam=-log(q)/4;self.J=J;self.L=L;self.taylor=taylor
        self.e,self.V=eigh_tridiagonal(np.zeros(L),np.sqrt((1-q**np.arange(1,L))/(1-q)))
        self.V*=np.where(self.V[0]<0,-1,1)[None,:]
        self.root=self.V[0];self.transfer=self.V.T@(q**np.arange(L)[:,None]*self.V)
        self.mass=self.root[:,None]*self.root[None,:]*self.transfer
        self.delta=self.e[:,None]-self.e[None,:];self.plus=self.e[:,None]+self.e[None,:]
        x=self.e[:,None];y=self.e[None,:];P={(0,0):np.ones((L,L))}
        for b in range(1,J):P[0,b]=y*P[0,b-1]-(float(qi(b-1,q))*P[0,b-2] if b>1 else 0)
        for a in range(1,J):
            for b in range(J-a):P[a,b]=x*P[a-1,b]-(float(qi(a-1,q))*P[a-2,b] if a>1 else 0)-(q**a*float(qi(b,q))*P[a-1,b-1] if b else 0)
        _,Xi,_,_=tower(J,q,self.lam)
        self.X=[np.zeros((L,L))]+[sum(float(c)*P[ab] for ab,c in Xi[j].items()) for j in range(1,J+1)]
        ac=np.zeros((taylor+1,taylor+1));ac[0,0]=1
        for k in range(taylor):
            for j in range(k+1):
                ac[k+1,j+1]+=ac[k,j]
                if j:ac[k+1,j-1]+=ac[k,j]*(1-q**j)/(1-q)
        self.ac=ac

    def base_block(self,t,ns=32):
        z,w=leggauss(ns);beta=np.zeros((self.L,self.L),dtype=complex)
        for s,weight in zip((z+1)*t/2,w*t/2):
            c=[sum((1j*s)**k*self.ac[k,j]/factorial(k) for k in range(j,self.taylor+1,2)) for j in range(self.J+1)]
            psi=sum(c[j]*self.X[j] for j in range(1,self.J+1))
            beta-=2*weight*(t-s)*np.sinc((t-s)*self.delta/(2*np.pi))*np.exp(-.5j*s*self.plus)*psi
        return beta

    def state(self,base,t,T):
        beta=np.exp(1j*T*self.delta)*base
        # sqrt(mu_i) times the Mehler conditional integral; no division by tiny weights.
        mroot=np.sum(self.root[None,:]*self.transfer*beta,axis=1)
        summands=-sqrt(self.q)*self.lam*self.mass*self.plus*beta
        mean=summands.sum()
        norm2=float(np.sum(self.mass*np.abs(beta)**2).real)
        return dict(beta=beta,mroot=mroot,mean=mean,t=t,T=T,
                    beta_norm=sqrt(max(0,norm2)),tail_beta=projected_tail(self.q,self.lam,t,self.J)/sqrt(self.lam),
                    analytic_beta_bound=beta_norm_bound(self.q,self.lam,t),
                    mean_cancellation=float(np.abs(summands).sum()/max(abs(mean),1e-300)))

    def excess(self,a,b):
        pref=sqrt(self.q)*self.lam**2*(1+self.q)
        gram=pref*np.vdot(a['mroot'],b['mroot'])
        means=np.conj(a['mean'])*b['mean'];connected=gram-means
        # Analytic majorants bound both full and truncated symbols; no numerical norm in this bound.
        da,db=a['tail_beta'],b['tail_beta'];na,nb=a['analytic_beta_bound'],b['analytic_beta_bound']
        product_tail=da*nb+na*db+da*db
        gram_tail=pref*product_tail
        mean_factor=sqrt(self.q)*self.lam*sqrt(2*(1+self.q))
        mean_tail=mean_factor**2*product_tail
        return dict(uncentered=gram,connected=connected,mean_product=means,
                    propagated_degree_bound=gram_tail+mean_tail,
                    subtraction_condition=float((abs(gram)+abs(means))/max(abs(connected),1e-300)))

    def C2(self,T,derivative=0):
        mass2=self.root[:,None]*self.root[None,:]*(self.V.T@(self.q**(2*np.arange(self.L))[:,None]*self.V))
        return np.sum(mass2*(1j*self.delta)**derivative*np.exp(1j*T*self.delta))
