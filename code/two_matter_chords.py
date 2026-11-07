"""Ordered two-matter chord Gram recursion (Lin--Stanford App. B.2).

The crossed form pairs the order (d1,d2) on the bra with (d2,d1)
on the ket. It is an overlap form, not a unitary scattering matrix.
"""
from functools import lru_cache
import mpmath as mp

def basis(n):return [(l,m,n-l-m) for l in range(n+1) for m in range(n-l+1)]

class TwoMatter:
    def __init__(self,h,d1,d2):
        self.h,self.d1,self.d2=map(mp.mpf,[h,d1,d2]);self.q=mp.exp(-self.h)
    def qi(self,n):return -mp.expm1(-self.h*n)/(-mp.expm1(-self.h))
    def lower_basis(self,t,side,reverse=False):
        l,m,r=t;a,b=(self.d2,self.d1) if reverse else (self.d1,self.d2)
        if side=='L':weights=[1,self.q**(a+l),self.q**(a+b+l+m)]
        else:weights=[self.q**(a+b+m+r),self.q**(b+r),1]
        out={}
        for j,n in enumerate(t):
            if n:
                u=list(t);u[j]-=1;out[tuple(u)]=self.qi(n)*weights[j]
        return out
    @lru_cache(None)
    def gram(self,u,v,cross=False):
        if sum(u)!=sum(v):return mp.mpf(0)
        for j,side in [(0,'L'),(2,'R')]:
            if u[j]:
                a=list(u);a[j]-=1
                return mp.fsum(c*self.gram(tuple(a),b,cross) for b,c in self.lower_basis(v,side,cross).items())
        for j,side in [(0,'L'),(2,'R')]:
            if v[j]:
                b=list(v);b[j]-=1
                return mp.fsum(c*self.gram(a,tuple(b),cross) for a,c in self.lower_basis(u,side,False).items())
        n=u[1];fac=mp.fprod(self.qi(j) for j in range(1,n+1))
        return fac*(self.q**(self.d1*self.d2+n*(self.d1+self.d2)) if cross else 1)
    def matrix(self,n,cross=False):
        states=basis(n)
        return mp.matrix([[self.gram(u,v,cross) for v in states] for u in states])
    def primary(self,k):
        if not k:return mp.matrix([1])
        states,below=basis(k),basis(k-1)
        M=mp.matrix(2*len(below),len(states))
        for j,v in enumerate(states):
            for s,side in enumerate(['L','R']):
                for u,c in self.lower_basis(v,side).items():M[s*len(below)+below.index(u),j]=c
        mid=states.index((0,k,0));cols=[j for j in range(len(states)) if j!=mid]
        A=mp.matrix([[M[i,j] for j in cols] for i in range(M.rows)])
        rhs=mp.matrix([-M[i,mid] for i in range(M.rows)])
        # mpmath's unpivoted Householder routine hits an exact zero pivot
        # for k=3 with this sparse ordering. Small high-precision normal
        # equations are used here; the full residual and precision change
        # are checked separately, not inferred from solver success.
        sol=mp.lu_solve(A.T*A,A.T*rhs)
        v=mp.matrix(len(states),1);v[mid]=1
        for j,z in zip(cols,sol):v[j]=z
        v/=mp.sqrt((v.T*self.matrix(k)*v)[0])
        overlap=(v.T*self.matrix(k))[mid]
        if overlap<0:v=-v
        assert mp.norm(M*v)<mp.mpf('1e-60')
        return v
