"""Exact finite-N single contractions of the quartic memory core.

The independent-copy covariance is the full embedded GUE, including overlaps.
Matrices use dyadic Gaussian integers, checked before Fraction Gram sums.
This is ordinary coupling chaos, not the noncommutative q-Wick product.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import hashlib, json, platform, time
import numpy as np
import scipy
from scipy.sparse import csr_matrix
import sympy as sp
from singlet_wick_moments import gue_basis

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-quadratic-contraction'


def adj(x): return x.conj().T
def anti(x,y): return x@y+y@x


class Contraction:
    def __init__(self,N,p,m):
        self.N,self.p,self.m=N,p,m
        cs,vs,rs=gue_basis(N,p)
        charge=np.array([i.bit_count() for i in range(1<<N)])
        ix=np.flatnonzero(charge==m); iz=np.flatnonzero(charge==m-1)
        self.hs=vs[:,ix[:,None],ix];self.ks=vs[:,iz[:,None],iz]
        self.cs=np.array([c[np.ix_(iz,ix)] for c in cs]);self.rs=rs
        self.d=len(ix);self.dm=len(iz)
        self.C=comb(m,p)*comb(N-m+p,p) if m>=p else 0
        self.Cm=comb(m-1,p)*comb(N-m+1+p,p) if m>p else 0
        self.a=(2*p-m)*self.C+m*self.Cm
        self.phis={}
        for label,left,right in [('hh',self.hs,self.hs),('kk',self.ks,self.ks),('hk',self.hs,self.ks)]:
            # twice Phi, with row-major input/output flattening
            ten=np.einsum('a,aij,akl->iljk',rs,left,right,optimize=True)
            self.phis[label]=csr_matrix(ten.reshape(left.shape[1]*right.shape[1],-1))
        assert np.array_equal(self.phi('hh',np.eye(self.d)),self.C*np.eye(self.d))
        assert np.array_equal(self.phi('kk',np.eye(self.dm)),self.Cm*np.eye(self.dm))

    def phi(self,label,x):
        return (self.phis[label]@x.ravel()).reshape(x.shape)/2

    def rem(self,x): return sum(adj(c)@x@c for c in self.cs)

    def gamma(self,x,y):
        return sum(self.phi('hk',adj(c)@x)@y@c for c in self.cs)

    def A(self,h,k): return self.rem(k@k)+(2*self.p-self.m)*h@h

    def core(self,h,k):
        # Independent commutator definition, not the formula used by contract.
        ds=[k@c-c@h for c in self.cs]
        fs=[k@d-d@h for d in ds]
        a=sum(adj(d)@d for d in ds)
        return sum(adj(f)@f for f in fs)+(h@h@a-2*h@a@h+a@h@h)/4

    def raw_word_core(self,h,k):
        h2=h@h;k2=k@k;r2=self.rem(k2);r3=self.rem(k2@k)
        return self.rem(k2@k2)-2*anti(h,r3)+1.25*anti(h2,r2)+3.5*h@r2@h+(4*self.p-3*self.m)*h2@h2

    def contract(self,h,k):
        """Q(H)=E_K [t^2] E_core(H+tK)=Delta E_core/2."""
        h2=h@h;k2=k@k;im=np.eye(self.dm);r2=self.rem(k2)
        ph=self.phi('hh',h);pk=self.phi('kk',k)
        u=self.gamma(k,im)+self.gamma(im,k)
        v=self.gamma(im,k2)+self.gamma(k,k)+self.gamma(k2,im)
        uh=sum(self.phi('hk',h@adj(c)@k)@c+self.phi('hk',h@adj(c))@k@c for c in self.cs)
        r4=self.rem(3*self.Cm*k2+anti(k,pk)+self.phi('kk',k2))
        hr3=anti(h,self.rem(2*self.Cm*k+pk))+v+adj(v)
        h2r2=2*self.C*r2+2*self.m*self.Cm*h2+h@u+uh+adj(u)@h+adj(uh)
        hr2h=self.phi('hh',r2)+u@h+h@adj(u)+self.m*self.Cm*h2
        h4=3*self.C*h2+anti(h,ph)+self.phi('hh',h2)
        return r4-2*hr3+1.25*h2r2+3.5*hr2h+(4*self.p-3*self.m)*h4

    def tlinear(self,h,k):
        im=np.eye(self.dm)
        u=self.gamma(k,im)+self.gamma(im,k)+(2*self.p-self.m)*(self.phi('hh',h)+self.C*h)
        return self.rem(2*self.Cm*k+self.phi('kk',k))+(3*self.p-self.m)*(2*self.C*h+self.phi('hh',h))-3*(self.a*h+(u+adj(u))/2)

    def direct_contract(self,h,k):
        value=np.zeros_like(h)
        base=self.core(h,k)
        for r,hh,kk in zip(self.rs,self.hs,self.ks):
            # Quartic even polarization removes the t^4 term exactly.
            value+=r*(self.core(h+hh,k+kk)+self.core(h-hh,k-kk)-2*base-2*self.core(hh,kk))/4
        return value


def integer_array(x,scale=8):
    a=scale*x
    assert np.max(abs(a))<2**35
    assert np.array_equal(a.real,np.rint(a.real)) and np.array_equal(a.imag,np.rint(a.imag))
    return a.real.astype(np.int64),a.imag.astype(np.int64)


def inner_exact(a,b):
    ar,ai=integer_array(a);br,bi=integer_array(b)
    # Python integers prevent integer overflow in the exact Gram sums.
    return F(sum(int(x)*int(y) for x,y in zip(ar.flat,br.flat))+sum(int(x)*int(y) for x,y in zip(ai.flat,bi.flat)),64*a.shape[0])


def evaluate(N,p,m):
    obj=Contraction(N,p,m);count=len(obj.rs)
    x=np.array([(7*i+3)%11-5 for i in range(count)])
    h=np.einsum('a,aij->ij',x,obj.hs);k=np.einsum('a,aij->ij',x,obj.ks)
    assert np.array_equal(obj.core(h,k),obj.raw_word_core(h,k))
    assert np.array_equal(obj.contract(h,k),obj.direct_contract(h,k))
    tvalues=[obj.tlinear(hh,kk) for hh,kk in zip(obj.hs,obj.ks)]
    tbar=sum(F(int(r),2)*inner_exact(hh,t) for r,hh,t in zip(obj.rs,obj.hs,tvalues))/obj.C
    # J_H F is included after clearing tbar's exact denominator.
    den=tbar.denominator;num=tbar.numerator
    def values(hh,kk):
        ff=den*obj.tlinear(hh,kk)-num*hh
        return [hh@hh,obj.A(hh,kk),(hh@ff+ff@hh)/2,obj.contract(hh,kk)]
    diagonal=[values(hh,kk) for hh,kk in zip(obj.hs,obj.ks)]
    gram=[[F(0) for _ in range(4)] for _ in range(4)]
    eq=np.zeros_like(h);maxentry=0
    for a in range(count):
        eq+=obj.rs[a]*diagonal[a][-1]/2
        for b in range(a,count):
            if a==b:bs=diagonal[a]
            else:
                mixed=values(obj.hs[a]+obj.hs[b],obj.ks[a]+obj.ks[b])
                bs=[(z-za-zb)/2 for z,za,zb in zip(mixed,diagonal[a],diagonal[b])]
            factor=F(int(obj.rs[a]*obj.rs[b]),2)*(1 if a==b else 2)
            for i in range(4):
                maxentry=max(maxentry,float(np.max(abs(8*bs[i]))))
                for j in range(i,4):gram[i][j]+=factor*inner_exact(bs[i],bs[j])
    for i in range(4):
        for j in range(i):gram[i][j]=gram[j][i]
    # Undo the deliberate denominator clearing of F.
    for i in range(4):
        for j in range(4):gram[i][j]/=den**(int(i==2)+int(j==2))
    assert np.array_equal(eq,np.trace(eq)/obj.d*np.eye(obj.d))
    mean=F(float(np.trace(eq).real),1) if False else F(int(np.trace(eq).real*8),8*obj.d)/2
    G=sp.Matrix([[sp.Rational(z.numerator,z.denominator) for z in row] for row in gram])
    results=[]
    for ids in [[0],[0,1],[0,1,2]]:
        small=G.extract(ids,ids);g=G.extract(ids,[3])
        if small.det()==0:
            results.append(dict(basis=ids,singular=True));continue
        coeff=small.inv()*g;norm=sp.factor(G[3,3]-(g.T*coeff)[0])
        assert norm>=0
        results.append(dict(basis=ids,coefficients=list(map(str,coeff)),residual_norm_squared=str(norm)))
    q=obj.contract(h,k)-2*float(mean)*np.eye(obj.d)
    comm=h@q-q@h
    witness=inner_exact(comm,comm)
    return dict(N=N,p=p,m=m,real_gaussian_coordinates=count,quadratic_coefficients=count*(count+1)//2,
                covariance_w=1,core_mean=str(mean),tbar=str(tbar),
                exact_gram=[[str(z) for z in row] for row in gram],
                basis=['H^2-C','A-a','J_H F','W2 core'],projections=results,
                deterministic_commutator_norm_squared=str(witness),
                deterministic_couplings=x.tolist(),max_scaled_coefficient=maxentry,
                raw_identity_and_independent_polarization=True)


def main():
    OUT.mkdir(exist_ok=True);target=OUT/'finite-contractions.json'
    if target.exists():raise FileExistsError(target)
    rows=[];start=time.time()
    for args in [(3,1,2),(4,2,3),(5,2,3)]:
        row=evaluate(*args);rows.append(row)
        print(json.dumps(dict(system=args,projections=row['projections'],commutator=row['deterministic_commutator_norm_squared'])),flush=True)
    paths=['code/memory_quadratic_contraction.py','code/singlet_wick_moments.py','code/complex_syk_dictionary.py']
    result=dict(status='passed',rows=rows,elapsed_seconds=time.time()-start,
                conventions='Finite N; fixed charge; annealed normalized trace; standard CAR; full overlapping embedded GUE; w=1. No DS extrapolation.',
                method='Exact dyadic coordinate coefficients and Fraction Gram. Independent quartic polarization of commutator-defined core.',
                environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__),
                input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in paths})
    target.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
