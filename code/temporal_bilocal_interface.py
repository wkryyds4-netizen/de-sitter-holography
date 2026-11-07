"""Exact temporal CAR derivatives and ensemble Majorana-degree weights.

Wick degree here means Hilbert-Schmidt/Clifford degree, not coupling chaos
or the U(N) irreducible covariance label used in earlier notes.
"""
from pathlib import Path
from itertools import combinations
from fractions import Fraction
from math import comb, factorial
import hashlib, json, platform
import numpy as np
import mpmath as mp
from complex_syk_dictionary import fermions, creators, comm

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-bilocal-interface'

def norm(a):
    v=np.sum(abs(a)**2)
    assert abs(v-round(v))<1e-6 and v<2**48
    return Fraction(round(v),a.shape[0])

def grades(H,cs):
    gammas=[g for c in cs for g in [c+c.conj().T,1j*(c.conj().T-c)]]
    dim=len(H);out=[Fraction(0) for _ in range(len(cs)+1)]
    for k in range(0,len(gammas)+1,2):
        for inds in combinations(range(len(gammas)),k):
            g=np.eye(dim,dtype=complex)
            for i in inds:g=g@gammas[i]
            g*=1j**(k*(k-1)//2)
            a=np.trace(g@H)
            assert abs(a.imag)<1e-8 and abs(a.real-round(a.real))<1e-8
            out[k//2]+=Fraction(round(a.real),dim)**2
    assert sum(out)==norm(H)
    return out

def tensors(H,cs):
    mixed=Fraction(0);anomalous=Fraction(0);scalar=Fraction(0);first=Fraction(0)
    for c in cs:
        d=comm(H,c);first+=norm(d)
        for b in cs:
            t=d@b.conj().T+b.conj().T@d
            mixed+=norm(t);scalar+=Fraction(round(np.trace(t).real),len(H))**2+Fraction(round(np.trace(t).imag),len(H))**2
            anomalous+=norm(d@b+b@d)
    return dict(first=first,mixed=mixed,anomalous=anomalous,scalar=scalar)

def algebra_checks():
    rows=[];rng=np.random.default_rng(110630)
    for N,p in [(3,1),(3,2),(4,2),(4,3),(5,2)]:
        cs=fermions(N);basis=creators(cs,p);L=len(basis)
        raw=rng.integers(-2,3,(L,L))+1j*rng.integers(-2,3,(L,L));J=raw+raw.conj().T
        H=sum(J[a,b]*(A@B.conj().T) for a,A in enumerate(basis) for b,B in enumerate(basis))
        w=grades(H,cs);actual=tensors(H,cs)
        expect=dict(first=sum(2*k*v for k,v in enumerate(w)),mixed=sum(4*k*k*v for k,v in enumerate(w)),
                    anomalous=sum(4*k*(k-1)*v for k,v in enumerate(w)),scalar=4*w[1])
        assert actual==expect,(N,p,actual,expect)
        # Full four-anticommutator identity, with odd operators involving H.
        a,b,c,d=cs[0].conj().T,comm(H,cs[1]),comm(H,cs[0].conj().T),cs[-1]
        anti=lambda x,y:x@y+y@x
        rhs=a@anti(b,c)@d-c@anti(d,a)@b+c@a@anti(d,b)-anti(a,c)@b@d
        assert np.array_equal(comm(a@b,c@d),rhs)
        rows.append(dict(N=N,p=p,weights=list(map(str,w)),actual={k:str(v) for k,v in actual.items()}))
    return rows

def combinatorial(N,p):
    # E||H_k||^2 / v, with H=2^p sum J_IJ C_I^dag C_J and E JJ=v swap.
    return [comb(N,p)*comb(p,k)*comb(N-p+k,k) for k in range(p+1)]

def ensemble_checks():
    rows=[]
    for N in range(2,13):
        for p in range(1,N+1):
            weights=combinatorial(N,p)
            via_overlap=[]
            for k in range(p+1):
                val=0
                for r in range(max(0,2*p-N),p+1):
                    if 0<=p-k<=r:
                        val+=comb(p,r)*comb(N-p,p-r)*comb(r,p-k)
                via_overlap.append(comb(N,p)*val)
            assert weights==via_overlap
            t=sum(weights);overlap=sum(comb(N,p)*comb(p,r)*comb(N-p,p-r)*2**r for r in range(max(0,2*p-N),p+1))
            assert t==overlap
            rows.append(dict(N=N,p=p,weights_over_v=list(map(str,weights))))
    # Direct ensemble Hessians from each non-Hermitian monomial; GUE swap
    # turns the quadratic expectation into a sum of Hilbert-Schmidt norms.
    direct=[]
    for N,p in [(3,1),(3,2),(4,2)]:
        cs=fermions(N);bs=creators(cs,p);values={k:Fraction(0) for k in ['first','mixed','anomalous','scalar']}
        for a in bs:
            for b in bs:
                T=2**p*a@b.conj().T
                for k,v in tensors(T,cs).items():values[k]+=v
        weights=combinatorial(N,p)
        expected=dict(first=sum(2*k*v for k,v in enumerate(weights)),mixed=sum(4*k*k*v for k,v in enumerate(weights)),
                      anomalous=sum(4*k*(k-1)*v for k,v in enumerate(weights)),scalar=4*weights[1])
        assert values==expected,(N,p,values,expected)
        direct.append(dict(N=N,p=p,monomials=len(bs)**2,values_over_v={k:str(v) for k,v in values.items()}))
    return rows,direct

def limits():
    mp.mp.dps=80;rows=[]
    for lam in ['0.25','0.5','1']:
        l=mp.mpf(lam)
        for p in [8,16,32,64,128]:
            N=int(mp.mpf(p*p)/l)
            ints=combinatorial(N,p);w=[mp.mpf(v) for v in ints]
            t=sum(w);mu1=sum(k*v for k,v in enumerate(w))/t;mu2=sum(k*k*v for k,v in enumerate(w))/t
            v=mp.mpf(N)/(p*p)*(mp.factorial(p)/mp.mpf(N)**p)**2
            h2=v*t
            rows.append(dict(N=N,p=p,lambda_value=lam,H2=str(h2),H2_limit=str(1/l),
                mixed_per_N=str(4*h2*mu2/N),anomalous_per_N=str(4*h2*(mu2-mu1)/N),
                centered_mixed_per_N=str(4*v*(sum(k*k*z for k,z in enumerate(w))-w[1])/N),
                scalar_fraction=str(w[1]/sum(k*k*z for k,z in enumerate(w))),
                mixed_scaled=str(mu2/p**2),anomalous_scaled=str((mu2-mu1)/p**2),
                size_deficit=str(p-mu1),size_deficit_limit=lam))
    return rows

def main():
    OUT.mkdir(parents=True,exist_ok=True);target=OUT/'micro.json'
    if target.exists():raise FileExistsError(target)
    algebra=algebra_checks();ensemble,direct=ensemble_checks();asymptotic=limits()
    files=['code/temporal_bilocal_interface.py','code/complex_syk_dictionary.py']
    result=dict(status='passed',seed=110630,algebra=algebra,combinatorial=ensemble,direct_ensemble=direct,
        double_scaling=asymptotic,environment=dict(python=platform.python_version(),numpy=np.__version__,mpmath=mp.__version__),
        error_scope='Exact finite CAR and integer combinatorics; 80-digit evaluations of finite combinatorial formulas. Only zero-time curvature is taken to the DS limit, not a uniform finite-time expansion.',
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print('Exact operator cases',len(algebra),'combinatorial cases',len(ensemble),'direct ensembles',len(direct),'DS rows',len(asymptotic))
    for r in asymptotic:
        if r['p']==128:print(r['lambda_value'],r['mixed_per_N'],r['anomalous_per_N'],flush=True)

if __name__=='__main__':main()
