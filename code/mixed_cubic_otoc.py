"""Unequal-probe BINT crossed series, with q^6 restored; Js=1."""
from pathlib import Path
from itertools import combinations
from math import exp,log,sqrt,factorial
from collections import Counter
import json,hashlib,platform
import numpy as np
import sympy as s
from singlet_coefficient import pairings
from singlet_multimark import crossing
from singlet_otoc import taylor_tail
ROOT=Path(__file__).resolve().parent.parent
q=s.Symbol('q')

def moment(arcs):
    word=''.join(c+'H'*k for c,k in zip('XRXR',arcs))
    hs=tuple(i for i,c in enumerate(word) if c=='H')
    xx=tuple(i for i,c in enumerate(word) if c=='X');rr=tuple(i for i,c in enumerate(word) if c=='R')
    out=Counter()
    for pa in pairings(hs):
        chords=[(c,1) for c in pa]+[(xx,2),(rr,3)]
        power=sum(a*b*crossing(c,d) for (c,a),(d,b) in combinations(chords,2))
        out[power]+=1
    return sum(v*q**k for k,v in out.items())

def coefficient(n):
    out=0;terms=0
    for a in range(n+1):
        for b in range(n-a+1):
            for c in range(n-a-b+1):
                d=n-a-b-c;arcs=(a,b,c,d)
                out+=s.Rational((-1)**(a+c+n//2),factorial(a)*factorial(b)*factorial(c)*factorial(d))*moment(arcs)
                terms+=1
    return s.factor(out),terms

def cutoff(lam,tmax,d1,d2,tol=1e-11):
    qv=exp(-4*lam);a=qv**d1;b=qv**d2;Q=qv**(d1*d2)
    C=1/np.prod([1-qv**j for j in range(1,2*d2)]);K=float(C/(1-a))
    z=2*tmax/sqrt(lam*(1-qv))
    inner=max(1,int(log(tol/(100*(1+K)))/log(b))+1)
    L=inner+max(12,int(z)+2)
    while True:
        epsv=2*taylor_tail(z,L);epsp=2*taylor_tail(z,L-inner)
        delta=epsv+epsp+2*b**inner;bound=2*Q*K*delta
        if bound<tol:break
        L+=1
    return dict(L=L,inner=inner,K=K,RH=2/sqrt(lam*(1-qv)),vacuum_bound=epsv,propagator_bound=epsp,
        source_tail=b**inner,otoc_truncation_bound=bound)

class Kernel:
    def __init__(self,lam,L,d1=2,d2=3):
        self.q=exp(-4*lam);self.Q=self.q**(d1*d2);self.L=L
        n=np.arange(1,L);off=np.sqrt(-np.expm1(-4*lam*n)/(-np.expm1(-4*lam)*lam))
        self.e,self.u=np.linalg.eigh(np.diag(off,1)+np.diag(off,-1))
        self.an=np.exp(-4*lam*d1*np.arange(L))
        logf=np.r_[0.,np.cumsum(np.log(-np.expm1(-4*lam*n)))];i=np.arange(L)
        logci=sum(np.log(-np.expm1(-4*lam*(i+j))) for j in range(1,2*d2))-sum(log(-np.expm1(-4*lam*j)) for j in range(1,2*d2))-4*lam*d1*i
        l=i[:,None];ii=i[None,:];index=l+ii;mask=index<L;safe=np.minimum(index,L-1)
        self.index=safe
        self.weight=np.where(mask,np.exp(.5*logci[None,:]-4*lam*d2*l+.5*(logf[safe]-logf[l])),0.)
    def at(self,t):
        minus=(self.u*np.exp(-1j*t*self.e)[None,:])@self.u.T;vac=minus[:,0]
        evolved=minus.conj()@(self.weight*vac[self.index])
        val=self.Q*float(np.sum(self.an[:,None]*abs(evolved)**2))
        return dict(time=t,otoc=val,commutator_squared=2*(1-val),vacuum_residual=float(abs(np.vdot(vac,vac).real-1)))

def main():
    target=ROOT/'results/susskind-mixed-cubic/otoc.json'
    if target.exists():raise FileExistsError(target)
    coefs=[]
    for n in (0,2,4,6):
        c,count=coefficient(n);coefs.append(dict(order=n,coefficient_in_sigma_t=str(c),arc_terms=count))
        print(n,c,flush=True)
    c2=s.sympify(coefs[1]['coefficient_in_sigma_t']);assert s.expand(c2+2*q**6*(1-q*q)*(1-q**3))==0
    rows=[]
    for lam in [.2,.1,.05]:
        bounds=[cutoff(lam,1,*ds) for ds in [(2,3),(3,2)]];L=max(b['L'] for b in bounds)
        k=Kernel(lam,L);ref=Kernel(lam,L+16);swap=Kernel(lam,L,3,2)
        points=[]
        for t in [0.,.01,.25,.5,1.]:
            v=k.at(t);other=ref.at(t);sv=swap.at(t)
            v.update(cutoff_difference=abs(v['otoc']-other['otoc']),swap_difference=abs(v['otoc']-sv['otoc']))
            if t<=.01:
                series=sum(float(s.sympify(c['coefficient_in_sigma_t']).subs(q,exp(-4*lam)))*(t/sqrt(lam))**c['order'] for c in coefs)
                v['sixth_order_series']=series;v['series_difference']=abs(v['otoc']-series)
                assert v['series_difference']<2e-12
            assert v['cutoff_difference']<1e-11 and v['swap_difference']<1e-11
            assert 0<=v['commutator_squared']<=2
            points.append(v)
        rows.append(dict(lambda_value=lam,cutoff=L,refined_cutoff=L+16,bounds=bounds,points=points))
        print(lam,L,points[-1],flush=True)
    inputs=['code/mixed_cubic_otoc.py','code/singlet_otoc.py','code/singlet_multimark.py','code/singlet_coefficient.py']
    data=dict(status='passed',conventions='Js=1; beta=mu=0; DS first; normalized X and Rhat; only one crossed channel',
      coefficients=coefs,records=rows,environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__),
      error_scope='Analytic Jacobi truncation bound; floating point assessed by cutoff, species-swap and independent chord Taylor comparisons, not interval certified',
      input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    target.write_text(json.dumps(data,indent=2)+'\n')
if __name__=='__main__':main()
