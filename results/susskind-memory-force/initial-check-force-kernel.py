"""Isolated-tag CAR, exact chord identities and independent force-kernel evaluation."""
from pathlib import Path
from itertools import combinations
from functools import lru_cache
from math import sqrt,exp,pi
import json,hashlib,platform
import numpy as np
import sympy as s
import mpmath as mp
from complex_syk_dictionary import fermions,comm
from singlet_coefficient import pairings
from singlet_otoc import taylor_tail
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-force'
q,r=s.symbols('q r')

@lru_cache(None)
def moment(u,v):
    out=0
    for pa in pairings(tuple(range(u+v))):
        cross=sum(a<c<b<d or c<a<d<b for j,(a,b) in enumerate(pa) for c,d in pa[j+1:])
        bridge=sum(a<u<=b for a,b in pa)
        out+=q**cross*r**bridge
    return s.expand(out)

def tagged_car():
    N,p=5,2;cs=fermions(N);d=2**N;rng=np.random.default_rng(110634)
    subsets=list(combinations(range(1,N),p));tagged=[(0,j) for j in range(1,N)]
    def create(I):
        a=np.eye(d,dtype=complex)
        for i in I:a=a@cs[i].conj().T
        return a
    def wordtrace(word):
        a=np.eye(d,dtype=complex)
        for x in word:a=a@x
        return np.trace(a)
    rows=[];nonzero=0;negative=0
    for length in [4,6,8]:
        for case in range(40):
            gap=1+case%(length-1);pair=(0,gap)
            remaining=tuple(i for i in range(length) if i not in pair)
            ps=list(pairings(remaining));pa=ps[case%len(ps)]
            I=subsets[rng.integers(len(subsets))];J=tagged[rng.integers(len(tagged))]
            A=create(I)@create(J).conj().T
            word=[None]*length;word[0]=A;word[gap]=A.conj().T
            for a,b in pa:
                I=subsets[rng.integers(len(subsets))];J=subsets[rng.integers(len(subsets))]
                B=create(I)@create(J).conj().T
                word[a]=B;word[b]=B.conj().T
            base=wordtrace(word);force=list(word)
            force[0]=-comm(word[0],cs[0].conj().T);force[gap]=comm(word[gap],cs[0])
            actual=wordtrace(force);assert actual==2*base
            nonzero+=int(base!=0);negative+=int(base.real<0)
            rows.append(dict(length=length,gap=gap,base_trace=[base.real,base.imag],force_trace=[actual.real,actual.imag]))
    assert nonzero>=10
    return dict(N=N,p=p,seed=110634,rows=rows,nonzero=nonzero,negative=negative)

def chord_identities():
    rows=[]
    for total in [0,2,4,6,8]:
        for u in range(total+1):
            v=total-u
            weighted=moment(u+2,v)+moment(u,v+2)-2*moment(u+1,v+1)
            assert s.expand(weighted.subs(r,1))==0
            deriv=s.expand(-s.diff(weighted,r).subs(r,1))
            expected=2*moment(u,v).subs(r,q)
            assert s.expand(deriv-expected)==0
            rows.append(dict(u=u,v=v,weighted_derivative=str(deriv),q_probe=str(s.expand(expected))))
    return rows

def mehler_limit():
    mp.mp.dps=65;rows=[]
    for qv in ['0.2','0.6','0.9']:
        qq=mp.mpf(qv)
        for aa,bb in [('0.4','1.1'),('1.2','2.1'),('0.8','0.83')]:
            a,b=mp.mpf(aa),mp.mpf(bb)
            def kern(rr):
                return mp.qp(rr*rr,qq)/abs(mp.qp(rr*mp.exp(1j*(a+b)),qq)*mp.qp(rr*mp.exp(1j*(a-b)),qq))**2
            de2=4*(mp.cos(a)-mp.cos(b))**2/(1-qq)
            target=2*kern(qq);errors=[]
            for eps in ['1e-5','1e-8','1e-11']:
                e=mp.mpf(eps);val=de2*kern(1-e)/e
                errors.append(float(abs(val/target-1)))
            assert errors[-1]<1e-8 and errors[-1]<errors[0]
            rows.append(dict(q=qv,theta=aa,phi=bb,sigma2=1,epsilons=['1e-5','1e-8','1e-11'],relative_errors=errors))
    return rows

def jacobi(lam,times,L):
    qq=exp(-4*lam);k=np.arange(1,L)
    off=np.sqrt(-np.expm1(-4*lam*k)/(-np.expm1(-4*lam)*lam))
    en,U=np.linalg.eigh(np.diag(off,1)+np.diag(off,-1))
    psi=U@(U[0,:,None]*np.exp(-1j*en[:,None]*np.array(times)[None,:]))
    return np.sum(np.exp(-4*lam*np.arange(L))[:,None]*abs(psi)**2,axis=0)

def quadrature(lam,times,nodes):
    qq=exp(-4*lam);x,w=np.polynomial.legendre.leggauss(nodes)
    th=(x+1)*pi/2;w=w*pi/2
    rho=np.ones(nodes)/(2*pi);num=1.;den=np.ones((nodes,nodes));k=0
    while qq**k/(1-qq)>1e-18:
        rho*=abs(1-qq**k*np.exp(2j*th))**2*(1-qq**(k+1))
        v=qq**(k+1);num*=1-qq**(k+2)
        den*=(1-2*v*np.cos(th[:,None]-th[None,:])+v*v)*(1-2*v*np.cos(th[:,None]+th[None,:])+v*v)
        k+=1
    measure=(rho*w)[:,None]*(rho*w)[None,:]*num/den
    en=2*np.cos(th)/sqrt(lam*(1-qq));de=en[:,None]-en[None,:]
    return dict(nodes=nodes,product_terms=k,values=[float(np.sum(measure*np.cos(t*de))) for t in times],
      mass=float(measure.sum()),second_moment=float(np.sum(measure*de*de)))

def curves():
    times=[0,.25,.5,1];rows=[]
    for lam in [.25,.1,.025,.01]:
        B=2/sqrt(lam*(1-exp(-4*lam)));L=max(8,int(B)+2)
        while 4*taylor_tail(B,L)>1e-13:L+=1
        a=jacobi(lam,times,L);b=jacobi(lam,times,L+12)
        err=float(np.max(abs(a-b)));assert err<2e-12
        row=dict(lam=lam,L=L,refined_L=L+12,analytic_truncation_bound=4*taylor_tail(B,L),
           refinement_difference=err,times=times,normalized_memory_values=b.tolist(),
           small_lambda_values=(1/np.cosh(2*np.array(times))**2).tolist())
        if lam in [.25,.1]:
            lo=quadrature(lam,times,128);hi=quadrature(lam,times,192)
            diff=float(np.max(abs(b-hi['values'])))
            assert max(diff,abs(hi['mass']-1),abs(hi['second_moment']-2/lam*(1-exp(-4*lam))))<2e-10
            row.update(mehler_coarse=lo,mehler_fine=hi,independent_difference=diff)
        rows.append(row)
    return rows

def main():
    target=OUT/'kernel-checks.json'
    if target.exists():raise FileExistsError(target)
    record=dict(status='passed',tagged_car=tagged_car(),chord_identities=chord_identities(),
      mehler_density_limit=mehler_limit(),curves=curves(),
      environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__,mpmath=mp.__version__),
      error_scope='Exact isolated-tag CAR and chord polynomials. Off-diagonal Mehler limit at 65 digits; equality is an almost-everywhere density limit, not pointwise on the diagonal. Jacobi truncation bounded analytically; roundoff and quadrature checked by refinement, not interval-certified.',
      input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['code/check_force_kernel.py','code/complex_syk_dictionary.py','code/singlet_coefficient.py','code/singlet_otoc.py']})
    target.write_text(json.dumps(record,indent=2)+'\n')
    print('tagged CAR',len(record['tagged_car']['rows']),'nonzero',record['tagged_car']['nonzero'],
      'exact polynomials',len(record['chord_identities']),'Mehler limits',len(record['mehler_density_limit']))
    for r in record['curves']:print(r['lam'],r['normalized_memory_values'])

if __name__=='__main__':main()
