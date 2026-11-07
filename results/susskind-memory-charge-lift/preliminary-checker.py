"""Exact obstruction to scalar product filters and charge-removal norm identity."""
from pathlib import Path
from fractions import Fraction
from math import comb
import hashlib,json,platform,time
import numpy as np
import sympy as sp
from complex_syk_dictionary import fermions
from memory_quadratic_contraction import Contraction,inner_exact
from memory_rank_blocks import casimir_matrix,projectors,exact_integer

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-charge-lift'

def family_value(x,p):
    value=1
    for a in range(p):value*=((x>>(2*a))&1)-((x>>(2*a+1))&1)
    return value

def family_A(x,p):
    t=[((x>>(2*a))&1)+((x>>(2*a+1))&1) for a in range(p)]
    b=[int(z==1) for z in t]
    return sum(t[a]*np.prod([b[j] for j in range(p) if j!=a],dtype=int) for a in range(p))

def diagonal_casimir(N,states,values):
    lookup=dict(zip(states,values));result=[]
    for x in states:
        result.append(2*sum(lookup[x]-lookup[x^(1<<i)^(1<<j)]
                    for i in range(N) if x&(1<<i)
                    for j in range(N) if not x&(1<<j)))
    return result

def main():
    OUT.mkdir(exist_ok=True);target=OUT/'checks.json'
    if target.exists():raise FileExistsError(target)
    start=time.time();rows=[]
    for N in range(4,13,2):
        m=N//2;states=[x for x in range(1<<N) if x.bit_count()==m];mask=(1<<N)-1
        for p in range(1,m):
            h=[family_value(x,p) for x in states]
            a=[int(family_A(x,p)) for x in states]
            direct=[sum((family_value(x^(1<<i),p)-family_value(x,p))**2
                        for i in range(N) if x&(1<<i)) for x in states]
            assert a==direct
            assert diagonal_casimir(N,states,h)==[2*p*(N+1-p)*v for v in h]
            odd=[]
            for x,v,av in zip(states,h,a):
                assert family_value(mask^x,p)==(-1)**p*v
                predicted=sum((((x>>(2*j))&1)+((x>>(2*j+1))&1)-1)*
                    np.prod([int((((x>>(2*k))&1)+((x>>(2*k+1))&1))==1)
                        for k in range(p) if k!=j],dtype=int) for j in range(p))
                assert av-int(family_A(mask^x,p))==2*predicted
                odd.append(int(predicted))
            norm=Fraction(sum(v*v for v in odd),len(states));assert norm>0
            rows.append(dict(N=N,p=p,m=m,pure_input_rank=p,
                        odd_A_norm_squared=str(norm),states=len(states)))
    o=Contraction(6,2,3)
    def mat(m):return np.diag([family_value(x,2) for x in range(64) if x.bit_count()==m]).astype(complex)
    h,k=mat(3),mat(2)
    assert np.array_equal(casimir_matrix(6,2)@k.ravel(),20*k.ravel())
    assert np.array_equal(casimir_matrix(6,3)@h.ravel(),20*h.ravel())
    pr=projectors(casimir_matrix(6,3),6,3)
    pn,pd=pr[1]
    h21=(pn@(h@h).ravel()).reshape(o.d,o.d)
    a1=(pn@o.A(h,k).ravel()).reshape(o.d,o.d)/pd
    q=o.contract(h,k);q1=(pn@q.ravel()).reshape(o.d,o.d)/pd
    D=np.diag([(x&15).bit_count()-2 for x in range(64) if x.bit_count()==3])
    assert not np.any(h21)
    assert np.array_equal(a1,D) and np.array_equal(q1,99*D)
    assert np.array_equal(q,o.direct_contract(h,k))
    assert np.array_equal(o.core(h,k),o.raw_word_core(h,k))
    exact_integer(8*q)
    # Transpose-particle-hole map on the full charge-three matrix algebra.
    states=[x for x in range(64) if x.bit_count()==3];pos={x:i for i,x in enumerate(states)}
    star=np.zeros((20,20),dtype=int)
    for j,x in enumerate(states):
        occupied=[i for i in range(6) if x&(1<<i)]
        star[pos[63^x],j]=(-1)**sum(i-a for a,i in enumerate(occupied))
    tau=np.kron(star,star)@np.eye(400,dtype=int).reshape(20,20,400).transpose(1,0,2).reshape(400,400)
    assert np.array_equal(tau@tau,np.eye(400,dtype=int))
    for r,(pn,pd) in enumerate(pr):assert np.array_equal(tau@pn,(-1)**r*pn)
    # Primitive traceless tensor and exact HS singular value of Rem.
    cs=fermions(6);lift=[]
    for m in [2,3,4]:
        obj=o if m==3 else Contraction(6,2,m)
        ix=[x for x in range(64) if x.bit_count()==m]
        iz=[x for x in range(64) if x.bit_count()==m-1]
        for r in range(min(m-1,6-m)+1):
            T=np.eye(64,dtype=complex)
            for j in range(r):T=T@cs[j].conj().T
            for j in reversed(range(r,2*r)):T=T@cs[j]
            plus=T[np.ix_(ix,ix)];minus=T[np.ix_(iz,iz)]
            rem=obj.rem(minus)
            assert np.array_equal(rem,(m-r)*plus)
            ell=Fraction(m*(m-r)*(7-m-r),7-m)
            assert inner_exact(rem,rem)==ell*inner_exact(minus,minus)
            lift.append(dict(N=6,m=m,rank=r,squared_singular_value=str(ell)))
    names=['code/check_memory_charge_lift.py','code/probe_memory_half_filling.py',
           'code/memory_quadratic_contraction.py','code/memory_rank_blocks.py',
           'code/complex_syk_dictionary.py','code/singlet_wick_moments.py']
    result=dict(status='passed',family_cases=rows,charge_lift_checks=lift,
        witness=dict(N=6,p=2,m=3,w=1,input_rank=2,output_rank=1,
            h='(n1-n2)(n3-n4)',D='n1+n2+n3+n4-2I',P1_h_squared='0',
            P1_A='D',P1_Q_core='99D',A_norm_squared=str(inner_exact(a1,a1)),
            Q_norm_squared=str(inner_exact(q1,q1)),
            independent_full_Gaussian_polarization_coordinates=len(o.rs)),
        full_operator_particle_hole_rank_projectors_verified=4,
        scope='Counterexample to general scalar filtering by P_r H^2; not a double-scaled stochastic norm or a new dynamical mode.',
        elapsed_seconds=time.time()-start,
        environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=sp.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='passed',family_cases=len(rows),lift_checks=len(lift),witness=result['witness'])))

if __name__=='__main__':main()
