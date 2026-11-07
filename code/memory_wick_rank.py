"""All-degree rank-deficit bounds for the highest coupling chaos.

Fresh output only. Exact Racah inversion and symbolic identities are checked
against CAR channel composition, not finite-size fits of a claimed limit.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import hashlib
import json
import platform
import time
import numpy as np
import sympy as sp
import mpmath as mp
from cubic_linear_dynamics import crossing
from complex_syk_dictionary import fermions, creators

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'results/susskind-memory-wick-rank'


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def mu(N, p, m, s):
    return F(choose(m-s, p)*choose(N-m+p-s, p),
             choose(m, p)*choose(N-m+p, p))


def parameters(N, p, m):
    L = N-m+p
    e = F(m*(N-m)*(N+1), N)
    b = F((N-1)*(N-2), m*(m-1)*(N-m)*(N-m-1)*N*(N+1))
    d = F(2*(N-1), N)/e+b*N
    u, v = mu(N,p,m,1), mu(N,p,m,2)
    c = e*(1-u)
    B = d/b*(u-v)-2*u*c
    V1 = (v-1+d*c)/b-c*c
    return e,b,d,u,v,c,B,V1


def distribution(N, p, m, k):
    R=min(m,N-m); dim=comb(N,m)
    dims=[choose(N,s)**2-choose(N,s-1)**2 for s in range(R+1)]
    prob=[F(dims[r],dim**2)*sum(dims[s]*mu(N,p,m,s)**k*crossing(N,m,r,s)
          for s in range(R+1)) for r in range(R+1)]
    assert sum(prob)==1 and min(prob)>=0
    assert all(prob[r]==0 for r in range(k*p+1,R+1))
    return prob


def exact_checks():
    rows=[]; identities=0
    for N,p,m in [(5,1,2),(6,2,3),(8,2,3),(8,2,4),(10,2,5),(12,2,4)]:
        e,b,d,u,v,c,B,V1=parameters(N,p,m)
        R=min(m,N-m); dim=comb(N,m)
        dims=[choose(N,s)**2-choose(N,s-1)**2 for s in range(R+1)]
        transition=[]
        for a in range(R+1):
            line=[F(dims[r],dim**2)*sum(dims[s]*crossing(N,m,a,s)*mu(N,p,m,s)*crossing(N,m,r,s)
                  for s in range(R+1)) for r in range(R+1)]
            assert min(line)>=0 and sum(line)==1
            x=a*(N+1-a)
            mean=sum(pr*r*(N+1-r) for r,pr in enumerate(line))
            second=sum(pr*(r*(N+1-r))**2 for r,pr in enumerate(line))
            assert mean==u*x+c
            assert second-mean*mean==(v-u*u)*x*x+B*x+V1
            transition.append(line);identities+=4
        prev=[F(1)]+[F(0)]*R
        for k in range(1,9):
            pr=distribution(N,p,m,k)
            assert pr==[sum(prev[a]*transition[a][r] for a in range(R+1)) for r in range(R+1)]
            mean=e*(1-u**k);second=(v**k-1+d*mean)/b
            assert mean==sum(pr[r]*r*(N+1-r) for r in range(R+1))
            assert second==sum(pr[r]*(r*(N+1-r))**2 for r in range(R+1))
            prev=pr;identities+=3
        rows.append(dict(N=N,p=p,m=m,k_checked=list(range(1,9)),last_probabilities=list(map(str,prev))))
    return dict(rows=rows,exact_equalities=identities)


def symbolic_checks():
    N,p,m,x=sp.symbols('N p m x');L=N-m+p
    e=m*(N-m)*(N+1)/N
    b=(N-1)*(N-2)/(m*(m-1)*(N-m)*(N-m-1)*N*(N+1))
    d=2*(N-1)/(N*e)+b*N
    u=(m-p)*(N-m)/(m*L)
    v=u*(m-p-1)*(N-m-1)/((m-1)*(L-1));c=e*(1-u)
    P=N**3-3*N**2*m+N**2*p+4*N*m**2-4*N*m*p+2*N*p-2*N-2*m**3+2*m**2*p-2*m*p+2*m
    Bref=p*(N-m)*(m-p)*P/(m*(m-1)*L**2*(L-1))
    Vref=p**2*(N-m)**2*(N+1)*(m-p+1)/(L**2*(L-1))
    checks=[sp.cancel(c-p*(N-m)*(N+1)/L),
            sp.cancel(d/b*(u-v)-2*u*c-Bref),
            sp.cancel((v-1+d*c)/b-c*c-Vref)]
    assert checks==[0,0,0]
    conditional=(v*(1-d*x+b*x*x)-1+d*(u*x+c))/b-(u*x+c)**2
    assert sp.cancel(conditional-((v-u*u)*x*x+Bref*x+Vref))==0
    lam=sp.symbols('lam',positive=True);nn=p*p/lam
    limits=[]
    for k in range(1,7):
        mean=e*(1-u**k);second=(v**k-1+d*mean)/b
        top=k*p*(N+1-k*p)
        expressions=[(top-mean)/N,(top*top-2*top*mean+second)/N**2]
        vals=[]
        for expr in expressions:
            z=sp.cancel(expr.subs(m,N/2).subs(N,nn));num,den=z.as_numer_denom()
            a=sp.Poly(num,p);bb=sp.Poly(den,p)
            assert a.degree()<=bb.degree()
            vals.append(sp.factor(a.LC()/bb.LC()) if a.degree()==bb.degree() else sp.Integer(0))
        assert sp.simplify(vals[0]-lam*k*k)==0
        assert sp.simplify(vals[1]-(lam*k*k+lam**2*k**4))==0
        limits.append(dict(k=k,scaled_casimir_deficit_moments=list(map(str,vals))))
    return dict(general_rational_identities=4,B=str(Bref),V1=str(Vref),
                fixed_k_symbolic_checks=limits,
                scope='General-k limit follows the written induction, not these six k checks.')


def rational_bounds():
    count=0
    for N in range(8,81,4):
        for p in range(1,N//8+1):
            for m in sorted(set([N//4,N//2,3*N//4])):
                e,b,d,u,v,c,B,V1=parameters(N,p,m)
                assert 0<=v<=u*u<=1
                assert abs(B)<=102400*p and V1<=512*p*p*N
                assert c<=8*p*N and 1-u<=F(16*p,N)
                for k in [1,2,3,5,11,2*N//p+1]:
                    mean=e*(1-u**k);second=(v**k-1+d*mean)/b
                    var=second-mean*mean
                    assert var<=410112*k*k*p*p*N
                    if k*p<=N//8:
                        delta=k*p*(N+1-k*p)-mean
                        assert 0<=delta<=73*k*k*p*p
                    count+=1
    return count


def car_channel_checks():
    rows=[]
    for N,p,m in [(4,1,2),(4,2,2),(5,2,3)]:
        cs=fermions(N);bs=creators(cs,p)
        ids=[j for j in range(2**N) if j.bit_count()==m];dim=len(ids)
        C=choose(m,p)*choose(N-m+p,p)
        channel=np.zeros((dim**2,dim**2),complex)
        for a in bs:
            for bb in bs:
                T=(a@bb.conj().T)[np.ix_(ids,ids)]
                channel+=np.kron(T,T.conj())/C
        casimir=np.zeros_like(channel);I=np.eye(dim)
        E=[[ (a.conj().T@bb)[np.ix_(ids,ids)] for bb in cs] for a in cs]
        for i in range(N):
            for j in range(N):
                left=np.kron(E[i][j],I)-np.kron(I,E[i][j].T)
                right=np.kron(E[j][i],I)-np.kron(I,E[j][i].T)
                casimir+=left@right
        ev,U=np.linalg.eigh(casimir)
        projectors=[U[:,abs(ev-2*r*(N+1-r))<1e-8] for r in range(min(m,N-m)+1)]
        power=np.eye(dim**2);err=0.
        for k in range(1,6):
            power=power@channel
            covariance=power.reshape(dim,dim,dim,dim).transpose(0,2,1,3).reshape(dim**2,dim**2)
            probs=[float(np.trace(P.conj().T@covariance@P).real/dim) for P in projectors]
            exact=distribution(N,p,m,k)
            err=max(err,max(abs(a-float(bb)) for a,bb in zip(probs,exact)))
        assert err<2e-11
        rows.append(dict(N=N,p=p,m=m,product_lengths=[1,2,3,4,5],max_probability_error=err))
    return rows


def large_parameter_checks():
    rows=[]
    for precision in [80,120]:
        mp.mp.dps=precision
        values=[]
        for lam,p,k in [('0.25',16,2),('0.25',64,5),('0.5',128,3),('1',64,10),('0.5',32,100)]:
            N=int(mp.mpf(p)**2/mp.mpf(lam));m=N//2
            e,b,d,u,v,c,B,V1=[mp.mpf(z.numerator)/z.denominator for z in parameters(N,p,m)]
            mean=e*(1-u**k);second=(v**k-1+d*mean)/b;var=second-mean*mean
            top=k*p*(N+1-k*p)
            values.append((N,p,k,(top-mean)/N,(top*top-2*top*mean+second)/N**2,var/N**2))
        if precision==80:
            coarse=values
        else:
            for old,new in zip(coarse,values):
                err=max(abs(a-bb) for a,bb in zip(old[3:],new[3:]))
                assert err<mp.mpf('1e-60')
                N,p,k,a,bb,var=new
                rows.append(dict(N=N,p=p,m=N//2,k=k,lambda_value=str(mp.mpf(p)**2/N),
                                 scaled_casimir_deficit_mean=str(a),scaled_casimir_deficit_second=str(bb),
                                 scaled_casimir_variance=str(var),precision_difference=str(err),
                                 interpretation='Fixed-k DS check only when k is held fixed; large-k row tests stable exact moment evaluation.'))
    return rows


if __name__=='__main__':
    OUT.mkdir(parents=True,exist_ok=True);dest=OUT/'rank-chain-checks.json'
    if dest.exists():raise FileExistsError(dest)
    start=time.time()
    result=dict(exact=exact_checks(),symbolic=symbolic_checks(),
                rational_bound_instances=rational_bounds(),CAR_channels=car_channel_checks(),
                high_precision=large_parameter_checks(),
                environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=sp.__version__,mpmath=mp.__version__),
                scope='Independent-copy rank law and a positive highest-coupling-chaos bound. Not a finite-time bound for the full ordinary propagator.',
                input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in
                              ['code/memory_wick_rank.py','code/cubic_linear_dynamics.py','code/complex_syk_dictionary.py']})
    result['elapsed_seconds']=time.time()-start
    dest.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='passed',exact_equalities=result['exact']['exact_equalities'],
                         rational_bound_instances=result['rational_bound_instances'],CAR_systems=len(result['CAR_channels']))))
