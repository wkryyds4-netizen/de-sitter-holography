"""Finite full-Hermite checks of creation commutators and rank transitions."""
from pathlib import Path
from math import comb,factorial
from itertools import product
import hashlib,json,platform
import numpy as np
from complex_syk_dictionary import fermions,creators
from memory_wick_rank import parameters,mu
from cubic_linear_dynamics import crossing

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-wick-rank'


def moment_norm(vec):
    return sum(np.prod([factorial(a) for a in key])*np.vdot(x,x).real/x.shape[0] for key,x in vec.items())


def main():
    dest=OUT/'creation-commutator-checks.json'
    if dest.exists():raise FileExistsError(dest)
    rows=[]
    for N,p,m in [(3,1,2),(4,1,2)]:
        cs=fermions(N);bs=creators(cs,p);ids=[x for x in range(2**N) if x.bit_count()==m]
        d=len(ids);C=comb(m,p)*comb(N-m+p,p);L=len(bs);basis=[]
        T=[[ (a@b.conj().T)[np.ix_(ids,ids)]/np.sqrt(C) for b in bs] for a in bs]
        for i in range(L):
            basis.append(T[i][i])
            for j in range(i):basis.extend([(T[i][j]+T[j][i])/np.sqrt(2),1j*(T[i][j]-T[j][i])/np.sqrt(2)])
        B=np.array(basis);gdim=len(B);I=np.eye(d)
        assert np.linalg.norm(sum(x@x for x in B)-I)<1e-12
        E=[[ (a.conj().T@b)[np.ix_(ids,ids)] for b in cs] for a in cs];cas=np.zeros((d*d,d*d),complex)
        for i,j in product(range(N),repeat=2):
            ca=np.kron(E[i][j],I)-np.kron(I,E[i][j].T)
            cb=np.kron(E[j][i],I)-np.kron(I,E[j][i].T)
            cas+=ca@cb
        ev,U=np.linalg.eigh(cas);R=min(m,N-m)
        pr=[U[:,abs(ev-2*r*(N+1-r))<1e-8] for r in range(R+1)]
        def project(x,r):
            v=x.reshape(-1);return (pr[r]@(pr[r].conj().T@v)).reshape(d,d)
        def rankop(x):return sum(r*project(x,r) for r in range(R+1))
        def defect(vec,k):return {a:k*p*x-rankop(x) for a,x in vec.items()}
        def step(vec,creation_only=False):
            out={}
            for key,x in vec.items():
                for j,b in enumerate(B):
                    up=list(key);up[j]+=1;up=tuple(up)
                    out[up]=out.get(up,0)+b@x
                    if key[j] and not creation_only:
                        down=list(key);down[j]-=1;down=tuple(down)
                        out[down]=out.get(down,0)+key[j]*b@x
            return out
        vec={(0,)*gdim:I.astype(complex)}
        cases=[]
        for power in range(5):
            for k in range(power%2,power+1,2):
                psi={a:x for a,x in vec.items() if sum(a)==k}
                if not psi:continue
                for x in psi.values():
                    for r in range(k*p+1,R+1):assert np.linalg.norm(project(x,r))<1e-10
                created=step(psi,True);left=defect(created,k+1);right=step(defect(psi,k),True)
                comm={a:left[a]-right[a] for a in left};value=moment_norm(comm)
                raw=0.
                for alpha,x in psi.items():
                    weight=np.prod([factorial(v) for v in alpha])
                    for r in range(R+1):
                        xr=project(x,r)
                        for b in B:
                            for t in range(R+1):
                                y=project(b@xr,t)
                                raw+=weight*(p+r-t)**2*np.vdot(y,y).real/d
                bound=(k+1)*raw
                assert value<=bound+2e-10*max(1,bound)
                # Deterministic rank cross terms of the summed insertion map.
                crossmax=0.
                for x in list(psi.values())[:3]:
                    for r,rp,t in product(range(R+1),repeat=3):
                        if r==rp:continue
                        cross=sum(np.vdot(project(b@project(x,r),t),project(b@project(x,rp),t)) for b in B)
                        crossmax=max(crossmax,float(abs(cross)))
                assert crossmax<1e-10
                cases.append(dict(power=power,chaos_degree=k,coefficient_matrices=len(psi),
                                  commutator_norm_squared=float(value),raw_tensor_bound=float(bound),
                                  rank_cross_max=crossmax))
            if power<4:vec=step(vec)
        rows.append(dict(N=N,p=p,m=m,C_normalized=1,Gaussian_coordinates=gdim,cases=cases))
    # Exact one-step second moments, including the large-degree branch.
    exact=0
    for N,p,m in [(16,1,4),(16,2,8),(24,2,6),(24,3,12),(32,2,24)]:
        R=min(m,N-m);dim=comb(N,m)
        dims=[comb(N,s)**2-(comb(N,s-1)**2 if s else 0) for s in range(R+1)]
        from fractions import Fraction as F
        lam=F(p*p,N)
        for r in range(R+1):
            transition=[F(dims[t],dim*dim)*sum(dims[s]*crossing(N,m,r,s)*mu(N,p,m,s)*crossing(N,m,t,s)
                        for s in range(R+1)) for t in range(R+1)]
            assert all(v==0 for t,v in enumerate(transition) if abs(t-r)>p)
            k=(r+p-1)//p
            value=sum(v*(p+r-t)**2 for t,v in enumerate(transition))
            assert value<=1000000*(lam*(k+1)+lam**2*(k+1)**2)
            exact+=1
    result=dict(systems=rows,exact_step_cases=exact,
                scope='Complete Hermite coefficient recurrences through H^4, including lower chaoses; creation-block and deterministic rank-cross checks. No truncated exponential is used as proof of the real-time bound.',
                environment=dict(python=platform.python_version(),numpy=np.__version__),
                input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in
                              ['code/check_memory_chaos_creation.py','code/memory_wick_rank.py','code/cubic_linear_dynamics.py','code/complex_syk_dictionary.py']})
    dest.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='passed',Hermite_cases=sum(len(x['cases']) for x in rows),exact_step_cases=exact)))


if __name__=='__main__':main()
