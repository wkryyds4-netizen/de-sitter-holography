"""Full Gaussian-Hermite coefficient checks for the finite Jacobi bridge."""
from pathlib import Path
from math import factorial,sqrt
import hashlib,json,platform,time
import numpy as np
from memory_charged_transfer import system
from memory_wick_rank import mu,distribution

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-uniform-bound'


def add(a,b,scale=1):
    out={k:v.copy() for k,v in a.items()}
    for k,v in b.items():out[k]=out.get(k,0)+scale*v
    return out


def norm2(v):
    return sum(float(np.prod([factorial(x) for x in k]))*np.vdot(x,x).real/x.shape[0] for k,x in v.items())


def step(v,B,mode='full'):
    out={}
    for key,x in v.items():
        for i,b in enumerate(B):
            bx=b@x
            if mode!='annihilation':
                up=list(key);up[i]+=1;up=tuple(up)
                out[up]=out.get(up,0)+bx
            if key[i] and mode!='creation':
                down=list(key);down[i]-=1;down=tuple(down)
                out[down]=out.get(down,0)+key[i]*bx
    return out


def raw_product(a,b):
    out={}
    for ka,va in a.items():
        for kb,vb in b.items():
            key=tuple(x+y for x,y in zip(ka,kb))
            out[key]=out.get(key,0)+va@vb
    return out


def check(N,p,m,max_n):
    s=system(N,p,m);B=s['up'];d=s['dp'];key=(0,)*len(B)
    Z=[{key:np.eye(d,dtype=complex)}]
    for n in range(max_n+1):Z.append(step(Z[-1],B,'creation'))
    Qprev={};Q=Z[0];rows=[]
    lipschitz=p/(m-p+1)+p/(N-m+1)
    for n in range(max_n+1):
        E=add(Z[n],Q,-1)
        top=max((np.linalg.norm(x) for k,x in E.items() if sum(k)>n-2),default=0.)
        assert top<3e-11
        if n:
            a=sum(float(mu(N,p,m,j*p)) for j in range(n))
            AZ=step(Z[n],B,'annihilation')
            direct={}
            for j in range(n):
                Phi={k:sum(b@x@b for b in B) for k,x in Z[j].items()}
                direct=add(direct,raw_product(Phi,Z[n-1-j]))
            identity=sqrt(max(0,norm2(add(AZ,direct,-1))))
            assert identity<2e-9
            residual=add(AZ,Z[n-1],-a)
            observed=sqrt(max(0,norm2(residual)))
            total=0.
            for j in range(n):
                pr=distribution(N,p,m,j)
                total+=sqrt(float(sum(value*(j*p-r)**2 for r,value in enumerate(pr))))
            bound=lipschitz*sqrt(factorial(n-1))*total
            assert observed<=bound+2e-9
            rows.append(dict(n=n,chaos_matrices=len(Z[n]),a_n=a,
                annihilation_identity_error=identity,residual_norm=observed,
                independent_tensor_bound=bound,polynomial_bridge_error=sqrt(max(0,norm2(E))),
                cancelled_top_degree_error=float(top)))
        if n<max_n:
            an=sum(float(mu(N,p,m,j*p)) for j in range(n))
            Qnext=add(step(Q,B),Qprev,-an)
            Qprev,Q=Q,Qnext
    return dict(N=N,p=p,m=m,C=1,Gaussian_coordinates=len(B),max_degree=max_n,cases=rows)


if __name__=='__main__':
    start=time.time();OUT.mkdir(parents=True,exist_ok=True)
    dest=OUT/'jacobi-recurrence-checks.json'
    if dest.exists():raise FileExistsError(dest)
    rows=[]
    for args in [(3,1,2,5),(4,1,2,4),(4,2,3,3)]:
        row=check(*args);rows.append(row)
        print(json.dumps(dict(system=args,status='passed')),flush=True)
    names=['code/memory_jacobi_bridge.py','code/memory_charged_transfer.py',
           'code/memory_wick_rank.py','code/cubic_linear_dynamics.py',
           'code/complex_syk_dictionary.py']
    result=dict(systems=rows,
                scope='Complete Hermite coefficients, independent covariance-word reconstruction of annihilation, one factorial normalization and finite positive-rank tensor bounds. Small systems verify finite identities, not the uniform double-scaling theorem.',
                environment=dict(python=platform.python_version(),numpy=np.__version__),
                input_sha256={x:hashlib.sha256((ROOT/x).read_bytes()).hexdigest() for x in names},
                elapsed_seconds=time.time()-start)
    with dest.open('x',encoding='utf-8') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='passed',cases=sum(len(x['cases']) for x in rows),elapsed_seconds=result['elapsed_seconds'])))
