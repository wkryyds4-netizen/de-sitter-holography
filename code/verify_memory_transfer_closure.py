"""Verify full reverse orthogonality and exact source-variance normalization."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,platform
import numpy as np
from memory_charged_transfer import system
from memory_wick_rank import parameters

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-charged-transfer'


def check(N,p,m):
    s=system(N,p,m);T=s['T'];S=T.conj().T
    dp,dm,n,alpha=s['dp'],s['dm'],s['n'],s['alpha']
    cross=np.zeros((dp*dp,dp*dp),complex)
    full=np.zeros_like(cross)
    for u,v in zip(s['up'],s['low']):
        Lu=np.kron(u,np.eye(dp));Lv=np.kron(v,np.eye(dm))
        Su=(S@(u.reshape(-1)/np.sqrt(dp))).reshape(dm,dm)*np.sqrt(dm)
        V=v-alpha*Su
        source=np.kron(V,np.eye(dm))@S
        defect=S@Lu-np.kron(Su,np.eye(dm))@S
        cross+=source.conj().T@defect
        B=Lv@S-alpha*S@Lu
        full+=B.conj().T@B
    nu=p*p*(m-p)*N*(N+1)/(m**3*n*(n+p-1))
    eta=alpha**2*p/(m*n*n*(n+p-1))
    prediction=nu*(np.eye(dp*dp)-s['Qp']/(m*n))+eta*s['Qp']
    errors=dict(complex_orthogonality=float(np.linalg.norm(cross)),
                complete_reverse_gram=float(np.linalg.norm(full-prediction)))
    assert max(errors.values())<3e-10
    return dict(N=N,p=p,m=m,nu0=nu,errors=errors)


if __name__=='__main__':
    dest=OUT/'closure-checks.json'
    if dest.exists():raise FileExistsError(dest)
    rows=[check(*a) for a in [(3,1,2),(4,1,2),(4,2,3),(5,1,2),(5,2,3),(5,2,4)]]
    exact=0
    for N in (16,24,32,40):
        for p in range(1,N//8+1):
            for m in range(N//4,3*N//4+1):
                n=N-m+1;alpha=F(m-p,m)
                e,b,d2,u,v,c,B,V1=parameters(N,p,m-1)
                second=(v-1+d2*c)/b
                source=alpha*F(n+p,n)*second/F(m*n)**2
                target=F(p*p*(m-p)*N*(N+1),m**3*n*(n+p-1))
                assert source==target
                exact+=1
    names=['code/verify_memory_transfer_closure.py','code/memory_charged_transfer.py',
           'code/memory_wick_rank.py','code/cubic_linear_dynamics.py',
           'code/cubic_linear_channels.py','code/cubic_linear_limit.py',
           'code/complex_syk_dictionary.py']
    result=dict(systems=rows,exact_source_variance_cases=exact,
                scope='Complete complex source-defect orthogonality, full reverse Gram, and rational Casimir second-moment reduction. No inference from finite fits to a uniform propagation rate.',
                environment=dict(python=platform.python_version(),numpy=np.__version__),
                input_sha256={x:hashlib.sha256((ROOT/x).read_bytes()).hexdigest() for x in names})
    with dest.open('x',encoding='utf-8') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='passed',complete_cross_systems=len(rows),exact_variance_cases=exact)))
