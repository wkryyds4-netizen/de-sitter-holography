"""Independent reverse-source and polar superoperator checks, without sampling."""
from pathlib import Path
import hashlib,json,platform,time
import numpy as np
import sympy as sp
from memory_charged_transfer import system, scalar_on_ranks

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-charged-transfer'


def checks(N,p,m):
    s=system(N,p,m)
    T=s['T'];S=T.conj().T;Qp=s['Qp'];Qm=s['Qm']
    dp,dm,n,alpha=s['dp'],s['dm'],s['n'],s['alpha']
    ev,V=np.linalg.eigh(T.conj().T@T)
    keep=ev>1e-10
    inv=(V[:,keep]/np.sqrt(ev[keep]))@V[:,keep].conj().T
    U=T@inv
    Pm=V[:,keep]@V[:,keep].conj().T
    centered=np.zeros((dp*dp,dp*dp),complex)
    polar=np.zeros((dm*dm,dm*dm),complex)
    adjoint=np.zeros((dp*dp,dp*dp),complex)
    for u,v in zip(s['up'],s['low']):
        Lu=np.kron(u,np.eye(dp));Lv=np.kron(v,np.eye(dm))
        vecv=v.reshape(-1)/np.sqrt(dm);vecu=u.reshape(-1)/np.sqrt(dp)
        source=(vecv-alpha*S@vecu).reshape(dm,dm)*np.sqrt(dm)
        B=Lv@S-alpha*S@Lu-np.kron(source,np.eye(dm))@S
        centered+=B.conj().T@B
        A=U@Lv-Lu@U
        polar+=A.conj().T@A
        R=Lv@U.conj().T-U.conj().T@Lu
        adjoint+=R.conj().T@R
    coefficient=alpha**2*p/(m*n*n*(n+p-1))
    err=float(np.linalg.norm(centered-coefficient*Qp))
    assert err<3e-10
    errors=dict(centered_reverse_gram=err,
                polar_partial_isometry=float(np.linalg.norm(U.conj().T@U-Pm)),
                polar_casimir_intertwining=float(np.linalg.norm(Qp@U-U@Qm)),
                polar_identity=float(np.linalg.norm(U@(np.eye(dm).reshape(-1)/np.sqrt(dm))-np.eye(dp).reshape(-1)/np.sqrt(dp))))
    assert max(errors.values())<3e-10
    return dict(N=N,p=p,m=m,C_upper=1,common_HS_dimension=int(keep.sum()),errors=errors,
                centered_reverse_coefficient=coefficient,
                polar_forward_ranks=scalar_on_ranks(polar,Qm,N),
                polar_reverse_ranks=scalar_on_ranks(adjoint,Qp,N))


def symbolic():
    m,n,p,r=sp.symbols('m n p r',positive=True)
    alpha=(m-p)/m
    sq=lambda x:(m-x)*(n-x)/(m*n)
    x=r*(m+n-r)
    rhs=p*(alpha*(m-n)+2*r*(n-r)/m-p*x/m**2)/(m*n)
    assert sp.factor(alpha**2*sq(r)-sq(r+p)-rhs)==0
    return dict(exact_squared_amplitude_identity=True,
                alpha_squared_Sr_squared_minus_Srp_squared=str(sp.factor(rhs)))


if __name__=='__main__':
    start=time.time();dest=OUT/'polar-and-reverse-checks.json'
    if dest.exists():raise FileExistsError(dest)
    rows=[]
    for args in [(3,1,2),(4,1,2),(4,2,3),(5,1,2),(5,2,3),(5,2,4)]:
        rows.append(checks(*args))
        print(json.dumps(dict(system=args,status='passed')),flush=True)
    names=['code/check_memory_polar_transfer.py','code/memory_charged_transfer.py','code/complex_syk_dictionary.py']
    result=dict(systems=rows,symbolic=symbolic(),
                scope='Complete centered reverse-defect Grams and polar partial-isometry maps at finite N; exact amplitude polynomial identity. No Monte Carlo or truncated propagator approximation.',
                environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=sp.__version__),
                input_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names},
                elapsed_seconds=time.time()-start)
    with dest.open('x',encoding='utf-8') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='passed',systems=len(rows),elapsed_seconds=result['elapsed_seconds'])))
