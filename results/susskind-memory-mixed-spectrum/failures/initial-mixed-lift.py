"""Upper-charge covariance and direct reverse-CAR lift, plus quartic marks."""
from pathlib import Path
import inspect,json,hashlib,platform
import numpy as np
from complex_syk_dictionary import fermions,creators
import singlet_filtered_bridge as old

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-mixed-spectrum'


def direct(N,p,m):
    cs=fermions(N);basis=creators(cs,p);lower=creators(cs,p-1)
    dim=1<<N;Q=np.diag([i.bit_count() for i in range(dim)]);I=np.eye(dim)
    ids=[i for i in range(dim) if i.bit_count()==m];d=m-p+1
    D=lambda X:Q@X-sum(c.conj().T@X@c for c in cs)
    Db=lambda X:(N*I-Q)@X-sum(c@X@c.conj().T for c in cs)
    M=lambda U,V,reverse: -sum((U@c-c@U)@(V@c.conj().T-c.conj().T@V) if reverse
                               else (U@c.conj().T-c.conj().T@U)@(V@c-c@V) for c in cs)
    restrict=lambda X:X[np.ix_(ids,ids)]
    rng=np.random.default_rng(110688+N*100+p*10+m)
    def sample(bs):
        raw=rng.integers(-2,3,(len(bs),len(bs)))+1j*rng.integers(-2,3,(len(bs),len(bs)))
        J=raw+raw.conj().T
        return sum(J[a,b]*aa@bb.conj().T for a,aa in enumerate(bs) for b,bb in enumerate(bs))
    rows=[]
    for trial in range(3):
        U=sample(lower);V=sample(basis);lift=(Q-(p-1)*I)@U
        assert np.max(abs(restrict(U)))>0 and np.max(abs(restrict(V)))>0
        checked=0;nonzero=0;maxentry=0
        def check(lhs,rhs):
            nonlocal checked,nonzero,maxentry
            lhs=restrict(lhs);rhs=restrict(rhs)
            for mat in [lhs,rhs]:
                assert np.all(mat.real==np.rint(mat.real)) and np.all(mat.imag==np.rint(mat.imag))
                maxentry=max(maxentry,float(np.max(abs(mat))));assert maxentry<2**48
            assert np.array_equal(lhs,rhs),(N,p,m,trial,checked)
            checked+=1;nonzero+=int(np.any(lhs))
        # Multiplication by d clears the scalar lift denominator.
        check(M(lift,V,True),(d+1)*M(U,V,True)-U@Db(V))
        check(M(V,lift,True),(d+1)*M(V,U,True)-Db(V)@U)
        check(M(lift,V,False),(d-1)*M(U,V,False)+U@D(V))
        check(M(V,lift,False),(d-1)*M(V,U,False)+D(V)@U)
        # General reverse superoperator identity; it is NOT a natural body mark.
        check(Db(lift),(d+1)*Db(U)-(N-m)*U)
        factors=[sample(basis) for _ in range(4)]
        for at in range(4):
            plain=I.copy();lifted=I.copy()
            for j,H in enumerate(factors):
                plain=plain@(U if j==at else H)
                lifted=lifted@(lift if j==at else H)
            check(4*p*lifted-D(lifted),(d-1)*((4*p-1)*plain-D(plain))+3*p*plain)
        assert nonzero==checked
        rows.append(dict(trial=trial,identities=checked,nonzero=nonzero,max_cleared_entry=maxentry))
    print('direct lift',N,p,m,'passed',flush=True)
    return dict(N=N,p=p,m=m,seed=110688+N*100+p*10+m,rows=rows)


def main():
    target=OUT/'charged-lift.json'
    if target.exists():raise FileExistsError(target)
    # Reuse the exact covariance tensors, explicitly extending their charge range.
    src=inspect.getsource(old.finite)
    assert src.count('for n in [m,m-1,m-2]:')==1
    src=src.replace('for n in [m,m-1,m-2]:','for n in [m+1]:')
    env=vars(old).copy();exec(compile(src,'<upper-charge covariance>','exec'),env)
    systems=[(4,2,2),(5,2,3),(6,2,3)]
    cov=[env['finite'](*x) for x in systems]
    rows=[direct(*x) for x in systems]
    names=['code/memory_mixed_lift.py','code/singlet_filtered_bridge.py',
           'code/complex_syk_dictionary.py','code/singlet_wick_moments.py',
           'code/cubic_covariance_bridge.py','code/cubic_linear_limit.py']
    data=dict(status='passed',upper_charge_covariance=cov,direct_CAR=rows,
        direct_identities=sum(z['identities'] for r in rows for z in r['rows']),
        covariance_entries=sum(z['entries'] for r in cov for z in r['covariance_blocks']),
        algorithm_extension='Existing exact rank-Casimir charged tensor routine, now at n=m+1; its prior normal product controls also retained.',
        modified_covariance_function_sha256=hashlib.sha256(src.encode()).hexdigest(),
        reverse_identity='d M_plus(gU,V)=(d+1)M_plus(U,V)-U Dbar(V), and d M_plus(V,gU)=(d+1)M_plus(V,U)-Dbar(V)U at m, g(Q)=(Q-p+1)/d.',
        environment=dict(python=platform.python_version(),numpy=np.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print('passed',data['direct_identities'],'direct identities',data['covariance_entries'],'covariance entries',flush=True)


if __name__=='__main__':main()
