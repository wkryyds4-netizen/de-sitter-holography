"""Complete CAR coefficient checks for adjacent-charge Gaussian transfer.

Matrices act on row-vectorized operators with normalized HS coordinates.
No external scripts, sample covariance estimates or rank truncations are used.
"""
from pathlib import Path
from math import comb
import hashlib
import json
import platform
import time
import numpy as np
import sympy as sp
from complex_syk_dictionary import fermions, creators

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'results/susskind-memory-charged-transfer'


def system(N, p, m):
    cs = fermions(N)
    bs = creators(cs, p)
    ip = [x for x in range(2**N) if x.bit_count() == m]
    im = [x for x in range(2**N) if x.bit_count() == m-1]
    dp, dm = len(ip), len(im)
    Cbase = comb(m, p)*comb(N-m+p, p)
    all_basis = []
    for i, a in enumerate(bs):
        all_basis.append(a @ a.conj().T)
        for b in bs[:i]:
            z = a @ b.conj().T
            all_basis.extend([(z+z.conj().T)/np.sqrt(2),
                              1j*(z-z.conj().T)/np.sqrt(2)])
    up = [b[np.ix_(ip, ip)]/np.sqrt(Cbase) for b in all_basis]
    low = [b[np.ix_(im, im)]/np.sqrt(Cbase) for b in all_basis]
    cc = [c[np.ix_(im, ip)] for c in cs]
    T = sum(np.kron(c.conj().T, c.T) for c in cc)/m*np.sqrt(dm/dp)
    def casimir(ids):
        d = len(ids)
        E = [[(a.conj().T@b)[np.ix_(ids, ids)] for b in cs] for a in cs]
        Q = np.zeros((d*d, d*d), complex)
        for i in range(N):
            for j in range(N):
                A = np.kron(E[i][j], np.eye(d))-np.kron(np.eye(d), E[i][j].T)
                B = np.kron(E[j][i], np.eye(d))-np.kron(np.eye(d), E[j][i].T)
                Q += A@B/2
        return Q
    Qp, Qm = casimir(ip), casimir(im)
    return dict(N=N, p=p, m=m, n=N-m+1, alpha=(m-p)/m,
                dp=dp, dm=dm, up=up, low=low, T=T, Qp=Qp, Qm=Qm)


def scalar_on_ranks(A, Q, N):
    ev, U = np.linalg.eigh(Q)
    rows = []
    for r in range(N//2+1):
        basis = U[:, abs(ev-r*(N+1-r)) < 1e-8]
        if not basis.shape[1]:
            continue
        block = basis.conj().T@A@basis
        value = np.trace(block).real/basis.shape[1]
        error = np.linalg.norm(A@basis-value*basis)
        assert error < 2e-10
        rows.append(dict(rank=r, dimension=basis.shape[1], eigenvalue=float(value),
                         scalar_block_error=float(error)))
    return rows


def check_system(N, p, m):
    s = system(N, p, m)
    T, Qp, Qm = s['T'], s['Qp'], s['Qm']
    dp, dm, n, alpha = s['dp'], s['dm'], s['n'], s['alpha']
    Tp = T.conj().T
    forward = np.zeros((dm*dm, dm*dm), complex)
    reverse = np.zeros((dp*dp, dp*dp), complex)
    vacp = np.eye(dp).reshape(-1)/np.sqrt(dp)
    vacm = np.eye(dm).reshape(-1)/np.sqrt(dm)
    source_error = 0.
    coefficient_error = 0.
    raw_source = 0.
    for u, v in zip(s['up'], s['low']):
        Lu = np.kron(u, np.eye(dp)); Lv = np.kron(v, np.eye(dm))
        A = T@Lv-alpha*Lu@T
        B = Lv@Tp-alpha*Tp@Lu
        forward += A.conj().T@A
        reverse += B.conj().T@B
        source = B@vacp
        exact_source = Qm@(v.reshape(-1)/np.sqrt(dm))/(m*n)
        source_error = max(source_error, float(np.linalg.norm(source-exact_source)))
        raw_source += float(np.vdot(source, source).real)
        coefficient_error = max(coefficient_error, float(np.linalg.norm(T@Lv@vacm-alpha*Lu@vacp)))
    gamma = p*(m-p)/(m**3*(m-1)*n)
    errors = dict(
        lower_channel_gram=float(np.linalg.norm(Tp@T-np.eye(dm*dm)+Qm/(m*n))),
        upper_channel_gram=float(np.linalg.norm(T@Tp-np.eye(dp*dp)+Qp/(m*n))),
        casimir_intertwining=float(np.linalg.norm(Qp@T-T@Qm)),
        forward_gram=float(np.linalg.norm(forward-gamma*Qm)),
        vacuum_source=source_error, Hamiltonian_intertwining=coefficient_error)
    assert max(errors.values()) < 2e-10
    return dict(N=N,p=p,m=m,C_upper=1,dimensions=[dp,dm],
                real_Gaussian_coordinates=len(s['up']),errors=errors,
                forward_ranks=scalar_on_ranks(forward,Qm,N),
                reverse_ranks=scalar_on_ranks(reverse,Qp,N),
                reverse_vacuum_norm_squared=raw_source,
                gamma=gamma)


def symbolic_checks():
    m,n,p,C=sp.symbols('m n p C', positive=True)
    alpha=(m-p)/m
    Cm=C*alpha*(n+p)/n
    c=p*n*(m+n)/(n+p)
    kappa=C*(m-p)*(m-p-1)/(m*(m-1))
    gamma=C*p*(m-p)/(m**3*(m-1)*n)
    assert sp.factor(Cm*(1-c/(m*n))-alpha**2*C)==0
    assert sp.factor((alpha**2*C-kappa)/(m*n)-gamma)==0
    return dict(exact_identities=2,forward_gram_coefficient=str(gamma))


if __name__ == '__main__':
    started=time.time()
    OUT.mkdir(parents=True,exist_ok=True)
    dest=OUT/'transfer-checks.json'
    if dest.exists(): raise FileExistsError(dest)
    rows=[]
    for args in [(3,1,2),(4,1,2),(4,2,3),(5,1,2),(5,2,3),(5,2,4)]:
        row=check_system(*args);rows.append(row)
        print(json.dumps(dict(system=args,status='passed',max_error=max(row['errors'].values()))),flush=True)
    names=['code/memory_charged_transfer.py','code/complex_syk_dictionary.py']
    result=dict(systems=rows,symbolic=symbolic_checks(),
                scope='Complete deterministic superoperator Grams over every Gaussian coordinate and every matrix input. Reverse scalar blocks are measured, not presumed zero. These finite checks do not prove a uniform memory-field bound.',
                environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=sp.__version__),
                input_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names},
                elapsed_seconds=time.time()-started)
    with dest.open('x',encoding='utf-8') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='passed',systems=len(rows),elapsed_seconds=result['elapsed_seconds'])))
