"""Exact-rank neutral coefficient identities and retained finite slope errors."""
from pathlib import Path
from math import comb
from fractions import Fraction as F
import datetime
import hashlib
import json
import platform
import mpmath as mp
import numpy as np
from memory_charged_transfer import system
from memory_wick_rank import mu

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-two-component'


def hs(X):return float(np.linalg.norm(X)/np.sqrt(len(X)))


def finite_checks():
    records=[]
    for N,p,m in [(3,1,2),(4,1,2),(4,2,3),(5,2,3)]:
        sys=system(N,p,m);basis=np.array(sys['up']);d=sys['dp']
        def phi(X):return np.sum(basis@X@basis,axis=0)
        ev,U=np.linalg.eigh(sys['Qp'])
        ranks={r:U[:,abs(ev-r*(N+1-r))<1e-8] for r in range(min(m,N-m)+1)}
        def proj(X,r):
            W=ranks[r]
            return (W@(W.conj().T@X.reshape(-1))).reshape(d,d)
        covariance_error=max(hs(phi(B)-sum(float(mu(N,p,m,r))*proj(B,r)
                for r in ranks)) for B in basis)
        # C=1. bar_kappa is the exact annealed h-regression, evaluated over
        # every Gaussian coordinate, not from the fixed sample below.
        bar_kappa=float(sum(np.trace(B@phi(B)).real/d for B in basis))
        alpha=(m-p)/m;eps=1-alpha*alpha;k2=float(mu(N,p,m,2*p))
        for seed in [17,29,43]:
            h=np.einsum('a,aij->ij',np.random.default_rng(seed).normal(size=len(basis)),basis)
            Z=h@h-np.eye(d);BZ=phi(h)-bar_kappa*h
            raw1=alpha*phi(h);proxy1=alpha*BZ
            raw2=alpha*alpha*phi(h@h)-alpha*(h@phi(h)+phi(h)@h)
            proxy2=alpha*alpha*(phi(Z)-k2*Z)-alpha*(h@BZ+BZ@h)
            conserved2=(alpha*alpha*k2-2*alpha*bar_kappa)*(h@h)+alpha*alpha*(1-k2)*np.eye(d)
            e1=hs(raw1-proxy1-alpha*bar_kappa*h)
            e2=hs(raw2-proxy2-conserved2)
            assert max(e1,e2,covariance_error)<3e-12
            records.append(dict(N=N,p=p,m=m,seed=seed,C=1.,coordinates=len(basis),
                alpha=alpha,epsilon=eps,m_epsilon=m*eps,bar_kappa=bar_kappa,kappa_2p=k2,
                first_conserved_difference_error=e1,second_conserved_difference_error=e2,
                exact_rank_covariance_error=covariance_error))
    return records


def slope_checks():
    rows=[]
    for lam in [F(1,4),F(1,2),F(1)]:
        for p in [16,32,64,128]:
            N=int(p*p/lam)
            for x in [-2,0,2]:
                m=N//2+x*p
                for k in [1,2]:
                    for deficit in [0,1,2,3,6]:
                        value=p*(mu(N,p,m,k*p-deficit)-mu(N,p,m,k*p))
                        assert F(0)<=value<=12*lam*deficit
                        vals=[]
                        for precision in [80,120]:
                            mp.mp.dps=precision
                            actual=mp.mpf(value.numerator)/value.denominator
                            ll=mp.mpf(lam.numerator)/lam.denominator
                            expected=4*ll*mp.exp(-4*ll*k)*deficit
                            vals.append((str(actual),str(expected)))
                        precision_difference=abs(mp.mpf(vals[0][0])-mp.mpf(vals[1][0]))
                        assert precision_difference<mp.mpf('1e-70')
                        actual,expected=map(mp.mpf,vals[1])
                        rows.append(dict(lambda_value=str(lam),N=N,p=p,m=m,x=x,k=k,deficit=deficit,
                            finite_slope=str(actual),limiting_slope=str(expected),
                            signed_error=str(actual-expected),
                            relative_error=None if not deficit else str((actual-expected)/expected),
                            precision_difference=str(precision_difference),
                            dominating_bound=str(12*lam*deficit)))
    return rows


if __name__=='__main__':
    OUT.mkdir(parents=True,exist_ok=True);dest=OUT/'neutral-coefficient-checks.json'
    if dest.exists():raise FileExistsError(dest)
    finite=finite_checks();slopes=slope_checks()
    names=['code/check_neutral_memory_coefficients.py','code/memory_charged_transfer.py',
           'code/memory_wick_rank.py','code/cubic_linear_dynamics.py',
           'code/complex_syk_dictionary.py']
    result=dict(finite_CAR=finite,rank_slopes=slopes,
        created_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        environment=dict(python=platform.python_version(),numpy=np.__version__,mpmath=mp.__version__),
        scope='12 complete-coordinate CAR coefficient/covariance identities and 360 exact-rational slope inequalities with 80/120-digit comparisons. Finite slope errors, including relative errors, are retained; they are not a proof or fit of the DS limit. No floating energy pinching needed for the conserved-polynomial difference.',
        input_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names})
    with dest.open('x',encoding='utf-8') as stream:stream.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='passed',CAR_cases=len(finite),slope_cases=len(slopes),
        max_CAR_error=max(max(row[key] for key in ['first_conserved_difference_error',
                    'second_conserved_difference_error','exact_rank_covariance_error']) for row in finite))))
