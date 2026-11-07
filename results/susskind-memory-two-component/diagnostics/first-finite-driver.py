"""Finite CAR four-copy force integrals and nonzero finite-size pinching.

Uses all Gaussian coordinates, all three Wick pairings, and noncommuting
fixed Hamiltonians. Quadrature order comparisons are not interval certificates.
"""
from pathlib import Path
import datetime
import hashlib
import json
import platform
import numpy as np
from numpy.polynomial.legendre import leggauss
from memory_charged_transfer import system

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-two-component'


def phi(B,X):
    return np.einsum('aij,...jk,akl->...il',B,X,B,optimize=True)


def wick4(B,A,D,E):
    first=phi(B,A)@D@phi(B,E)
    second=phi(B,A@phi(B,D)@E)
    inside=D@B@E
    third=np.sum((B@A)@phi(B,inside),axis=0)
    return first+second+third


def hs(X):
    return float(np.linalg.norm(X)/np.sqrt(len(X)))


def integral(B,e,s,order):
    nodes,weights=leggauss(order)
    x=(nodes+1)/2;wx=weights/2
    out=np.zeros((len(e),len(e)),complex)
    for a,wa in zip(x,wx):
        yy=a*(nodes+1)/2;wy=a*weights/2
        for b,wb in zip(yy,wy):
            mats=[np.diag(np.exp(1j*s*length*e)) for length in [1-a,a-b,b]]
            out+=s*s*wa*wb*wick4(B,*mats)
    return out


if __name__=='__main__':
    OUT.mkdir(parents=True,exist_ok=True)
    dest=OUT/'four-copy-finite-time-checks.json'
    if dest.exists():raise FileExistsError(dest)
    records=[]
    for N,p,m in [(3,1,2),(4,1,2),(4,2,3)]:
        sys=system(N,p,m);original=np.array(sys['up']);d=sys['dp']
        assert np.linalg.norm(sum(b@b for b in original)-np.eye(d))<1e-13
        for seed in [17,29]:
            g=np.random.default_rng(seed).normal(size=len(original))
            h=np.einsum('a,aij->ij',g,original)
            e,U=np.linalg.eigh(h)
            B=U.conj().T@original@U
            I=np.eye(d)
            fourth=wick4(B,I,I,I)
            mu4=float(np.trace(fourth).real/d)
            assert hs(fourth-mu4*I)<2e-12
            values={}
            for s in [-.5,.5,1.0]:
                low=integral(B,e,s,8);high=integral(B,e,s,12)
                diff=hs(high-low)
                assert diff<2e-10
                # The continuous Gaussian fixtures have simple spectra;
                # exact pinching is diagonal in these fixtures only.
                mingap=float(min(abs(e[j]-e[i]) for i in range(d) for j in range(i+1,d)))
                assert mingap>1e-8
                off=high-np.diag(np.diag(high))
                assert np.linalg.norm(high,2)<=1.5*s*s+2e-12
                values[s]=high
                records.append(dict(N=N,p=p,m=m,seed=seed,C=1.,s=s,
                    Gaussian_coordinates=len(B),time_orders=[8,12],
                    time_order_difference=diff,integral_HS_norm=hs(high),
                    pinched_residual_HS_norm=hs(off),
                    integral_operator_norm=float(np.linalg.norm(high,2)),
                    analytic_operator_bound=1.5*s*s,
                    fourth_moment=mu4,min_spectral_gap=mingap))
            assert hs(values[-.5]-values[.5].conj().T)<2e-12
    # Finite-N regression is not an exact conserved identity.
    assert max(row['pinched_residual_HS_norm'] for row in records)>1e-5
    names=['code/check_four_copy_finite_time.py','code/memory_charged_transfer.py',
           'code/complex_syk_dictionary.py']
    result=dict(records=records,
        created_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        environment=dict(python=platform.python_version(),numpy=np.__version__),
        scope='18 finite CAR word integrals using every Gaussian coordinate and all Wick pairings; two time quadrature orders, adjoint symmetry, scalar zero-time fourth moment, and finite-N nonconservation. No DS fit or certified quadrature error; asymptotic vanishing is an analytic result.',
        input_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names})
    with dest.open('x',encoding='utf-8') as stream:stream.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='passed',finite_time_rows=len(records),
        max_order_difference=max(x['time_order_difference'] for x in records),
        max_finite_N_dynamic_norm=max(x['pinched_residual_HS_norm'] for x in records))))
