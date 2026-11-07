"""Independent q-Hermite, force, centering and time-coefficient checks.

Exact rational polynomial identities and floating spectral quadrature have
separate scopes. Output is exclusive-create; no historical main is used.
"""
from datetime import datetime, timezone
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path
import hashlib
import json
import platform
import time

import mpmath as mp
import numpy as np
import scipy
from scipy.linalg import eigh_tridiagonal

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-joint-regression/joint-series-checks.json'
PARAMETERS=dict(rational_q=['1/5','1/3','3/4'], rational_C=['1/2','2'],
                maximum_degree=14, force_degree=8,
                spectral_energy_pairs=[['-2','3'],['1/3','1/3'],['0','0']],
                physical_lambda=[0.125,0.25,0.5], Js=1.0,
                times=[-1.5,-0.4,0.0,0.7,1.5], spectral_degrees=list(range(9)),
                derivative_orders=[0,1,2], quadrature_sizes=[48,72],
                floating_comparison_tolerance=2e-10,
                scope='Rational q,C checks are formal algebra. Quadrature uses q=exp(-4 lambda), C=Js^2/lambda.')


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def digest(x):
    return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def add(a,b,factor=F(1)):
    out=dict(a)
    for d,x in b.items():
        out[d]=out.get(d,F(0))+factor*x
    return {d:x for d,x in out.items() if x}


def qinteger(n,q):
    return sum((q**k for k in range(n)),F(0))


def expansions(q,C,maximum):
    # Ordinary monomials expanded in monic q-Hermites.
    a=[{0:F(1)}]
    # Independently, q-Hermites expanded in ordinary monomials.
    h=[{0:F(1)},{1:F(1)}]
    for n in range(maximum):
        out={}
        for d,x in a[-1].items():
            out=add(out,{d+1:x})
            if d:
                out=add(out,{d-1:C*qinteger(d,q)*x})
        a.append(out)
    for n in range(1,maximum):
        h.append(add({d+1:x for d,x in h[-1].items()},h[-2],-C*qinteger(n,q)))
    for n,row in enumerate(a):
        recovered={}
        for d,x in row.items():
            recovered=add(recovered,h[d],x)
        assert recovered=={n:F(1)}
        assert all(x>0 for x in row.values())
    return a,h


def scalar_convolution(a,b,maximum):
    # Exponential generating coefficients in powers (is)^n/n!.
    return [sum((F(comb(n,k))*a[k]*b[n-k] for k in range(n+1)),F(0))
            for n in range(maximum+1)]


def exact_checks():
    rows=[]
    force_count=0
    for qstr in PARAMETERS['rational_q']:
        for Cstr in PARAMETERS['rational_C']:
            q,C=F(qstr),F(Cstr)
            maximum=PARAMETERS['maximum_degree']
            a,h=expansions(q,C,maximum)
            # A finite scalar Hilbert-Schmidt eigencoordinate with left energy x,
            # right energy y independently checks Q1 and center phases.
            checked=[]
            for xs,ys in PARAMETERS['spectral_energy_pairs']:
                x,y=F(xs),F(ys)
                midpoint=(x+y)/2
                nmax=PARAMETERS['force_degree']
                for d in range(nmax+3):
                    c=[row.get(d,F(0)) for row in a]
                    force=[c[n+2]-(x+y)*c[n+1]+x*y*c[n] for n in range(nmax+1)]
                    centered=scalar_convolution([(-midpoint)**n for n in range(nmax+1)],force,nmax)
                    # Apply factorized -(partial_s-i*x)(partial_s-i*y)
                    # to c(s), then multiply each half-time phase separately.
                    right_factor=[c[n+1]-y*c[n] for n in range(nmax+2)]
                    independent=[right_factor[n+1]-x*right_factor[n] for n in range(nmax+1)]
                    independent=scalar_convolution([(-x/2)**n for n in range(nmax+1)],independent,nmax)
                    independent=scalar_convolution([(-y/2)**n for n in range(nmax+1)],independent,nmax)
                    assert independent==centered
                    force_count+=len(centered)
                    if d>=2:
                        assert centered[0]==int(d==2)
                        assert centered[1]==F(int(d==3))-3*midpoint*int(d==2)
                        expect=F(int(d==4))-4*midpoint*int(d==3)
                        expect+=(C*(3+2*q+q*q)+x*y+5*midpoint*midpoint)*int(d==2)
                        assert centered[2]==expect
                    checked.append(dict(energies=[xs,ys],d=d,
                                        centered_first_three=[str(v) for v in centered[:3]]))
            rows.append(dict(q=qstr,C=Cstr,degree_identities=len(a),
                             a_rows=[{str(k):str(v) for k,v in r.items()} for r in a],
                             force_checks=checked))
    return dict(rows=rows,ordinary_Hermite_identities=sum(r['degree_identities'] for r in rows),
                centered_force_coefficient_identities=force_count,
                low_coefficients='A0=J2; A1=J3-3 e J2; A2=J4-4 e J3+[C(3+2q+q^2)+EEprime+5e^2]J2, e=(E+Eprime)/2; J0=J1=0')


def quadrature_coefficients(lam,size):
    C=1/lam
    q=np.exp(-4*lam)
    diagonal=np.zeros(size)
    off=np.sqrt(C*(1-q**np.arange(1,size))/(1-q))
    energy,U=eigh_tridiagonal(diagonal,off)
    weights=U[0]**2
    R=2*np.sqrt(C/(1-q))
    assert np.max(np.abs(energy))<=R*(1+1e-14)
    norm=[1.0]
    for d in range(1,9):
        norm.append(norm[-1]*C*(1-q**d)/(1-q))
    poly=[np.ones_like(energy),energy]
    for d in range(1,8):
        poly.append(energy*poly[-1]-C*(1-q**d)/(1-q)*poly[-2])
    result={}
    for t in PARAMETERS['times']:
        phase=np.exp(1j*t*energy)
        for d in range(9):
            for j in PARAMETERS['derivative_orders']:
                result[t,d,j]=np.dot(weights,phase*(1j*energy)**j*poly[d])/norm[d]
    return result,dict(q=float(q),C=C,spectral_radius=float(R),
                       weight_sum=float(weights.sum()),maximum_energy=float(np.max(np.abs(energy))))


def analytic_truncation_bound(lam,size,t,d,j):
    # Only polynomial/Taylor quadrature truncation, not floating roundoff.
    with mp.workdps(80):
        lam=mp.mpf(str(lam)); C=1/lam; q=mp.exp(-4*lam)
        R=2*mp.sqrt(C/(1-q)); bounds=[mp.mpf(1),R]; norm=mp.mpf(1)
        for k in range(1,d+1):
            norm*=C*(1-q**k)/(1-q)
        for k in range(1,d):
            bounds.append(R*bounds[-1]+C*(1-q**k)/(1-q)*bounds[-2])
        order=2*size-d-j
        value=2*R**j*bounds[d]/norm*(abs(mp.mpf(str(t)))*R)**order/mp.factorial(order)
        return mp.nstr(value,20)


def numerical_checks():
    rows=[]
    maximum=0.0
    for lam in PARAMETERS['physical_lambda']:
        small,meta=quadrature_coefficients(lam,PARAMETERS['quadrature_sizes'][0])
        large,_=quadrature_coefficients(lam,PARAMETERS['quadrature_sizes'][1])
        for key,x in small.items():
            y=large[key]
            err=float(abs(x-y)); maximum=max(maximum,err)
            assert err<PARAMETERS['floating_comparison_tolerance']
            t,d,j=key
            # e^{isH} coefficients respect c_d(-s)=conjugate(c_d(s)).
            rows.append(dict(lambda_value=lam,time=t,degree=d,derivative=j,
                             value_48=[float(x.real),float(x.imag)],
                             value_72=[float(y.real),float(y.imag)],absolute_difference=err,
                             analytic_Taylor_quadrature_tail_48=analytic_truncation_bound(lam,48,t,d,j)))
    return dict(rows=rows,cases=len(rows),maximum_absolute_difference=maximum,
                scope='Two floating Jacobi quadratures. Analytic bounds cover Taylor/quadrature truncation only; floating roundoff is not interval-certified. No finite-N dynamics is simulated.')


def main():
    if OUT.exists():
        raise FileExistsError(OUT)
    started=time.perf_counter()
    exact=exact_checks()
    numerical=numerical_checks()
    record=dict(schema_version=1,status='passed',created_utc=datetime.now(timezone.utc).isoformat(),
                parameters=PARAMETERS,parameter_sha256=digest(PARAMETERS),exact=exact,numerical=numerical,
                counts=dict(ordinary_Hermite_identities=exact['ordinary_Hermite_identities'],
                            centered_force_coefficient_identities=exact['centered_force_coefficient_identities'],
                            quadrature_coefficients=numerical['cases']),
                input_sha256={'code/check_memory_joint_series.py':sha(Path(__file__))},
                environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                                 mpmath=mp.__version__,platform=platform.platform()),
                elapsed_seconds=time.perf_counter()-started,
                scope='Exact ordinary/spectral polynomial and centered-force coefficient identities, plus independent numerical evaluation of c_d and its first two derivatives. This does not prove the charged pinching theorem or substitute for common-J mixed Gram checks.')
    OUT.parent.mkdir(parents=True,exist_ok=True)
    with OUT.open('x',encoding='utf-8') as stream:
        json.dump(record,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(status='passed',counts=record['counts'],
                          maximum_quadrature_difference=numerical['maximum_absolute_difference'],
                          output_sha256=sha(OUT),elapsed_seconds=record['elapsed_seconds'])))


if __name__=='__main__':
    main()
