"""Independent bilocal spectral-filter checks; exclusive scientific output.

Exact rational series compare the Volterra ODE with centered charged
coefficients and neutral exponential divided differences. Numerical checks
compare time quadrature with a triangular-matrix exponential at ordinary,
repeated and near-repeated real energies. No prior scientific main is run.
"""
from datetime import datetime, timezone
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path
import hashlib
import json
import platform
import time

import numpy as np
import scipy
from scipy.linalg import eigh_tridiagonal, expm
from numpy.polynomial.legendre import leggauss

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-bilocal-reconstruction/spectral-filter-checks.json'
PARAMETERS=dict(rational_q=['1/3','3/4'],rational_C=['1/2','2'],
    exact_degree=14,field_degrees=list(range(1,9)),
    rational_energies=[['-2','3'],['1/3','1/3'],['0','0']],
    physical_lambda=[0.25,0.5],Js=1.0,times=[-1.25,0.0,0.4,1.0],
    numerical_field_degrees=list(range(1,7)),quadrature_sizes=[[40,24],[64,36]],
    relative_energy_pairs=[[-0.4,0.3],[0.2,0.2],[0.2,0.200000000001]],
    floating_tolerance=3e-10,scalar_repeated_nodes=[[-1.0,-1.0,-1.0],
        [0.0,0.0,0.0],[0.7,0.7,0.7],[-0.5,-0.5,0.9],[-0.5,0.9,0.9]],
    scope='Exact rational q,C are formal. Floating spectra use q=exp(-4lambda), C=Js^2/lambda. No finite-N limit is inferred from these checks.')

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def digest(value):return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def enc(z):return [float(np.real(z)),float(np.imag(z))]

def ordinary_coefficients(q,C,nmax):
    rows=[{0:F(1)}]
    for n in range(nmax):
        row={}
        for d,v in rows[-1].items():
            row[d+1]=row.get(d+1,F(0))+v
            if d:row[d-1]=row.get(d-1,F(0))+v*C*sum(q**k for k in range(d))
        rows.append(row)
    return rows

def center(values,e):
    return [sum((F(comb(n,k))*(-e)**(n-k)*values[k] for k in range(n+1)),F(0))
            for n in range(len(values))]

def exact_checks():
    rows=[]; charged=neutral=0
    nmax=PARAMETERS['exact_degree']
    for qs in PARAMETERS['rational_q']:
        for Cs in PARAMETERS['rational_C']:
            q,C=F(qs),F(Cs); a=ordinary_coefficients(q,C,nmax)
            for xs,ys in PARAMETERS['rational_energies']:
                x,y=F(xs),F(ys);e=(x+y)/2;omega2=(y-x)**2/4
                for d in PARAMETERS['field_degrees']:
                    c=[r.get(d,F(0)) for r in a];g=center(c,e)
                    force=center([c[n+2]-(x+y)*c[n+1]+x*y*c[n] for n in range(nmax-1)],e)
                    ode_c=[F(0),F(0)];ode_n=[F(0),F(0)]
                    for n in range(nmax-1):
                        ode_c.append(force[n]+omega2*ode_c[n])
                        ode_n.append(2*C*g[n]+omega2*ode_n[n])
                    # Independent polynomial divided difference: (z^n)[x,y,z]=h_(n-2)(x,y,z).
                    uncentered=[F(0),F(0)]
                    for n in range(2,nmax+1):
                        uncentered.append(2*C*sum((x**u*y**v*a[n-2-u-v].get(d,F(0))
                            for u in range(n-1) for v in range(n-1-u)),F(0)))
                    divided=center(uncentered,e)
                    assert ode_n==divided
                    neutral+=len(divided)
                    if d>=2:
                        assert ode_c==g
                        charged+=len(g)
                    assert ode_n[2]==0 and ode_n[3]==2*C*int(d==1)
                    rows.append(dict(q=qs,C=Cs,energies=[xs,ys],degree=d,
                        charged_checked=d>=2,charged=[str(z) for z in ode_c] if d>=2 else None,
                        neutral=[str(z) for z in ode_n],exact=True))
    return dict(rows=rows,charged_coefficient_identities=charged,
                neutral_coefficient_identities=neutral)

def spectral(lam,n):
    C=1/lam;q=np.exp(-4*lam);r=2*np.sqrt(C/(1-q))
    x,U=eigh_tridiagonal(np.zeros(n),np.sqrt(C*(1-q**np.arange(1,n))/(1-q)))
    weights=U[0]**2;poly=[np.ones(n),x];norm=[1.,C]
    for d in range(1,max(PARAMETERS['numerical_field_degrees'])):
        poly.append(x*poly[-1]-C*(1-q**d)/(1-q)*poly[-2])
        norm.append(norm[-1]*C*(1-q**(d+1))/(1-q))
    return C,q,r,x,weights,np.array(poly),np.array(norm)

def divdiff_expm(t,nodes):
    T=np.diag(np.asarray(nodes,dtype=complex))+np.diag(np.ones(2),1)
    return expm(1j*t*T)[0,2]

def one_case(lam,t,pair,spectral_size,time_size):
    C,q,r,x,w,poly,norm=spectral(lam,spectral_size)
    E,Ep=r*np.asarray(pair);e=(E+Ep)/2;omega=(Ep-E)/2
    z,g=leggauss(time_size);s=t*(z+1)/2;tw=t*g/2
    phase=np.exp(1j*s[:,None]*x[None,:])
    c=(phase*w)@poly.T/norm
    cp=(phase*w*(1j*x))@poly.T/norm
    cpp=(phase*w*(-x*x))@poly.T/norm
    # np.sinc is normalized; valid for exact zero and near-zero differences.
    green=(t-s)*np.sinc(omega*(t-s)/np.pi)
    centerphase=np.exp(-1j*e*s)
    charged=-np.sum((tw*green*centerphase)[:,None]*(-cpp+2j*e*cp+E*Ep*c),axis=0)
    neutral=-2*C*np.sum((tw*green*centerphase)[:,None]*c,axis=0)
    target_charged=np.exp(-1j*e*t)*((np.exp(1j*t*x)*w)@poly.T/norm)
    dd=np.array([divdiff_expm(t,[E,Ep,xx]) for xx in x])
    target_neutral=2*C*np.exp(-1j*e*t)*((w*dd)@poly.T/norm)
    result={}
    for d in PARAMETERS['numerical_field_degrees']:
        result[d]=dict(charged=charged[d],charged_target=target_charged[d],
                       neutral=neutral[d],neutral_target=target_neutral[d])
    return result

def numerical_checks():
    rows=[];largest_identity=largest_refinement=0.
    for lam in PARAMETERS['physical_lambda']:
        for t in PARAMETERS['times']:
            for pair in PARAMETERS['relative_energy_pairs']:
                low=one_case(lam,t,pair,*PARAMETERS['quadrature_sizes'][0])
                high=one_case(lam,t,pair,*PARAMETERS['quadrature_sizes'][1])
                for d in PARAMETERS['numerical_field_degrees']:
                    identity=[abs(high[d]['neutral']-high[d]['neutral_target'])]
                    if d>=2:identity.append(abs(high[d]['charged']-high[d]['charged_target']))
                    refinement=[abs(high[d][k]-low[d][k]) for k in high[d]]
                    a=float(max(identity));b=float(max(refinement))
                    largest_identity=max(largest_identity,a);largest_refinement=max(largest_refinement,b)
                    assert max(a,b)<PARAMETERS['floating_tolerance']
                    rows.append(dict(lambda_value=lam,time=t,relative_energies=pair,degree=d,
                        charged_checked=d>=2,coarse={k:enc(v) for k,v in low[d].items()},
                        fine={k:enc(v) for k,v in high[d].items()},
                        identity_error=a,refinement_error=b))
    # Independent Hermite--Genocchi simplex integral checks repeated scalar nodes.
    z,w=leggauss(30);u=(z+1)/2;w=w/2;special=[]
    for nodes in PARAMETERS['scalar_repeated_nodes']:
        for t in PARAMETERS['times']:
            value=0j
            for a,wa in zip(u,w):
                b=(1-a)*u
                point=a*nodes[0]+b*nodes[1]+(1-a-b)*nodes[2]
                value+=wa*(1-a)*np.dot(w,-t*t*np.exp(1j*t*point))
            exact=divdiff_expm(t,nodes);error=float(abs(value-exact))
            assert error<PARAMETERS['floating_tolerance']
            if nodes[0]==nodes[1]==nodes[2]:
                assert abs(exact+t*t*np.exp(1j*t*nodes[0])/2)<1e-13
            assert abs(exact)<=t*t/2+1e-13
            special.append(dict(nodes=nodes,time=t,expm=enc(exact),simplex=enc(value),error=error))
    return dict(rows=rows,cases=len(rows),repeated_nodes=special,
        maximum_identity_error=largest_identity,maximum_refinement_error=largest_refinement,
        maximum_repeated_node_error=max(row['error'] for row in special),
        scope='Floating double-resolution quadratures and independent 3x3 matrix exponential. '
              'Near-coincident energies are evaluated without singular divisions. These are not interval certificates.')

def main():
    if OUT.exists():raise FileExistsError(OUT)
    start=time.perf_counter();exact=exact_checks();numerical=numerical_checks()
    parents=['results/susskind-memory-joint-regression/local-evidence.json',
             'code/check_memory_joint_series.py']
    result=dict(schema_version=1,status='passed',created_utc=datetime.now(timezone.utc).isoformat(),
        parameters=PARAMETERS,parameter_sha256=digest(PARAMETERS),exact=exact,numerical=numerical,
        counts=dict(charged_coefficient_identities=exact['charged_coefficient_identities'],
            neutral_coefficient_identities=exact['neutral_coefficient_identities'],
            floating_field_cases=numerical['cases'],repeated_scalar_cases=len(numerical['repeated_nodes'])),
        input_sha256={p:sha(ROOT/p) for p in ['code/check_bilocal_spectral_filter.py']+parents},
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                         platform=platform.platform()),elapsed_seconds=time.perf_counter()-start,
        scope='Spectral Volterra amplitudes, their charged cancellation and neutral divided-difference form. '
              'Does not establish the microscopic strong-return limit, higher-moment bounds or bulk mass identification.')
    OUT.parent.mkdir(parents=True,exist_ok=True)
    with OUT.open('x',encoding='utf-8') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(status='passed',counts=result['counts'],
        maximum_identity_error=numerical['maximum_identity_error'],
        maximum_refinement_error=numerical['maximum_refinement_error'],
        elapsed_seconds=result['elapsed_seconds'],sha256=sha(OUT))))

if __name__=='__main__':main()
