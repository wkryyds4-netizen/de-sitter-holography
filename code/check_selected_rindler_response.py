"""Fresh exact Taylor and finite-h packet checks; refuses sealed output."""
from pathlib import Path
from datetime import datetime,timezone
from fractions import Fraction
from math import factorial,ceil,log
import hashlib,json,platform
import numpy as np
import scipy
import sympy as sp
from scipy.linalg import eigh_tridiagonal
from selected_rindler_response import gap_difference,exponential_coefficient,packets,response,packet_tail_log

ROOT=Path(__file__).resolve().parents[1]
DIR=ROOT/'results/susskind-selected-rindler-response'
OUT=DIR/'response-checks.json'
if OUT.exists():raise FileExistsError(OUT)
DIR.mkdir(exist_ok=True)
u=sp.symbols('u')
pred=sp.series((u*sp.tanh(u)-5*sp.tanh(u)**2)/sp.cosh(u)**4,u,0,16).removeO().expand()
exact=[];zeros=0
for n in range(0,16,2):
    poly=gap_difference(n);k=n//2
    for j in range(k+1):
        assert exponential_coefficient(poly,j)==0,(n,j)
        zeros+=1
    coefficient=(-1)**k*exponential_coefficient(poly,k+1)/factorial(n)
    expected=Fraction(str(pred.coeff(u,n)))
    assert coefficient==expected,(n,coefficient,expected)
    exact.append(dict(time_degree=n,q_polynomial={str(k):v for k,v in sorted(poly.items())},response_coefficient=str(coefficient)))

times=np.linspace(0,1.5,7);rows=[]
for h in (.08,.04,.02,.01):
    L=ceil(8*times[-1]/h)+64
    c0,diff,cond,norm=packets(h,L,times)
    c02,diff2,cond2,norm2=packets(h,L+64,times)
    assert max(np.max(abs(c0-c02)),np.max(abs(diff-diff2)))<2e-11
    assert np.max(abs(diff2.imag))<2e-11
    assert norm2<2e-11
    bounds=[None if t==0 else packet_tail_log(h,L,float(t))/log(10) for t in times]
    assert max(x for x in bounds if x is not None)<-25
    entries=[]
    for i,t in enumerate(times):
        target=response(float(t));observed=float(diff2[i].real/h)
        entries.append(dict(u=float(t),trace=float(c02[i]),selected=float(c02[i]+diff2[i].real),
            difference=float(diff2[i].real),response_over_h=observed,prediction=target,
            coefficient_error=observed-target,relative_coefficient_error=None if not target else abs((observed-target)/target),
            cutoff_response_difference=float(abs(diff2[i]-diff[i])/h),
            phase_condition=None if t==0 else float(cond2[i]),tail_error_log10=bounds[i]))
    independent=None
    if h==.08:
        q=np.exp(-h);size=L+64
        energies,V=eigh_tridiagonal(np.zeros(size),np.sqrt(-np.expm1(-h*np.arange(1,size))/(-np.expm1(-h)*h)))
        phase=np.exp(1j*np.outer(energies,times));a=V@(V[0,:,None]*phase);b=V@(V[2,:,None]*phase)
        alternative=np.sqrt(q*(1+q))*np.sum(np.exp(-2*h*np.arange(size))[:,None]*np.conj(b)*a,axis=0)
        independent=float(np.max(abs(alternative-diff2)))
        assert independent<2e-11
    rows.append(dict(h=h,lambda_value=h/4,dimensions=[L,L+64],norm_residual=norm2,
                     independent_eigensolver_difference=independent,values=entries))
    print('h',h,'maximum response discrepancy',max(abs(x['coefficient_error']) for x in entries),flush=True)
result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),
    exact_taylor=exact,laurent_cancellations=zeros,rows=rows,
    parameters=dict(h='4 lambda',energy='K=H/(2 J_s), variance 1/h',time='u=2 J_s t',
                    limits='first N,p double scaling at fixed lambda, then h to zero at fixed u',
                    state='one selected excitation restricted to h,x; density 1-sqrt(q)+sqrt(q) h_normalized^2'),
    environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__),
    input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in [
        'code/selected_rindler_response.py','code/check_selected_rindler_response.py',
        'results/susskind-slater-selected-state/exact-checks.json',
        'manuscript/susskind-slater-selected-state/main.tex']},
    scope='Exact polynomial Taylor checks and finite-deformation empirical diagnostics of a semiclassical response prediction; no proof of uniform differentiated saddle remainder, no mass identification.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
print('passed',hashlib.sha256(OUT.read_bytes()).hexdigest(),flush=True)
