"""Finite-resolution dephasing of the covariant channel; fresh outputs only."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,math,platform
import numpy as np
import scipy
from scipy.linalg import eigh_tridiagonal
from scipy.integrate import simpson
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent.parent
DIR=ROOT/'results/susskind-slater-neutral-slope'
OUT=DIR/'dephasing-checks.json';FIG=DIR/'neutral-slope-dephasing.png'
if OUT.exists() or FIG.exists():raise FileExistsError('Preserve sealed science')

def spectrum(q,L):
    E,V=eigh_tridiagonal(np.zeros(L),np.sqrt((1-q**np.arange(1,L))/(1-q)))
    K=(V.T@(q**np.arange(L)[:,None]*V))*np.outer(V[0],V[0])
    return E,V,K

def diagonal_integral(q,points,terms):
    theta=np.linspace(0,math.pi,points);E=2*np.cos(theta)/math.sqrt(1-q)
    r=q**np.arange(1,terms+1)
    prod=np.exp(np.log1p(-r).sum()+np.log(1-2*np.cos(2*theta[:,None])*r+r*r).sum(axis=1))
    rho=math.sqrt(1-q)*np.sin(theta)*prod/math.pi
    jac=2*np.sin(theta)/math.sqrt(1-q)
    prev=np.zeros_like(E);current=np.ones_like(E);diag=np.ones_like(E)
    for n in range(terms):
        nxt=(E*current-math.sqrt((1-q**n)/(1-q))*prev)/math.sqrt((1-q**(n+1))/(1-q))
        diag+=q**(n+1)*nxt*nxt;prev,current=current,nxt
    return float(simpson(rho*jac,x=theta)),float(math.pi*simpson(rho*rho*diag*jac,x=theta))

rows=[];times=[1.,2.,4.,8.,12.,24.]
for q in (.0625,.25,.64):
    lam=-math.log(q)/4;s=math.sqrt(q);r=1-s*(1+lam)
    bound=r+4*s*(1+lam)/(1-q)
    mass,kappa=diagonal_integral(q,4097,192)
    mass2,kappa2=diagonal_integral(q,8193,288)
    assert abs(mass-1)<2e-12 and abs(mass2-1)<2e-12 and abs(kappa-kappa2)<2e-11
    evaluations=[]
    for L in (512,768):
        E,V,K=spectrum(q,L);delta=E[:,None]-E[None,:]
        vals=[]
        for T in times:
            terms=K*np.sinc(T*delta/math.pi)**2
            value=float(np.sum(terms));assert 0<value<=1+1e-12
            x=4*T/math.sqrt(1-q)
            logbound=(math.log(4)+x+L*math.log(x)-math.lgamma(L+1))/math.log(10)
            vals.append(dict(T=T,trace_variance=value,T_variance=T*value,
                asymptotic_ratio=T*value/kappa,one_selected_upper_bound=bound*value,
                finite_jacobi_error_bound_log10=logbound,
                signed_sum_condition=float(np.sum(abs(terms))/value)))
        evaluations.append(dict(L=L,time_rows=vals,finite_L_plateau=float(np.trace(K))))
    change=max(abs(a['trace_variance']-b['trace_variance']) for a,b in
               zip(evaluations[0]['time_rows'],evaluations[1]['time_rows']))
    assert change<2e-11
    rows.append(dict(q=q,lambda_value=lam,selected_primitive_norm=r,
        selected_density_norm=bound,kappa_u=kappa,density_mass=mass,
        diagonal_integral_refinement=abs(kappa-kappa2),cutoffs=evaluations,max_cutoff_difference=change))

# Independent triangular time integral of the one-tower correlation.
q=.25;T=4.;ts=np.linspace(0,2*T,2049);E,V,K=spectrum(q,192)
amp=V@(np.exp(1j*E[:,None]*ts)*V[0,:,None])
C=np.sum(q**np.arange(192)[:,None]*abs(amp)**2,axis=0)
time_value=float(simpson((1-ts/(2*T))*C,x=ts)/T)
ref=next(x['trace_variance'] for x in rows[1]['cutoffs'][1]['time_rows'] if x['T']==T)
assert abs(time_value-ref)<1e-10

DIR.mkdir(parents=True,exist_ok=True)
fig,axes=plt.subplots(1,2,figsize=(10,3.5),layout='constrained')
for row,color in zip(rows,('#3465a4','#c06c21','#2a8b66')):
    vals=row['cutoffs'][1]['time_rows']
    axes[0].loglog(times,[v['trace_variance'] for v in vals],'o-',color=color,label=f"q={row['q']}")
    axes[1].semilogx(times,[v['asymptotic_ratio'] for v in vals],'o-',color=color)
axes[0].set(xlabel='Averaging half-window T',ylabel='Trace dephasing variance')
axes[0].legend(frameon=False)
axes[1].axhline(1,color='.5',ls=':')
axes[1].set(xlabel='Averaging half-window T',ylabel='T variance / asymptotic coefficient')
for ax in axes:ax.grid(alpha=.2)
fig.savefig(FIG,dpi=180);plt.close(fig)
inputs=[Path(__file__),DIR/'rank-checks.json',DIR/'amplification-checks.json']
result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),rows=rows,
    independent_time_integral=dict(q=q,T=T,points=2049,Jacobi_dimension=192,
                                   absolute_difference=abs(time_value-ref)),
    parameters=dict(q=[.0625,.25,.64],times=times,Jacobi_dimensions=[512,768],
                    diagonal_quadratures=[[4097,192],[8193,288]],average='symmetric [-T,T]'),
    environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__),
    input_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
    figure_sha256={FIG.relative_to(ROOT).as_posix():hashlib.sha256(FIG.read_bytes()).hexdigest()},
    scope='Continuum finite-window dephasing diagnostics and one-selected state upper bound; not exact finite-N pinching. Analytic Jacobi error is exact-arithmetic, roundoff and density infinite-product tails are not interval enclosed.')
with OUT.open('x',encoding='utf-8') as stream:json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
print(json.dumps(dict(status='passed',physical_cases=3,time_rows=36,time_check_error=abs(time_value-ref),sha256=hashlib.sha256(OUT.read_bytes()).hexdigest())))
