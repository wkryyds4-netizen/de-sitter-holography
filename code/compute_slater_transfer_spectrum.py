"""q-Jacobi quadrature and selected-state spectral/preparation kernels."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,platform,math
import numpy as np
import scipy
from scipy.linalg import eigh_tridiagonal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent.parent
DIR=ROOT/'results/susskind-slater-transfer'
OUT=DIR/'spectral-checks.json';FIG=DIR/'selected-spectrum.png'
if OUT.exists() or FIG.exists():raise FileExistsError('Preserve existing scientific outputs')
def jacobi(q,L):
    off=np.sqrt((1-q**np.arange(1,L))/(1-q))
    E,V=eigh_tridiagonal(np.zeros(L),off)
    return E,V[0,:]**2
def density(q,theta):
    r=q**np.arange(1,257)
    logprod=np.log1p(-r).sum()
    factors=1-2*np.cos(2*theta[:,None])*r[None,:]+r[None,:]**2
    dtheta=np.exp(logprod+np.log(factors).sum(axis=1))*4*np.sin(theta)**2/(2*np.pi)
    E=2*np.cos(theta)/np.sqrt(1-q)
    return E,dtheta*np.sqrt(1-q)/(2*np.sin(theta)),dtheta
rows=[]
for q in [0.0625,0.25,0.64]:
    s=math.sqrt(q);R=2/math.sqrt(1-q)
    E,w=jacobi(q,160);E2,w2=jacobi(q,240)
    selected=w*(1-s+s*E*E);selected2=w2*(1-s+s*E2*E2)
    theta=np.linspace(1e-6,math.pi-1e-6,8193)
    Ed,rho,rhot=density(q,theta)
    mass=float(np.trapezoid(rhot,theta));second=float(np.trapezoid(Ed**2*rhot,theta))
    assert abs(mass-1)<1e-10 and abs(second-1)<1e-10
    assert abs(selected.sum()-1)<1e-12
    times=[]
    for t in [0,0.5,1,2,4]:
        vacuum=np.sum(w*np.exp(1j*t*E));value=np.sum(selected*np.exp(1j*t*E))
        other=np.sum(selected2*np.exp(1j*t*E2))
        difference=abs(value-other)
        assert difference<3e-12
        K=318
        logbound=math.log(2)+t*R+K*math.log(t*R)-math.lgamma(K+1) if t else -math.inf
        times.append(dict(t=t,vacuum_real=float(vacuum.real),selected_real=float(value.real),
            imaginary_roundoff=float(abs(value.imag)),L160_L240_difference=float(difference),
            exact_arithmetic_quadrature_tail_log10=logbound/math.log(10)))
    prep=[]
    for beta in [0,0.5,1,2]:
        filtered=selected*np.exp(-beta*E);Z=filtered.sum()
        prep.append(dict(beta=beta,Z=float(Z),mean_energy=float(np.sum(E*filtered)/Z),
                         energy_variance=float(np.sum(E*E*filtered)/Z-(np.sum(E*filtered)/Z)**2)))
    tails=[]
    for L in [8,16,32]:
        r=q*q
        upper=s*r**L*((L+1)**2/(1-r)+2*(L+1)*r/(1-r)**2+r*(1+r)/(1-r)**3)
        actual=sum(s*(n+1)**2*q**(2*n)/((1-q**(n+1))/(1-q)) for n in range(L,L+300))
        assert actual<=upper*(1+1e-14)
        tails.append(dict(L=L,HS_squared_upper=upper,sum_300_diagnostic=actual))
    rows.append(dict(q=q,gamma=q**0.25,support_radius=R,density_mass=mass,density_second=second,
        selected_mass=float(selected.sum()),time_rows=times,preparation_rows=prep,kernel_tails=tails))

q=.25;s=math.sqrt(q);theta=np.linspace(1e-6,math.pi-1e-6,1601)
E,rho,_=density(q,theta);nodes,weights=jacobi(q,240)
ts=np.linspace(0,6,301);waves=np.cos(ts[:,None]*nodes[None,:])
fig,axes=plt.subplots(1,2,figsize=(10.5,3.6),layout='constrained')
axes[0].plot(E[::-1],rho[::-1],label='Reference state',color='#3465a4',lw=2)
axes[0].plot(E[::-1],((1-s+s*E*E)*rho)[::-1],label='Selected excitation',color='#2a8b66',lw=2)
axes[0].set(xlabel='Dimensionless energy E',ylabel='Spectral density',title='Spectral redistribution, q = 0.25')
axes[0].legend(frameon=False,fontsize=9)
axes[1].plot(ts,waves@weights,label='R(t)',color='#3465a4',lw=2)
axes[1].plot(ts,waves@(weights*(1-s+s*nodes**2)),label='C_M(t)',color='#2a8b66',lw=2)
axes[1].axhline(0,color='.65',lw=.6)
axes[1].set(xlabel='Dimensionless time t',ylabel='Real return amplitude',title='Finite-time returns')
axes[1].legend(frameon=False,fontsize=9)
for ax in axes:ax.grid(alpha=.18)
fig.savefig(FIG,dpi=180);plt.close(fig)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
inputs=[Path(__file__),DIR/'exact-checks.json']
result=dict(status='passed',created_utc=datetime.now(timezone.utc).isoformat(),rows=rows,
 parameters=dict(q=[.0625,.25,.64],Jacobi_sizes=[160,240],density_theta_points=8193,
 q_product_terms=256,time=[0,.5,1,2,4],beta=[0,.5,1,2],figure_q=.25),
 environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__,dtype='float64'),
 input_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in inputs},figure_sha256={FIG.relative_to(ROOT).as_posix():sha(FIG)},
 scope='Floating-point spectral calculations with independent density normalization and quadrature-size comparisons. Analytic Gaussian-quadrature Taylor tail bounds apply in exact arithmetic; no interval enclosure of eigensolver roundoff is claimed. Euclidean filtering is vector preparation, not a KMS-state identification.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2,allow_nan=True);f.write('\n')
print(json.dumps(dict(status='passed',q_cases=len(rows),times=15,preparations=12,sha256=sha(OUT),figure_sha256=sha(FIG))))
