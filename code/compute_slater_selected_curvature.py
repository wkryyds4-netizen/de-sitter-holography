"""Fresh physical-curve curvature dynamics, two independent representations."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,math,platform
import numpy as np
import scipy
from scipy.linalg import eigh_tridiagonal
from scipy.special import roots_legendre
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent.parent
DIR=ROOT/'results/susskind-slater-selected-state'
OUT=DIR/'curvature-dynamics.json';FIG=DIR/'selected-curvature.png'
if OUT.exists() or FIG.exists():raise FileExistsError('Preserve successful output')
assert (DIR/'exact-checks.json').exists()

def jacobi(q,L,times):
    E,V=eigh_tridiagonal(np.zeros(L),np.sqrt((1-q**np.arange(1,L))/(1-q)))
    phases=np.exp(1j*E[:,None]*np.asarray(times)[None,:])
    a=V@(phases*V[0,:,None]);b=V@(phases*V[2,:,None]);r=q*q
    c0=np.sum(r**np.arange(L)[:,None]*abs(a)**2,axis=0)
    cross=np.sum(r**np.arange(L)[:,None]*np.conj(b)*a,axis=0)
    c1=c0+math.sqrt(q*(1+q))*cross
    assert np.max(abs(c1.imag))<2e-12
    return c0,c1.real,float(np.max(abs(np.sum(abs(a)**2,axis=0)-1)))

def quadrature(q,points,terms,times):
    x,weights=roots_legendre(points);theta=(x+1)*math.pi/2;weights*=math.pi/2
    logmu=np.log(4*np.sin(theta)**2/(2*math.pi))
    for n in range(1,terms+1):
        z=q**n;logmu+=math.log1p(-z)+np.log(1-2*z*np.cos(2*theta)+z*z)
    mu=weights*np.exp(logmu);E=2*np.cos(theta)/math.sqrt(1-q);r=q*q
    cp=np.cos(theta[:,None]+theta[None,:]);cm=np.cos(theta[:,None]-theta[None,:])
    logK=np.zeros((points,points))
    for n in range(terms):
        z=r*q**n
        logK+=math.log1p(-r*r*q**n)-np.log(1-2*z*cp+z*z)-np.log(1-2*z*cm+z*z)
    measure=np.exp(logK)*np.outer(mu,mu)
    selected=(1-math.sqrt(q)+math.sqrt(q)*E**2)[:,None]*measure
    delta=E[:,None]-E[None,:]
    result=[]
    for t in times:
        c0=np.sum(measure*np.exp(1j*t*delta));parts=selected*np.exp(1j*t*delta);c1=np.sum(parts)
        result.append(dict(t=float(t),trace=float(c0.real),selected=float(c1.real),
            imaginary_residual=float(abs(c1.imag)),
            phase_sum_condition=float(np.sum(abs(parts))/abs(c1))))
    log_tail=6*q**(terms+1)/((1-q)*(1-q**(terms+1)))+(r*r+4*r)*q**terms/((1-q)*(1-r*q**terms))
    return result,dict(mu_mass=float(mu.sum()),two_energy_mass=float(measure.sum()),
                       selected_mass=float(selected.sum()),uniform_product_log_tail=log_tail)

def log_tail(q,L,t):
    if not t:return None
    x=abs(t)*2/math.sqrt(1-q)
    # A common upper bound for the |0> and |2> packet errors.
    return (math.log(4*(1+math.sqrt(q*(1+q))))+x+(L-2)*math.log(x)-math.lgamma(L-1))/math.log(10)

times=np.array([0.,.5,1.,2.,4.,6.]);rows=[]
for q in (.0625,.25,.64):
    lam=-math.log(q)/4;s=math.sqrt(q);r=q*q
    first=jacobi(q,192,times);second=jacobi(q,256,times)
    quad1,mass1=quadrature(q,192,128,times)
    quad2,mass2=quadrature(q,288,192,times)
    assert max(abs(mass2[k]-1) for k in ['mu_mass','two_energy_mass','selected_mass'])<3e-12
    records=[]
    for i,t in enumerate(times):
        c0,c1=float(second[0][i]),float(second[1][i])
        cutoff=float(abs(c1-first[1][i]));quaderr=abs(c1-quad2[i]['selected'])
        refine=abs(quad1[i]['selected']-quad2[i]['selected'])
        assert max(cutoff,quaderr,refine)<3e-11
        records.append(dict(t=float(t),trace=c0,selected=c1,selected_minus_trace=c1-c0,
            cutoff_difference=cutoff,quadrature_difference=quaderr,quadrature_refinement=refine,
            relative_cutoff_difference=cutoff/max(abs(c1),1e-15),
            near_zero=bool(abs(c1)<1e-12),phase_sum_condition=quad2[i]['phase_sum_condition'],
            analytic_L256_error_log10=log_tail(q,256,t)))
    small=np.array([1e-3,5e-4]);_,csmall,_=jacobi(q,256,small)
    est=2*(1-csmall)/small**2;expected=2*(1-r)+s*(1+q)*(1-r)**2
    assert max(abs(est-expected))<6e-6
    rows.append(dict(q=q,lambda_value=lam,selected_density_lower_bound=1-s*(1+lam),
        reduced_density_lower_bound=1-s,
        density_nonstationarity_norm=2*q*lam*(1+lam)*(1+q)*(1-q*q),
        commutator_norm=expected,finite_difference_estimates=est.tolist(),
        wavepacket_norm_residual=second[2],quadrature_masses=[mass1,mass2],time_rows=records))
    print('physical q',q,'max kernel difference',max(x['quadrature_difference'] for x in records),flush=True)

ts=np.linspace(0,6,241);fig,axs=plt.subplots(1,2,figsize=(10,3.7),layout='constrained')
c0,c1,_=jacobi(.25,256,ts)
axs[0].plot(ts,c0,label='Trace state',color='#3465a4',lw=2)
axs[0].plot(ts,c1,label='One selected excitation',color='#c06c21',lw=2)
axs[0].set(xlabel='Dimensionless time t',ylabel='Normalized curvature correlation',title='q = 0.25')
axs[0].legend(frameon=False)
for q,color in zip((.0625,.25,.64),('#3465a4','#c06c21','#2a8b66')):
    a,b,_=jacobi(q,256,ts);axs[1].plot(ts,b-a,label=f'q = {q}',color=color,lw=2)
axs[1].set(xlabel='Dimensionless time t',ylabel='Selected minus trace correlation')
axs[1].legend(frameon=False)
for ax in axs:ax.axhline(0,color='.65',lw=.6);ax.grid(alpha=.2)
fig.savefig(FIG,dpi=180);plt.close(fig)
inputs=[Path(__file__),DIR/'exact-checks.json']
result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),rows=rows,
    parameters=dict(physical_q=[.0625,.25,.64],lambda_relation='lambda=-log(q)/4',
        Jacobi_dimensions=[192,256],quadratures=[[192,128],[288,192]],times=times.tolist(),
        small_times=[.001,.0005],figure_times=[0,6,241],probe_crossing='r=q^2'),
    environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                     matplotlib=matplotlib.__version__,dtype='float64/complex128'),
    input_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
    figure_sha256={FIG.relative_to(ROOT).as_posix():hashlib.sha256(FIG.read_bytes()).hexdigest()},
    scope='Physical-curve diagnostics of the exact limiting selected curvature kernel. Jacobi tail and product-log tail are exact-arithmetic bounds only; floating roundoff and quadrature discretization are not interval enclosed. No finite-N convergence rate or finite-G scattering is inferred.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
print(json.dumps(dict(status='passed',physical_cases=3,time_rows=18,sha256=hashlib.sha256(OUT.read_bytes()).hexdigest())),flush=True)
