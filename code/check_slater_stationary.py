"""Finite-resolution diagnostics for the stationary-state analysis.

No sealed generator is imported or rerun. Infinite-dimensional statements
are proved in the manuscript, not inferred from these finite matrices.
"""
from pathlib import Path
from datetime import datetime, timezone
from fractions import Fraction
from itertools import combinations
import hashlib, json, math, platform
import numpy as np
import scipy
from scipy.linalg import eigh_tridiagonal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
DIR = ROOT/'results/susskind-slater-stationary'
DIR.mkdir(parents=True, exist_ok=True)
OUT = DIR/'stationary-checks.json'
FIG = DIR/'finite-time-states.png'
if OUT.exists() or FIG.exists(): raise FileExistsError('Preserve successful outputs')

def spectrum(q, L):
    E, V = eigh_tridiagonal(np.zeros(L), np.sqrt((1-q**np.arange(1,L))/(1-q)))
    return E, V[0]**2

def density_data(q, points, terms):
    theta = np.linspace(0, math.pi, points)
    r = q**np.arange(1, terms+1)
    prod = np.exp(np.log1p(-r).sum()+np.log(1-2*np.cos(2*theta[:,None])*r+r*r).sum(axis=1))
    rho = math.sqrt(1-q)*np.sin(theta)*prod/math.pi
    jac = 2*np.sin(theta)/math.sqrt(1-q)
    return float(np.trapezoid(rho*jac, theta)), float(math.pi*np.trapezoid(rho*rho*jac, theta))

rows = []
for q in (.0625, .25, .64):
    mass, kappa = density_data(q, 8193, 256)
    mass2, kappa2 = density_data(q, 16385, 384)
    assert abs(mass-1)<2e-12 and abs(mass2-1)<2e-12
    assert abs(kappa-kappa2)<2e-12
    evaluations = []
    for L in (800, 1100):
        E, w = spectrum(q,L)
        delta = E[:,None]-E[None,:]
        ww = w[:,None]*w[None,:]
        times = []
        for T in (2.,4.,8.,16.,32.,64.):
            kernel = np.sinc(T*delta/math.pi)
            signal = float(np.sum(ww*kernel))
            purity = float(np.sum(ww*kernel**2))
            assert 0<signal<=1 and 0<purity<=1
            times.append(dict(T=T, vacuum_projection=signal, purity=purity,
                              T_times_projection=T*signal,T_times_purity=T*purity,
                              projection_relative_to_asymptotic=T*signal/kappa,
                              purity_relative_to_asymptotic=T*purity/kappa,
                              signed_sum_condition=float(np.sum(abs(ww*kernel))/signal)))
        evaluations.append(dict(L=L, time_rows=times, diagonal_pinching_plateau=float(np.sum(w*w))))
    differences = [max(abs(a[k]-b[k]) for k in ('vacuum_projection','purity'))
                   for a,b in zip(evaluations[0]['time_rows'],evaluations[1]['time_rows'])]
    assert max(differences)<2e-11
    rows.append(dict(q=q,kappa=kappa,density_mass=mass,density_refinement=abs(kappa-kappa2),
                     quadratures=evaluations,max_L_difference=max(differences)))

# Independent time-domain integral checks of sinc and triangular kernels.
E,w=spectrum(.25,1100)
time_checks=[]
for T in (2.,8.,32.):
    u=np.linspace(0,2*T,8193)
    R=np.cos(u[:,None]*E)@w
    projection=float(np.trapezoid(R[:4097]**2,u[:4097])/T)
    purity=float(np.trapezoid((1-u/(2*T))*R**2,u)/T)
    saved=next(r for r in rows[1]['quadratures'][1]['time_rows'] if r['T']==T)
    errors=dict(projection=abs(projection-saved['vacuum_projection']),purity=abs(purity-saved['purity']))
    assert max(errors.values())<4e-7
    time_checks.append(dict(T=T,errors=errors,points=8193))

# Enumerate occupancy projectors independently of the closed trace counts.
counts=[]
for N,p in ((4,1),(6,1),(6,2),(8,2),(10,3)):
    m=N//2; O=tuple(range(m)); empty=tuple(range(m,N))
    sector=[x for x in range(1<<N) if x.bit_count()==m]
    full_total=central_total=0
    labels=0
    for A in combinations(empty,p):
        for B in combinations(O,p):
            a=sum(1<<i for i in A);b=sum(1<<i for i in B)
            full_total+=sum((x&a)==0 and (x&b)==b for x in range(1<<N))
            central_total+=sum((x&a)==0 and (x&b)==b for x in sector)
            labels+=1
    full=Fraction(full_total,labels*(1<<N))
    central=Fraction(central_total,labels*len(sector))
    assert full==Fraction(1,4**p)
    assert central==Fraction(math.comb(N-2*p,m-p),math.comb(N,m))
    counts.append(dict(N=N,p=p,labels=labels,full=str(full),central=str(central)))

fig,axes=plt.subplots(1,2,figsize=(10,3.5),layout='constrained')
for row in rows:
    data=row['quadratures'][1]['time_rows'];T=[d['T'] for d in data]
    axes[0].semilogx(T,[d['T_times_projection'] for d in data],marker='o',label=f"q={row['q']}")
    axes[1].semilogx(T,[d['T_times_purity'] for d in data],marker='o',label=f"q={row['q']}")
    for ax in axes:ax.axhline(row['kappa'],ls=':',alpha=.5)
axes[0].set(xlabel='Averaging half-width T',ylabel='T times vacuum probability')
axes[1].set(xlabel='Averaging half-width T',ylabel='T times purity')
for ax in axes:ax.grid(alpha=.2);ax.legend(frameon=False)
fig.savefig(FIG,dpi=180);plt.close(fig)

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
inputs=[Path(__file__),ROOT/'results/susskind-slater-transfer/local-evidence.json',
        ROOT/'manuscript/susskind-slater-transfer/main.tex',ROOT/'notes/susskind-slater-ladder-proof.md']
result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),rows=rows,
            independent_time_checks=time_checks,occupancy_counts=counts,
            parameters=dict(q=[.0625,.25,.64],Jacobi_sizes=[800,1100],T=[2,4,8,16,32,64],
                            density_grids=[8193,16385],q_product_terms=[256,384]),
            environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__,dtype='float64'),
            input_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in inputs},
            figure_sha256={FIG.relative_to(ROOT).as_posix():sha(FIG)},
            scope='Finite quadrature/time-integral diagnostics and exact occupancy counts. No numerical proof of stationary-state absence, no interval certification, no long-time-uniform finite-N assertion.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
print(json.dumps(dict(status='passed',q_cases=3,time_rows=36,independent_time_checks=3,occupancy_systems=5,
                     maximum_quadrature_difference=max(r['max_L_difference'] for r in rows),
                     kappa=[r['kappa'] for r in rows],sha256=sha(OUT))))
