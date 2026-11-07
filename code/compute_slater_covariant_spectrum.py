"""Physical-parameter checks and stationary covariant-channel spectrum."""
from pathlib import Path
from datetime import datetime,timezone
import math,json,hashlib,platform
import numpy as np
import scipy
from scipy.linalg import eigh_tridiagonal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent.parent
DIR=ROOT/'results/susskind-slater-covariant'
OUT=DIR/'spectral-checks.json';FIG=DIR/'covariant-channel.png'
if OUT.exists() or FIG.exists():raise FileExistsError('Preserve sealed outputs')
def calc(q,L,times):
    E,V=eigh_tridiagonal(np.zeros(L),np.sqrt((1-q**np.arange(1,L))/(1-q)))
    amplitudes=V@(np.exp(1j*E[:,None]*np.asarray(times)[None,:])*V[0,:,None])
    return np.sum(q**np.arange(L)[:,None]*abs(amplitudes)**2,axis=0),np.sum(abs(amplitudes)**2,axis=0)
rows=[];times=np.array([0,.5,1,2,4,6.])
for q in (.0625,.25,.64):
    lam=-math.log(q)/4;g=q**.25;s=math.sqrt(q);d=math.sqrt(lam)
    residual=1-lam/math.expm1(2*lam)
    gram=np.array([[1,0,g],[0,1,-d*g],[g,-d*g,1]])
    assert np.linalg.eigvalsh(gram)[0]>0 and .5<residual<1
    a,an=calc(q,80,times);b,bn=calc(q,120,times)
    assert max(abs(a-b))<2e-12 and max(abs(bn-1))<2e-12
    Emax=2/math.sqrt(1-q)
    bounds=[]
    for t in times:
        bounds.append(None if not t else (math.log(4)+t*Emax+120*math.log(t*Emax)-math.lgamma(121))/math.log(10))
    small,_=calc(q,120,np.array([0.,1e-3,5e-4]))
    curvature=2*(1-small[1:])/np.array([1e-3,5e-4])**2
    assert max(abs(curvature-2*(1-q)))<4e-6
    rows.append(dict(q=q,lambda_value=lam,gamma=g,c_U=-d*g,new_direction_norm=residual,
                     first_Gram_eigenvalues=np.linalg.eigvalsh(gram).tolist(),
                     time_rows=[dict(t=float(t),correlator=float(v),cutoff_difference=float(abs(x-v)),
                                     analytic_L120_bound_log10=bound,relative_cutoff_difference=float(abs(x-v)/v))
                                for t,x,v,bound in zip(times,a,b,bounds)],
                     finite_difference_commutator_norm=curvature.tolist(),expected_commutator_norm=2*(1-q)))

lambdas=np.geomspace(.001,3,161)
residuals=1-lambdas/np.expm1(2*lambdas)
fig,axs=plt.subplots(1,2,figsize=(10,3.5),layout='constrained')
axs[0].semilogx(lambdas,residuals,color='#3465a4',lw=2)
axs[0].axhline(.5,color='.5',ls=':')
axs[0].set(xlabel='lambda',ylabel='Squared norm of the new direction',ylim=(.45,1.03))
ts=np.linspace(0,6,241)
for q,color in zip((.0625,.25,.64),('#3465a4','#c06c21','#2a8b66')):
    C,_=calc(q,120,ts);axs[1].plot(ts,C,label=f'q={q}',color=color,lw=2)
axs[1].set(xlabel='Dimensionless relative time t',ylabel='Stationary covariant correlator')
axs[1].legend(frameon=False)
for ax in axs:ax.grid(alpha=.2)
fig.savefig(FIG,dpi=180);plt.close(fig)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
inputs=[Path(__file__),DIR/'word-checks.json',DIR/'site-checks.json']
result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),rows=rows,
            parameters=dict(q=[.0625,.25,.64],physical_relation='lambda=-log(q)/4',Jacobi_sizes=[80,120],times=times.tolist(),figure_lambda_range=[.001,3]),
            environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__,dtype='float64'),
            input_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in inputs},figure_sha256={FIG.relative_to(ROOT).as_posix():sha(FIG)},
            scope='Physical-curve finite diagnostics; analytic truncation bound is exact-arithmetic only, eigensolver roundoff not enclosed. None at t=0 denotes an exactly zero bound. q-Mehler spectrum is an established formula applied to this microscopic channel.')
with OUT.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
print(json.dumps(dict(status='passed',physical_cases=3,time_rows=18,sha256=sha(OUT))))
