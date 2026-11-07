"""Time-dependent collapsed connected four-point function, with checks."""
from pathlib import Path
from math import exp,sqrt,pi
import hashlib,json,platform
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from singlet_multimark import probe_moment
from run_singlet_kernel import truncation_size

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-singlet-multimark'


def evolve(lam,times,L):
    n=np.arange(1,L);q=exp(-4*lam)
    off=np.sqrt(-np.expm1(-4*lam*n)/(-np.expm1(-4*lam)*lam))
    e,u=np.linalg.eigh(np.diag(off,1)+np.diag(off,-1))
    psi=u@(u[0,:,None]*np.exp(-1j*e[:,None]*np.asarray(times)[None,:]))
    probs=abs(psi)**2
    g=[np.sum(np.exp(-4*lam*d*np.arange(L))[:,None]*probs,axis=0) for d in [2,4]]
    connected=(1+q**4)*g[1]-2*g[0]**2
    return np.array(g),connected,float(np.max(abs(probs.sum(axis=0)-1)))


def integral(lam,times,dimension,nodes):
    q=exp(-4*lam);r=q**dimension
    x,w=np.polynomial.legendre.leggauss(nodes);th=(x+1)*pi/2;w*=pi/2
    count=1
    while q**count/(1-q)>1e-16:count+=1
    rho=np.ones(nodes)/(2*pi);den=np.ones((nodes,nodes));num=1.
    diff=th[:,None]-th[None,:];total=th[:,None]+th[None,:]
    for k in range(count):
        z=q**k
        rho*=abs(1-z*np.exp(2j*th))**2*(1-q*z);num*=1-r*r*z
        den*=(1-2*r*z*np.cos(diff)+(r*z)**2)*(1-2*r*z*np.cos(total)+(r*z)**2)
    weight=(w*rho)[:,None]*(w*rho)[None,:]*num/den
    e=2*np.cos(th)/sqrt(lam*(1-q));gap=e[:,None]-e[None,:]
    return dict(values=[float(np.sum(weight*np.cos(t*gap))) for t in times],
        moments=[float(np.sum(weight*gap**k)) for k in [0,2]],nodes=nodes,product_terms=count,
        product_tail_scale=q**count/(1-q))


def exact_collapsed():
    from collections import Counter
    from singlet_coefficient import pairings
    from singlet_multimark import crossing
    from itertools import combinations
    records=[]
    for n in [0,2,4,6,8]:
        for a in range(n+1):
            b=n-a
            actual=probe_moment('XX'+'H'*a+'XX'+'H'*b)
            out=Counter()
            for p in pairings(tuple(range(n))):
                cr=sum(crossing(x,y) for x,y in combinations(p,2))
                transfer=4*sum(x<a<=y for x,y in p)
                out[0,cr]+=1
                out[0,cr+transfer]+=1
                out[0,cr+transfer+4]+=1
            assert actual==dict(out),(a,b,actual,out)
            records.append(dict(a=a,b=b))
    return records


def main():
    target=OUT/'fourpoint-numerics.json'
    if target.exists():raise FileExistsError(target)
    exact=exact_collapsed();times=np.linspace(0,1.5,61);rows=[]
    fig,axes=plt.subplots(1,2,figsize=(9,3.7))
    for lam in [.25,.1,.05,.02,.01]:
        L,bound=truncation_size(lam,float(times[-1]),1e-13)
        g,k,norm=evolve(lam,times,L);gf,kf,_=evolve(lam,times,L+16)
        delta=abs(k-kf);scale=(1+exp(-16*lam))*abs(g[1])+2*g[0]**2
        q=exp(-4*lam);k0=np.expm1(-16*lam)
        assert abs(k[0]-k0)<1e-12 and max(delta)<1e-11 and norm<1e-12
        rows.append(dict(lambda_value=lam,L=L,kappa=k.tolist(),g2=g[0].tolist(),g4=g[1].tolist(),
            analytic_each_kernel_truncation_bound=bound,analytic_connected_truncation_bound=6*bound,
            refinement_max_absolute=float(max(delta)),refinement_max_relative=float(max(delta/np.maximum(abs(k),1e-300))),
            max_cancellation_condition=float(max(scale/np.maximum(abs(k),1e-300))),norm_residual=norm,
            exact_kappa_zero=float(k0),exact_curvature=2/lam*(3-4*q*q+q**8)))
        axes[0].plot(times,k,label=rf'$\lambda={lam}$');axes[1].plot(times,k/lam,label=rf'$\lambda={lam}$')
        print('lambda',lam,'L',L,'kappa(.5)',float(k[20]),'refinement',float(max(delta)),flush=True)
    for ax in axes:
        ax.set_xlabel(r'$\mathcal{J}_s t$');ax.grid(alpha=.2);ax.axhline(0,color='k',lw=.6)
    axes[0].set_ylabel(r'$\kappa_4(t)$');axes[1].set_ylabel(r'$\kappa_4(t)/\lambda$')
    axes[0].legend(fontsize=8);fig.tight_layout();fig.savefig(OUT/'connected-fourpoint.png',dpi=180);plt.close(fig)
    tests=[];testtimes=[0.,.2,.5,1.]
    for lam in [.25,.1]:
        L,_=truncation_size(lam,1.);g,k,_=evolve(lam,testtimes,L);q=exp(-4*lam)
        quadratures=[]
        for d in [2,4]:
            coarse=integral(lam,testtimes,d,96);fine=integral(lam,testtimes,d,192)
            diff=float(max(abs(np.array(fine['values'])-g[[2,4].index(d)])))
            err=abs(fine['moments'][1]/(2/lam*(1-q**d))-1)
            assert diff<1e-10 and err<1e-10
            quadratures.append(dict(dimension=d,coarse=coarse,fine=fine,value_error=diff,m2_relative_error=err))
        ki=(1+q**4)*np.array(quadratures[1]['fine']['values'])-2*np.array(quadratures[0]['fine']['values'])**2
        tests.append(dict(lambda_value=lam,times=testtimes,quadratures=quadratures,
            connected_max_absolute_error=float(max(abs(ki-k))),connected_max_relative_error=float(max(abs(ki-k)/abs(k)))))
    inputs=['code/run_singlet_fourpoint.py','code/singlet_multimark.py','code/singlet_coefficient.py','code/run_singlet_kernel.py']
    data=dict(status='passed',units='Js=1, beta=mu=0, fixed-positive-lambda double scaling before varying lambda',
        observable='kappa4(t)=<X(t)^2 X(0)^2>-1-2<X(t)X(0)>^2; not OTOC',
        exact_collapsed_moments=exact,times=times.tolist(),rows=rows,independent_quadrature=tests,
        error_scope='Taylor bound controls Jacobi truncation. IEEE float64 rounding and product/integration errors are empirical, not interval certified. No long-time extrapolation or uniform small-lambda remainder claimed.',
        environment=dict(python=platform.python_version(),numpy=np.__version__,matplotlib=matplotlib.__version__),
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    target.write_text(json.dumps(data,indent=2)+'\n');print('passed',flush=True)

if __name__=='__main__':main()
