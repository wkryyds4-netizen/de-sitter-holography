"""Evaluate the proven q^3 kernel by Jacobi propagation and Mehler quadrature."""
from pathlib import Path
import json,hashlib,platform
from math import exp,sqrt,log,lgamma,pi
import numpy as np
import sympy as s
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-cubic-kernel'
def cutoff(lam,tmax,tol=1e-13):
    q=exp(-4*lam);a=2*tmax/sqrt(lam*(1-q));n=max(12,int(a)+2)
    while log(4)+n*log(a)-lgamma(n+1)-log(1-a/(n+1))>log(tol):n+=1
    return n,4*exp(n*log(a)-lgamma(n+1))/(1-a/(n+1))
def jacobi(lam,times,n):
    q=exp(-4*lam);k=np.arange(1,n)
    off=np.sqrt(-np.expm1(-4*lam*k)/(-np.expm1(-4*lam)*lam))
    mat=np.diag(off,1)+np.diag(off,-1);energy,v=np.linalg.eigh(mat)
    psi=v@(v[0,:,None]*np.exp(-1j*energy[:,None]*np.asarray(times)[None,:]))
    values=np.sum(np.exp(-12*lam*np.arange(n))[:,None]*abs(psi)**2,axis=0)
    return values,float(np.max(abs(np.sum(abs(psi)**2,axis=0)-1)))
def mehler(lam,times,nodes):
    q=exp(-4*lam);r=q**3;x,w=np.polynomial.legendre.leggauss(nodes)
    th=(x+1)*pi/2;w=w*pi/2;terms=1
    while q**terms/(1-q)>1e-17:terms+=1
    rho=np.ones(nodes)/(2*pi);numerator=1.
    denominator=np.ones((nodes,nodes))
    difference=th[:,None]-th[None,:];total=th[:,None]+th[None,:]
    for k in range(terms):
        z=q**k;rho*=abs(1-z*np.exp(2j*th))**2*(1-q*z)
        numerator*=1-r*r*z
        denominator*=(1-2*r*z*np.cos(difference)+(r*z)**2)*(1-2*r*z*np.cos(total)+(r*z)**2)
    measure=(w*rho)[:,None]*(w*rho)[None,:]*numerator/denominator
    energy=2*np.cos(th)/sqrt(lam*(1-q));omega=energy[:,None]-energy[None,:]
    return dict(nodes=nodes,terms=terms,values=[float(np.sum(measure*np.cos(t*omega))) for t in times],
      moments=[float(np.sum(measure*omega**k)) for k in [0,2,4]],product_tail=q**terms/(1-q))
def main():
    target=OUT/'kernel.json'
    if target.exists():raise FileExistsError(target)
    times=np.linspace(0,1,201);limit=1/np.cosh(2*times)**6;records=[];comparisons=[]
    for lam in [.25,.1,.025,.01,.005]:
        n,tail=cutoff(lam,1);vals,norm=jacobi(lam,times,n);ref,norm2=jacobi(lam,times,n+16)
        absolute=float(np.max(abs(vals-ref)));relative=float(np.max(abs(vals-ref)/np.maximum(abs(ref),1e-300)))
        assert absolute<2e-12 and relative<2e-10
        records.append(dict(lam=lam,cutoff=n,refined_cutoff=n+16,analytic_truncation_bound=tail,
          max_absolute_refinement=absolute,max_relative_refinement=relative,norm_error=max(norm,norm2),
          values=list(map(float,ref)),small_lambda_max_absolute_deviation=float(np.max(abs(ref-limit))),
          small_lambda_max_relative_deviation=float(np.max(abs(ref-limit)/limit))))
        print('kernel',lam,n,'refinement',absolute,'limit-relative',records[-1]['small_lambda_max_relative_deviation'],flush=True)
        if lam in [.25,.1]:
            test_times=[0.,.25,.5,1.]
            a=mehler(lam,test_times,128);b=mehler(lam,test_times,192);j,_=jacobi(lam,test_times,n+16)
            err=float(np.max(abs(j-b['values'])));refine=float(np.max(abs(np.array(a['values'])-b['values'])))
            q=exp(-4*lam)
            exact=[1,2*(1-q**3)/lam,2*(1-q)**2*(1+q+q*q)*(3*q**3+6*q*q+6*q+5)/lam**2]
            merror=float(np.max(abs(np.array(b['moments'])/exact-1)))
            assert err<1e-11 and refine<1e-11 and merror<1e-11
            comparisons.append(dict(lam=lam,times=test_times,coarse=a,fine=b,jacobi=list(map(float,j)),
              max_absolute_difference=err,quadrature_refinement=refine,relative_moment_error=merror))
    t=s.Symbol('t',real=True)
    f=s.sech(2*t)**6
    assert s.diff(f,t,2).subs(t,0)==-24 and s.diff(f,t,4).subs(t,0)==1920
    fig,ax=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
    for row in records:
        ax[0].plot(times,row['values'],label=fr"$\lambda={row['lam']}$")
        ax[1].plot(times,np.array(row['values'])/limit-1)
    ax[0].plot(times,limit,'k--',label=r"$\mathrm{sech}^{6}(2t)$")
    ax[0].set(xlabel=r"$\mathcal{J}_s t$",ylabel="Normalized correlation",yscale='log')
    ax[1].set(xlabel=r"$\mathcal{J}_s t$",ylabel="Relative deviation from small-coupling limit")
    ax[0].legend(fontsize=8);ax[0].grid(alpha=.2);ax[1].grid(alpha=.2)
    fig.savefig(OUT/'correlation.png',dpi=180);plt.close(fig)
    data=dict(status='passed',units='Js=1; full beta=mu=0 trace; DS then fixed real time',times=list(map(float,times)),
      records=records,independent_mehler=comparisons,
      error_scope='Jacobi operator/Taylor truncation bound analytic; eigensolver roundoff and Mehler quadrature checked empirically by refinement, not interval-certified',
      small_lambda_fourier='S3(omega)=pi*omega*(omega^2+16Js^2)*(omega^2+64Js^2)/(7680Js^6*sinh(pi*omega/(4Js))); C=integral S3 exp(it omega) domega/(2pi)',
      environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__,matplotlib=matplotlib.__version__),
      input_sha256={'code/run_cubic_kernel.py':hashlib.sha256((ROOT/'code/run_cubic_kernel.py').read_bytes()).hexdigest()})
    target.write_text(json.dumps(data,indent=2)+'\n')
if __name__=='__main__':main()
