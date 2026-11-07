"""Evaluate the q probe kernel and independently check its d=1 limit."""
from pathlib import Path
import json,hashlib,platform
from math import exp,sqrt,pi
import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_cubic_kernel import cutoff
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-linear-kernel'
def jacobi(lam,times,n):
    k=np.arange(1,n);off=np.sqrt(-np.expm1(-4*lam*k)/(-np.expm1(-4*lam)*lam))
    E,V=np.linalg.eigh(np.diag(off,1)+np.diag(off,-1))
    psi=V@(V[0,:,None]*np.exp(-1j*E[:,None]*np.array(times)[None,:]))
    return np.sum(np.exp(-4*lam*np.arange(n))[:,None]*abs(psi)**2,axis=0),float(np.max(abs(np.sum(abs(psi)**2,axis=0)-1)))
def mehler(lam,times,nodes):
    q=exp(-4*lam);r=q;x,w=np.polynomial.legendre.leggauss(nodes)
    th=(x+1)*pi/2;w=w*pi/2;terms=1
    while q**terms/(1-q)>1e-17:terms+=1
    rho=np.ones(nodes)/(2*pi);num=1.;den=np.ones((nodes,nodes))
    diff=th[:,None]-th[None,:];total=th[:,None]+th[None,:]
    for k in range(terms):
        z=q**k;rho*=abs(1-z*np.exp(2j*th))**2*(1-q*z);num*=1-r*r*z
        den*=(1-2*r*z*np.cos(diff)+(r*z)**2)*(1-2*r*z*np.cos(total)+(r*z)**2)
    measure=(w*rho)[:,None]*(w*rho)[None,:]*num/den
    E=2*np.cos(th)/sqrt(lam*(1-q));om=E[:,None]-E[None,:]
    return dict(nodes=nodes,terms=terms,values=[float(np.sum(measure*np.cos(t*om))) for t in times],
      moments=[float(np.sum(measure*om**k)) for k in [0,2,4]],product_tail=q**terms/(1-q))
def main():
    OUT.mkdir(parents=True,exist_ok=True);path=OUT/'kernel.json'
    if path.exists():raise FileExistsError(path)
    times=np.linspace(0,1,201);limit=1/np.cosh(2*times)**2;rows=[];independent=[]
    for lam in [.25,.1,.025,.01,.005]:
        n,tail=cutoff(lam,1);a,an=jacobi(lam,times,n);b,bn=jacobi(lam,times,n+16)
        ae=float(np.max(abs(a-b)));re=float(np.max(abs(a-b)/abs(b)))
        assert ae<2e-12 and re<2e-10
        rows.append(dict(lam=lam,cutoff=n,refined_cutoff=n+16,analytic_truncation_bound=tail,
          max_absolute_refinement=ae,max_relative_refinement=re,norm_error=max(an,bn),values=list(map(float,b)),
          max_small_lambda_relative_deviation=float(np.max(abs(b-limit)/limit))))
        if lam in [.25,.1]:
            tt=[0,.25,.5,1];coarse=mehler(lam,tt,128);fine=mehler(lam,tt,192);j,_=jacobi(lam,tt,n+16)
            q=exp(-4*lam);exact=[1,2*(1-q)/lam,2*(1-q)**2*(3*q+5)/lam**2]
            err=float(np.max(abs(j-fine['values'])));qe=float(np.max(abs(np.array(coarse['values'])-fine['values'])))
            me=float(np.max(abs(np.array(fine['moments'])/exact-1)))
            assert max(err,qe,me)<1e-10
            independent.append(dict(lam=lam,times=tt,coarse=coarse,fine=fine,jacobi=list(map(float,j)),max_absolute_difference=err,quadrature_refinement=qe,relative_moment_error=me))
        print('kernel',lam,n,'relative limit difference',rows[-1]['max_small_lambda_relative_deviation'],flush=True)
    fourier=[]
    for omega in [0,1,4,10]:
        vals=[]
        for dps in [50,80]:
            mp.mp.dps=dps
            val=2*mp.quad(lambda t:mp.cos(omega*t)/mp.cosh(2*t)**2,[0,1,2,4,8,16,32,mp.inf])
            exact=mp.mpf(1) if omega==0 else mp.pi*omega/(4*mp.sinh(mp.pi*omega/4))
            assert abs(val/exact-1)<mp.mpf('1e-40');vals.append(mp.nstr(val,48))
        change=abs(mp.mpf(vals[0])/mp.mpf(vals[1])-1);assert change<mp.mpf('1e-40')
        fourier.append(dict(omega=omega,values=vals,relative_precision_change=mp.nstr(change,20)))
    fig,ax=plt.subplots(1,2,figsize=(10,3.6),layout='constrained')
    for row in rows:
        ax[0].plot(times,row['values'],label=fr"$\lambda={row['lam']}$")
        ax[1].plot(times,np.array(row['values'])/limit-1)
    ax[0].plot(times,limit,'k--',label=r'$\mathrm{sech}^{2}(2t)$')
    ax[0].set(xlabel=r'$\mathcal{J}_s t$',ylabel='Normalized F correlation',yscale='log')
    ax[1].set(xlabel=r'$\mathcal{J}_s t$',ylabel='Relative deviation from small-coupling limit')
    ax[0].legend(fontsize=8);ax[0].grid(alpha=.2);ax[1].grid(alpha=.2)
    fig.savefig(OUT/'correlation.png',dpi=180);plt.close(fig)
    names=['code/run_linear_kernel.py','code/run_cubic_kernel.py']
    path.write_text(json.dumps(dict(status='passed',units='Js=1; DS first, fixed real time',times=list(map(float,times)),
      records=rows,independent_mehler=independent,fourier=fourier,
      error_scope='Analytic Jacobi truncation bound. Float refinement and Mehler/Fourier checks are empirical, not interval-certified.',
      environment=dict(python=platform.python_version(),numpy=np.__version__,mpmath=mp.__version__,matplotlib=matplotlib.__version__),
      input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')
if __name__=='__main__':main()
