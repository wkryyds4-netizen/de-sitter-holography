"""q-Jacobi and independent q-Mehler evaluations of the singlet two-point kernel."""
from pathlib import Path
from math import exp,sqrt,log,lgamma,pi
import json,hashlib,platform
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-singlet-kernel'

def truncation_size(lam,tmax,tol=1e-12):
    q=exp(-4*lam);a=2*tmax/sqrt(lam*(1-q))
    # Both Jacobi matrices have norm <= Emax. Their first L Taylor
    # coefficients on |0> agree. Bound the two exponential tails, then
    # the difference of bounded-vertex expectations. Floating eigensolver
    # rounding is a separate empirical check, not covered by this bound.
    L=max(12,int(a)+2)
    while True:
        logarithm=log(4)+L*log(a)-lgamma(L+1)-log(1-a/(L+1)) if a else -float('inf')
        if logarithm<log(tol):return L,exp(logarithm)
        L+=1

def jacobi(lam,times,L):
    q=exp(-4*lam);n=np.arange(1,L)
    off=np.sqrt((-np.expm1(-4*lam*n))/(-np.expm1(-4*lam))/lam)
    T=np.diag(off,1)+np.diag(off,-1)
    e,u=np.linalg.eigh(T)
    psi=u@(u[0,:,None]*np.exp(-1j*e[:,None]*np.asarray(times)[None,:]))
    values=np.sum(np.exp(-8*lam*np.arange(L))[:,None]*abs(psi)**2,axis=0)
    return values,float(np.max(abs(np.sum(abs(psi)**2,axis=0)-1)))

def mehler(lam,times,nodes):
    q=exp(-4*lam);r=q*q
    x,w=np.polynomial.legendre.leggauss(nodes)
    th=(x+1)*pi/2;w=w*pi/2
    count=1
    while q**count/(1-q)>1e-16:count+=1
    rho=np.ones(nodes)/(2*pi)
    numerator=1.
    denom=np.ones((nodes,nodes))
    diff=th[:,None]-th[None,:];total=th[:,None]+th[None,:]
    for k in range(count):
        z=q**k
        rho*=abs(1-z*np.exp(2j*th))**2*(1-q*z)
        numerator*=1-r*r*z
        denom*=(1-2*r*z*np.cos(diff)+(r*z)**2)*(1-2*r*z*np.cos(total)+(r*z)**2)
    weights=(w*rho)[:,None]*(w*rho)[None,:]*numerator/denom
    e=2*np.cos(th)/sqrt(lam*(1-q));gap=e[:,None]-e[None,:]
    vals=[float(np.sum(weights*np.cos(t*gap))) for t in times]
    moments=[float(np.sum(weights*gap**k)) for k in [0,2,4,6]]
    return dict(nodes=nodes,product_terms=count,values=vals,moments=moments,
                product_tail_scale=q**count/(1-q))

def main():
    OUT.mkdir(exist_ok=True);target=OUT/'kernel-numerics.json'
    if target.exists():raise FileExistsError(target)
    times=np.linspace(0,1.5,61);rows=[]
    plt.figure(figsize=(7,4.2))
    for lam in [.25,.1,.05,.02,.01]:
        L,bound=truncation_size(lam,float(times[-1]))
        vals,norm=jacobi(lam,times,L)
        finer,_=jacobi(lam,times,L+16)
        diff=float(np.max(abs(vals-finer)))
        assert norm<1e-12 and diff<1e-11
        q=exp(-4*lam)
        WB=4*(1-q)**2/(1+q);WD=4*(1/lam+4*q/(1+q))
        rows.append(dict(lambda_value=lam,L=L,analytic_truncation_bound=bound,
                         eigensolver_norm_residual=norm,L_plus_16_difference=diff,
                         static_weight=WB,dynamic_weight=WD,static_fraction=WB/(WB+WD),
                         values=vals.tolist(),max_difference_from_sech4=float(np.max(abs(vals-1/np.cosh(2*times)**4)))))
        plt.plot(times,vals,label=rf'$\lambda={lam}$')
        print(f'lambda={lam} L={L} refinement={diff:.3g}',flush=True)
    plt.plot(times,1/np.cosh(2*times)**4,'k--',lw=1.4,label=r'$\mathrm{sech}^4(2\mathcal{J}_s t)$')
    plt.xlabel(r'$\mathcal{J}_s t$');plt.ylabel(r'$C_{\perp}(t)/C_{\perp}(0)$')
    plt.legend(fontsize=8);plt.grid(alpha=.2);plt.tight_layout()
    plt.savefig(OUT/'correlator.png',dpi=180);plt.close()
    testtimes=[0.,.2,.5,1.]
    checks=[]
    for lam in [.25,.1]:
        L,_=truncation_size(lam,1.);j,_=jacobi(lam,testtimes,L)
        coarse=mehler(lam,testtimes,96);fine=mehler(lam,testtimes,192)
        err=float(np.max(abs(j-np.array(fine['values']))))
        q=exp(-4*lam);m2=2/lam*(1-q*q)
        m4=m2/lam*(1-q)*(5+6*q+3*q*q)
        m6=m2/lam**2*(1-q)**2*(10*q**5+40*q**4+85*q**3+115*q*q+91*q+35)
        exact=np.array([1.,m2,m4,m6])
        relative=float(np.max(abs(np.array(fine['moments'])/exact-1)))
        assert err<1e-10 and relative<1e-9
        checks.append(dict(lambda_value=lam,coarse=coarse,fine=fine,jacobi=j.tolist(),
                           absolute_difference=err,moment_relative_residual=relative))
    data=dict(status='passed',units='Js=1; string Hamiltonian; fixed-positive-lambda double-scaled kernel',
              times=times.tolist(),rows=rows,independent_mehler=checks,
              precision='IEEE float64; analytic Taylor bound controls Jacobi truncation only. Refinement and quadrature residuals do not certify floating-point roundoff.',
              environment=dict(python=platform.python_version(),numpy=np.__version__,matplotlib=matplotlib.__version__),
              input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['code/run_singlet_kernel.py']})
    with target.open('x') as f:json.dump(data,f,indent=2);f.write('\n')
    print('Jacobi/Mehler values and moments passed',flush=True)

if __name__=='__main__':main()
