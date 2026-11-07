"""Conserved-subtracted cubic spectrum from the three diagonal Mehler measures."""
from pathlib import Path
from math import exp,sqrt,log,lgamma,pi
import json,hashlib,platform
import numpy as np
import sympy as s
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-filtered-spectrum'

def parameters(lam):
    q=exp(-4*lam);ss=1/lam;r=q*q
    a=-(1-q)**2/(1+q*q)
    D=(1+q*q)*(1+q+q*q)
    A1=4*(1-q*q)*(1+q*q)*(1+q+q*q)**2
    A2=16*q*q*(1+q)*(1+q*q)*(1+q+q*q)*(2*q*q-q+2)
    VR=ss**3*(lam*A1+lam*lam*A2)/D**2
    VF=64*q*q;W=4*(1/lam+4*q/(1+q))
    norm=VR+VF+W*ss*(a*a*(1+r)/2+4.5*(1-r))
    mu4=ss**2*(10+2*q-(16+8*q)*r+6*(1+q)*r*r)
    m2=VR*2*ss*(1-q**3)+VF*2*ss*(1-q)+W*(a*a*ss**2*(1+q)*(1-r*r)/2+2.25*mu4)
    return dict(q=q,sigma2=ss,a=a,VR=VR,VF=VF,W=W,norm=norm,mu2=m2/norm)

def cutoff(lam,tmax,tol=1e-12):
    z=parameters(lam);B=2/sqrt(lam*(1-z['q']));bt=B*tmax
    factor=(2*(z['VR']+z['VF'])+z['W']*B*B*(2*z['a']**2+18))/z['norm']
    n=max(8,int(bt)+4)
    while True:
        k=n-2
        delta=2*exp(k*log(bt)-lgamma(k+1))/(1-bt/(k+1))
        if delta*factor<tol:return n,delta*factor
        n+=1

def jacobi(lam,times,n):
    z=parameters(lam);q=z['q'];k=np.arange(1,n)
    off=np.sqrt(-np.expm1(-4*lam*k)/(-np.expm1(-4*lam)*lam))
    en,V=np.linalg.eigh(np.diag(off,1)+np.diag(off,-1))
    base=V[0,:,None]*np.exp(-1j*en[:,None]*np.asarray(times)[None,:])
    psi=V@base;one=V@(en[:,None]*base);two=V@(en[:,None]**2*base)
    pp=abs(psi)**2;nn=np.arange(n)[:,None]
    c1=np.sum(q**nn*pp,axis=0);c3=np.sum(q**(3*nn)*pp,axis=0)
    mixed=np.real(two*np.conj(psi));onesq=abs(one)**2;rw=q**(2*nn)
    jsum=np.sum(rw*(mixed+onesq),axis=0)/2
    lsum=2*np.sum(rw*(mixed-onesq),axis=0)
    pieces=np.array([z['VF']*c1,z['VR']*c3,z['W']*z['a']**2*jsum,z['W']*2.25*lsum])/z['norm']
    return np.sum(pieces,axis=0),pieces

def quadrature(lam,times,nodes):
    z=parameters(lam);q=z['q'];x,w=np.polynomial.legendre.leggauss(nodes)
    th=(x+1)*pi/2;w=w*pi/2;steps=1
    while q**steps/(1-q)>1e-18:steps+=1
    rho=np.ones(nodes)/(2*pi)
    for k in range(steps):rho*=abs(1-q**k*np.exp(2j*th))**2*(1-q**(k+1))
    E=2*np.cos(th)/sqrt(lam*(1-q));de=E[:,None]-E[None,:];su=(E[:,None]+E[None,:])/2
    measure=(rho*w)[:,None]*(rho*w)[None,:]
    def kernel(r):
        num=1.;den=np.ones((nodes,nodes))
        for k in range(steps):
            v=r*q**k;num*=1-r*r*q**k
            den*=(1-2*v*np.cos(th[:,None]-th[None,:])+v*v)*(1-2*v*np.cos(th[:,None]+th[None,:])+v*v)
        return num/den
    K1,K2,K3=[kernel(q**j) for j in [1,2,3]]
    full=measure*(z['VF']*K1+z['VR']*K3+z['W']*(z['a']*su+1.5*de)**2*K2)/z['norm']
    mixed=float(np.sum(measure*K2*su*de))
    return dict(nodes=nodes,terms=steps,values=[float(np.sum(full*np.cos(t*de))) for t in times],
        moments=[float(np.sum(full*de**j)) for j in [0,2]],energy_sum_difference_overlap=mixed,
        minimum_density_weight=float(np.min(full)))

def main():
    OUT.mkdir(parents=True,exist_ok=True);target=OUT/'spectrum.json'
    if target.exists():raise FileExistsError(target)
    t=np.linspace(0,1,201);f=1/np.cosh(2*t);limit=(47*f**6-36*f**4)/11
    rows=[];independent=[]
    for lam in [.25,.1,.025,.01,.005]:
        n,tail=cutoff(lam,1);a,_=jacobi(lam,t,n);b,parts=jacobi(lam,t,n+16)
        err=float(np.max(abs(a-b)));assert err<2e-11
        par=parameters(lam)
        assert abs(b[0]-1)<1e-11
        rows.append(dict(lam=lam,cutoff=n,refined_cutoff=n+16,analytic_normalized_truncation_bound=tail,
            max_absolute_refinement=err,parameters=par,values=list(map(float,b)),
            contributions=parts.tolist(),max_absolute_small_lambda_deviation=float(np.max(abs(b-limit))),
            equal_time_fractions=parts[:,0].tolist()))
        print('raw spectrum',lam,n,'error',err,'norm*lambda',lam*par['norm'],flush=True)
        if lam in [.25,.1]:
            tt=[0,.25,.5,1];lo=quadrature(lam,tt,128);hi=quadrature(lam,tt,192);j,_=jacobi(lam,tt,n+16)
            e=float(np.max(abs(j-hi['values'])));ref=float(np.max(abs(np.array(lo['values'])-hi['values'])))
            assert max(e,ref,abs(hi['moments'][0]-1),abs(hi['moments'][1]/par['mu2']-1))<1e-10
            assert hi['minimum_density_weight']>=0
            independent.append(dict(lam=lam,times=tt,coarse=lo,fine=hi,jacobi=list(map(float,j)),max_absolute_difference=e))
    x=s.Symbol('x',real=True);g=(47*s.sech(2*x)**6-36*s.sech(2*x)**4)/11
    assert g.subs(x,0)==1 and -s.diff(g,x,2).subs(x,0)==s.Rational(552,11)
    zero=float(np.arccosh(np.sqrt(47/36))/2)
    assert abs((47/np.cosh(2*zero)**6-36/np.cosh(2*zero)**4)/11)<1e-14
    fig,ax=plt.subplots(1,2,figsize=(10,3.7),layout='constrained')
    for row in rows:ax[0].plot(t,row['values'],label=fr"$\lambda={row['lam']}$")
    ax[0].plot(t,limit,'k--',label='Small-coupling limit');ax[0].axhline(0,color='.6',lw=.6)
    ax[0].set(xlabel=r'$\mathcal{J}_s t$',ylabel='Normalized subtracted cubic correlation');ax[0].legend(fontsize=8)
    small=rows[-1];labels=['F','R','Energy-sum descendant','Time descendant']
    for curve,label in zip(small['contributions'],labels):ax[1].plot(t,curve,label=label)
    ax[1].set(xlabel=r'$\mathcal{J}_s t$',ylabel=r'Contributions at $\lambda=0.005$');ax[1].legend(fontsize=8)
    for a in ax:a.grid(alpha=.2)
    fig.savefig(OUT/'spectrum.png',dpi=180);plt.close(fig)
    names=['code/run_filtered_spectrum.py']
    target.write_text(json.dumps(dict(status='passed',units='Js=1; DS first; U=M3-PN; correlator <U(t)U-dagger>',
        times=t.tolist(),records=rows,independent_mehler=independent,
        small_lambda=dict(normalized_correlation='(47 sech(2t)^6 - 36 sech(2t)^4)/11',
            norm_leading='176 Js^6/lambda',normalized_mu2='552 Js^2/11',first_positive_zero=zero,
            positive_spectrum='(47 (omega/Js)^2 + 128) S2(omega)/880'),
        error_scope='Taylor/Jacobi truncation analytic; eigensolver and Mehler precision empirical; relative errors not used across correlation zeros',
        environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__,matplotlib=matplotlib.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')
if __name__=='__main__':main()
