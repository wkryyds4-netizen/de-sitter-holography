"""Gaussian F curve and exponential-packet Q corner; no interval certification."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,sys,platform
import numpy as np
import scipy
from scipy.special import loggamma,roots_legendre,erf,erfc,gamma
import mpmath
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from logtime_limit import pair_coefficients
ROOT=Path(__file__).resolve().parents[1]
W_BREAKS=[-12,-6,-3,-1,0,1,3,6,12]
CHI_BREAKS=[-12,-8,-4,-2,-1,-.5,0,.5,1,2,4,8,12]
NU=.6;SIGMA=1.;BETA=1.

def grid(breaks,n):
    a,w=roots_legendre(n)
    return (np.concatenate([(lo+hi)/2+(hi-lo)*a/2 for lo,hi in zip(breaks[:-1],breaks[1:])]),
            np.concatenate([w*(hi-lo)/2 for lo,hi in zip(breaks[:-1],breaks[1:])]))
def f(w,chi,nu=NU,sigma=SIGMA):
    w=np.asarray(w);chi=np.asarray(chi,dtype=complex);a=.25+1j*chi/4
    lg=sum(loggamma(a+sg*1j*(w+ep*nu)/2) for sg in [-1,1] for ep in [-1,1])
    return np.exp(lg-np.pi*w-sigma*sigma*w*w)
def J(chi,n):
    w,ww=grid(W_BREAKS,n);chi=np.atleast_1d(chi).astype(complex)
    return f(w[None,:],chi[:,None])@ww
def amplitude(chi,n):
    chi=np.asarray(chi,dtype=complex)
    return np.exp(1j*chi*np.log(2)+loggamma(-1j*chi/2)-loggamma((1+1j*chi)/2)-SIGMA**2*chi**2/2)*J(chi,n)**2/(1024*np.pi**6.5)
def corner_difference(h,n):
    a,w=roots_legendre(n);a=(a+1)/2;w=w/2
    rho,s,eta=np.meshgrid(a,a,a,indexing='ij')
    wr,ws,we=np.meshgrid(w,w,w,indexing='ij')
    t=(1-s)*eta
    u=h*((rho-1)*s-rho+.5);v=h*((rho-1)*t-rho+.5)
    integrand=BETA*np.exp(-BETA*h*rho)*(rho-1)**2*(1-s)*f(u,-h)*f(v,-h)
    return h**3*np.sum(wr*ws*we*integrand)
def cj(z):return [float(np.real(z)),float(np.imag(z))]

def main():
    checks=[]
    def check(name,error,tol):checks.append({'name':name,'absolute_error':float(error),'tolerance':tol,'passed':bool(error<tol)})
    c=mpmath.mp.clone();c.dps=45;nu=c.mpf('.6')
    def fm(w,ch):
        a=c.mpf('.25')+c.j*ch/4
        return c.exp(-c.pi*w-w*w)*c.fprod(c.gamma(a+sg*c.j*(w+ep*nu)/2) for sg in [-1,1] for ep in [-1,1])
    identities=[]
    for off in [['.2','.5','.4','0'],['-.3','.2','.5','0'],['1.2','-.1','.4','0']]:
        x,y,z,d=map(c.mpf,off);chi=x-y+z;_,vm=pair_coefficients(c,off,nu)
        labs=[x,y,z,d]
        H=c.exp(-c.pi*y)/(1024*c.pi**5)*c.fprod(c.gamma(c.mpf('.25')+c.j*(labs[j]-labs[(j+1)%4]+ep*nu)/2) for j in range(4) for ep in [-1,1])
        Q=x*x+(y-x)**2+(y-z)**2+z*z
        direct=H*vm/c.pi**3*c.exp(-Q/2)
        reduced=2**(2*c.j*chi)*c.cosh(c.pi*chi/2)*c.exp(-chi*chi/2)*fm(x-chi/2,chi)*fm(z-chi/2,chi)/(512*c.pi**7)
        check('factorization-'+str(len(identities)),float(abs(direct-reduced)),1e-40)
        identities.append({'offsets':off,'direct':[str(c.re(direct)),str(c.im(direct))],'reduced':[str(c.re(reduced)),str(c.im(reduced))]})
    for ch in [c.mpf('.3'),c.mpf('.7')-c.j,c.mpf('-.2')-3*c.j]:
        direct=c.gamma(-c.j*ch)*c.cosh(c.pi*ch/2)
        simplified=2**(-c.j*ch-1)*c.sqrt(c.pi)*c.gamma(-c.j*ch/2)/c.gamma((1+c.j*ch)/2)
        check('even-pole-identity-'+str(ch),float(abs(direct-simplified)),1e-40)
    j0=c.quad(lambda w:fm(w,0),[-c.inf,-2,0,2,c.inf]);j2=c.quad(lambda w:fm(w,-2*c.j),[-c.inf,-2,0,2,c.inf])
    nref=(j0/(16*c.pi**3))**2;cref=c.exp(2)*j2*j2/(32*c.pi**6)
    qref=fm(0,0)**2/(768*c.pi**7)
    Tvals=[-8,-6,-4,-2,0,2,4,6];runs={};dense=np.linspace(-8,8,161)
    for n in [24,40,64]:
        jj0=J([0],n)[0];jj2=J([-2j],n)[0];N=(jj0/(16*np.pi**3))**2;C=np.exp(2)*jj2**2/(32*np.pi**6)
        x,w=grid(CHI_BREAKS,n);amp=amplitude(x,n)
        reg=amp-1j*N*np.exp(-x*x)/(2*np.pi*x)
        shifted=amplitude(x-3j,n)
        rows=[]
        for T in Tvals:
            X=N/2*(1-erf(T/2))+np.sum(w*np.exp(1j*T*x)*reg)
            remainder=np.sum(w*np.exp(1j*T*(x-3j))*shifted)
            Xcontour=N-C*np.exp(2*T)+remainder
            rows.append({'T':T,'X':cj(X),'D':cj(2*(N-X)),'contour_X':cj(Xcontour),'contour_absolute_difference':abs(X-Xcontour),'contour_cancellation_scale':float(abs(N)+abs(C*np.exp(2*T))+np.sum(abs(w*np.exp(1j*T*(x-3j))*shifted))),
                'early_ratio':cj((N-X)/(C*np.exp(2*T)))})
        runs[n]={'N':cj(N),'C':cj(C),'J0':cj(jj0),'Jminus2i':cj(jj2),'rows':rows}
        if n==64:
            curve=N/2*(1-erf(dense/2))+np.exp(1j*dense[:,None]*x[None,:])@(w*reg)
        print('fixed curve order',n,'N',N,'C',C,flush=True)
    for name,ref in [('N',nref),('C',cref)]:
        z=complex(*runs[64][name]);check('mp-reference-'+name,abs(z-complex(ref)),2e-13)
    for name in ['N','C']:check('order-stability-'+name,abs(complex(*runs[64][name])-complex(*runs[40][name])),2e-13)
    for lo,hi in zip(runs[40]['rows'],runs[64]['rows']):
        T=hi['T'];check('real-axis-order-T'+str(T),abs(complex(*lo['X'])-complex(*hi['X'])),2e-12)
        check('real-axis-reality-T'+str(T),abs(hi['X'][1]),2e-12)
        if T<=0:check('contour-T'+str(T),hi['contour_absolute_difference'],2e-12)
    corners=[];f00=f(0,0);expected=BETA*f00*f00/6
    for h in [.2,.1,.05,.02,.01,.005]:
        low=corner_difference(h,12);high=corner_difference(h,20)
        check('corner-order-h'+str(h),abs(low-high)/h**3,2e-10)
        corners.append({'h':h,'DeltaJ':cj(high),'scaled_cubic_ratio':cj(high/(h**3*expected)),'scaled_quadrature_difference':abs(low-high)/h**3})
    # Analytic expressions evaluated in double precision, not outward-rounded intervals.
    gauss_mass=np.exp(np.pi**2/4)*np.sqrt(np.pi)
    jbound=gamma(.25)**4*gauss_mass
    jtail=gamma(.25)**4*gauss_mass/2*(erfc(12+np.pi/2)+erfc(12-np.pi/2))
    outer_tail=jbound*jbound/(512*np.pi**7)*np.sqrt(np.pi/(24)*1/np.tanh(6*np.pi))*np.sqrt(2*np.pi)*erfc(12/np.sqrt(2))
    dest=ROOT/'results/stage3/preparation';dest.mkdir(parents=True,exist_ok=True)
    fig,axes=plt.subplots(1,2,figsize=(10.8,4.0),layout='constrained')
    N=complex(*runs[64]['N']).real;C=complex(*runs[64]['C']).real
    axes[0].plot(dense,curve.real/N,color='#275d8c',lw=2);axes[0].set(xlabel='T = t - log(2M)',ylabel='X_F(T) / N',title='Fixed threshold, Gaussian preparation',ylim=(-.02,1.04))
    mask=dense<=-1
    axes[1].semilogy(dense[mask],2*(N-curve.real[mask]),label='Limiting commutator',color='#275d8c',lw=2)
    axes[1].semilogy(dense[mask],2*C*np.exp(2*dense[mask]),'--',label='2 C exp(2T)',color='#bd6530')
    axes[1].set(xlabel='T',ylabel='D_F(T)',title='Early ordered-limit asymptotic');axes[1].legend(frameon=False)
    for ax in axes:ax.grid(alpha=.18)
    fig.savefig(dest/'fixed-curve.pdf');fig.savefig(dest/'fixed-curve.png',dpi=180);plt.close(fig)
    paths=['code/preparation_asymptotics.py','code/logtime_limit.py','code/crossing_assembly.py','manuscript/preparation.tex','manuscript/logtime.tex','sources/kolchmeyer-liu-2411.08090v1.html','sources/dlmf-5.11.html']
    out={'created_at_utc':datetime.now(timezone.utc).isoformat(),'status':'passed' if all(r['passed'] for r in checks) else 'failed','checks':checks,
        'factorization_rows':identities,'fixed_runs':runs,'corner_rows':corners,
        'mp_references':{'N':str(c.re(nref)),'C_F':str(c.re(cref)),'Q_commutator_cubic_beta1':str(c.re(qref)),'working_dps':45},
        'analytic_tail_expressions_evaluated':{'J_w_tail_at_12':float(jtail),'real_axis_A_outer_tail_at_12':float(outer_tail),'scope':'bounds for J and A separately; evaluated not outward rounded, not a total certified quadrature error'},
        'parameters':{'nu':NU,'sigma':SIGMA,'beta':BETA,'w_breaks':W_BREAKS,'chi_breaks':CHI_BREAKS,'orders':[24,40,64],'corner_cube_orders':[12,20],'contour_shift':3,'T_values':Tvals},
        'dense_fixed_curve':{'T':dense.tolist(),'X':[cj(v) for v in curve]},
        'versions':{'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'scipy':scipy.__version__,'mpmath':mpmath.__version__,'matplotlib':matplotlib.__version__},
        'limitations':['empirical quadrature; no interval certification','full Q time curve not computed; Q diagnostics test its local corner coefficient','lower-contour evaluation becomes cancellation-sensitive at positive T, retained but not used as primary values','no finite-M or finite-G claim; no original-priority claim'],
        'hashes':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}}
    (dest/'diagnostics.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(out['status'],len(checks),'checks','tails',jtail,outer_tail,flush=True)
    for r in corners:print('corner',r['h'],r['scaled_cubic_ratio'],flush=True)
    if out['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
