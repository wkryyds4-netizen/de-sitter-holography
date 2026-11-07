"""Calculate the target spectrum, thermal contour checks and radial examples."""
from pathlib import Path
import json, hashlib, platform
import numpy as np
import scipy
import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.linalg import expm
from rindler_singlet_packet import spectrum,target,midpoint,radial_profiles,radial_norm

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-rindler-packet'


def encode(z):
    z=complex(z)
    return [z.real,z.imag]


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    path=OUT/'study.json'
    if path.exists():raise FileExistsError(path)
    mp.mp.dps=40
    checks=[]
    def check(name,error,tol,**extra):
        checks.append(dict(name=name,error=float(abs(error)),tolerance=tol,passed=bool(abs(error)<tol),**extra))

    # Independent Fourier integral of the real-time target, including zero.
    for nu in [0,.3,1,2,4]:
        direct=2*mp.quad(lambda t: target(t)*mp.cos(nu*t),[0,2,6,12,24,48,mp.inf])
        check('Fourier transform',direct-spectrum(nu),1e-25,nu=nu)
    check('spectral variance',mp.quad(lambda n:n*n*spectrum(n)/mp.pi,[0,1,4,12,mp.inf])-1,1e-30)

    epsilons=[.5,.2,.1,.05,.01]
    kernel=[]
    for e in epsilons:
        c0=midpoint(0,e)
        ci=midpoint(mp.j*mp.pi/2,e)
        norm=mp.quad(lambda n: 2/mp.pi*n*mp.sinh(mp.pi*n)*spectrum(n)*mp.exp(-e*n*n),[0,1,4,12,24,mp.inf])
        check('KL spectral norm',norm/radial_norm(e)-1,1e-30,epsilon=e)
        times=[]
        for eta in [0,.5,1,2,4]:
            val=midpoint(eta,e)
            heat=mp.quad(lambda s: mp.exp(-s*s/(4*e))*target(eta-s)/mp.sqrt(4*mp.pi*e),[-mp.inf,-2,0,2,mp.inf])
            check('spectral vs heat convolution',val-heat,1e-25,epsilon=e,eta=eta)
            check('uniform regulator bound',max(mp.mpf(0),abs(val-target(eta))-e),1e-30,epsilon=e,eta=eta)
            times.append(dict(eta=eta,midpoint=float(val),target=float(target(eta)),absolute_difference=float(abs(val-target(eta)))))
        kernel.append(dict(epsilon=e,C0=float(c0),C_ipi2=float(mp.re(ci)),
            regulated_fourpoint_at_zero=float(mp.re(c0*c0+2*ci*ci)),normalized_fourpoint_at_zero=float(mp.re(1+2*(ci/c0)**2)),
            radial_norm_squared=float(norm),times=times))

    # A finite Gibbs trace: the same two-point dressing is tested at four points.
    H=np.diag([-.7,.3,1.1]); X=np.array([[.2,1,.3j],[1,-.4,.6],[-.3j,.6,.2]])
    X=X/np.sqrt(np.trace(X@X).real/3)
    beta=1.3; Z=np.trace(expm(-beta*H)).real/3
    R=lambda power:expm(-power*beta*H)/Z**power
    evo=lambda A,t:expm(1j*t*H)@A@expm(-1j*t*H)
    E=expm(beta*H/4);A2=np.sqrt(Z)*E@X@E
    E4=expm(beta*H/8);A4=Z**.25*E4@X@E4
    tr=lambda a:np.trace(a)/3
    matrix=[]
    for t in [0,.4,1.1]:
        XX=evo(X,t);AA=evo(A2,t);BB=evo(A4,t)
        c=tr(XX@X);c2=tr(R(.5)@AA@R(.5)@A2)
        f=tr(XX@X@XX@X)
        f2=tr(R(.25)@AA@R(.25)@A2@R(.25)@AA@R(.25)@A2)
        f4=tr(R(.25)@BB@R(.25)@A4@R(.25)@BB@R(.25)@A4)
        predicted=Z*tr(XX@E@X@E@XX@E@X@E)
        check('finite trace midpoint',c-c2,1e-12,t=t)
        check('fixed A2 fourpoint formula',f2-predicted,1e-12,t=t)
        check('n-specific A4 fourpoint',f4-f,1e-12,t=t)
        matrix.append(dict(t=t,trace_two=encode(c),midpoint_two=encode(c2),trace_four=encode(f),fixed_A2_four=encode(f2),A4_four=encode(f4)))

    # Bessel grid shared by three physical radial smears. Intermediate saved once.
    logx=np.linspace(-8,3,111)
    profiles,cond=radial_profiles(logx,[.5,.2,.1],order=128,cutoff=32,dps=25)
    np.savez(OUT/'radial-grid.npz',log_x=logx,epsilons=[.5,.2,.1],profiles=profiles,cancellation=cond)
    print('radial grid saved',flush=True)
    selected=np.array([-8,-4,-1,0,1,2,3],dtype=float)
    high,hcond=radial_profiles(selected,[.5,.2,.1],order=256,cutoff=40,dps=40)
    low,lcond=radial_profiles(selected,[.5,.2,.1],order=128,cutoff=32,dps=25)
    refinement=[]
    for i,x in enumerate(selected):
        for j,e in enumerate([.5,.2,.1]):
            check('radial refinement',high[i,j]-low[i,j],2e-7,log_x=float(x),epsilon=e)
            refinement.append(dict(log_x=float(x),epsilon=e,low=float(low[i,j]),high=float(high[i,j]),
                absolute_difference=float(abs(high[i,j]-low[i,j])),cancellation=float(hcond[i,j])))

    # Smooth spectral plot uses inexpensive float64 cosine quadrature.
    from scipy.special import roots_legendre
    n,w=roots_legendre(300);n=(n+1)*20;w=w*20
    sc=8*np.pi/3*n*(1+n*n)/np.sinh(np.pi*n)
    eta=np.linspace(0,5,150)
    fig,axs=plt.subplots(1,2,figsize=(10,3.7))
    axs[0].plot(eta,1/np.cosh(eta/2)**4,'k--',label='unregulated target')
    for e in [.5,.2,.1]:
        axs[0].plot(eta,np.cos(np.outer(eta,n))@(w*sc*np.exp(-e*n*n)/np.pi),label=fr'$\epsilon={e}$')
    axs[0].set(xlabel=r'boost time $\eta$',ylabel=r'$C_\epsilon(\eta)$',title='Midpoint correlator')
    axs[0].legend(fontsize=8)
    # One epsilon, three masses: horizontal shifts of the same radial function.
    for mass in [.5,1,2]:
        axs[1].plot(logx-np.log(mass),profiles[:,1],label=fr'$m/m_0={mass}$')
    axs[1].set(xlabel=r'$\log(m_0\rho)$',ylabel=r'$f_{m,0.2}(\rho)$',title='Same time kernel, different radial scales')
    axs[1].legend(fontsize=8)
    fig.tight_layout();fig.savefig(OUT/'rindler-packet.png',dpi=180);plt.close(fig)
    data=dict(status='passed' if all(c['passed'] for c in checks) else 'failed',
        model='1+1D free massive scalar, Minkowski vacuum, right-wedge boost beta=2pi, midpoint correlator; DS then small lambda target, eta=4 Js t as candidate normalization',
        parameters=dict(dps=40,radial_dps=25,radial_order=128,radial_cutoff=32,refinement_order=256,refinement_cutoff=40,log_x_range=[-8,3],radial_points=111),
        errors='Quadrature refinements are empirical, not interval bounds. Analytic uniform time-kernel regulator error <= epsilon. Full radial norm analytic; radial spatial tails not numerically certified.',
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__,matplotlib=matplotlib.__version__),
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['code/rindler_singlet_packet.py','code/run_rindler_singlet_packet.py']},
        kernel=kernel,matrix_dressing=matrix,radial_refinement=refinement,checks=checks)
    path.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(status=data['status'],checks=len(checks),failed=[c for c in checks if not c['passed']],max_radial_difference=max(r['absolute_difference'] for r in refinement))))
    if data['status']!='passed':raise SystemExit(1)


if __name__=='__main__':main()
