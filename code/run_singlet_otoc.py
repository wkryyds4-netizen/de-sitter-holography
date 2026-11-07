"""Finite-q OTOC and squared commutator, using existing BINT channels."""
from pathlib import Path
from math import exp,log,sqrt,pi
import hashlib,json,platform,time
import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from singlet_otoc import Kernel,cutoff

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-singlet-otoc'


def uncrossed_integral(lam,times,nodes):
    q=exp(-4*lam);r=q*q
    x,w=np.polynomial.legendre.leggauss(nodes);th=(x+1)*pi/2;w*=pi/2
    count=1
    while q**count/(1-q)>1e-17:count+=1
    rho=np.ones(nodes)/(2*pi);den=np.ones((nodes,nodes));num=1.
    diff=th[:,None]-th[None,:];total=th[:,None]+th[None,:]
    for k in range(count):
        z=q**k
        rho*=abs(1-z*np.exp(2j*th))**2*(1-q*z);num*=1-r*r*z
        den*=(1-2*r*z*np.cos(diff)+(r*z)**2)*(1-2*r*z*np.cos(total)+(r*z)**2)
    K=num/den;weight=w*rho;e=2*np.cos(th)/sqrt(lam*(1-q));out=[]
    for t in times:
        inner=K@(weight*np.exp(-1j*t*e))
        terms=weight*np.exp(2j*t*e)*inner**2;value=terms.sum()
        out.append(dict(time=t,value=float(value.real),imag=float(value.imag),
            cancellation_condition=float(sum(abs(terms))/max(abs(value),1e-300))))
    return dict(nodes=nodes,product_terms=count,records=out)


def main():
    OUT.mkdir(exist_ok=True);target=OUT/'otoc-numerics.json'
    if target.exists():raise FileExistsError(target)
    started=time.monotonic();times=np.linspace(0,1.5,31);groups=[]
    fig,axes=plt.subplots(1,2,figsize=(9,3.7))
    for lam in [.25,.1,.05,.02,.01]:
        cut=cutoff(lam,float(times[-1]));kernel=Kernel(lam,cut['L'])
        rows=[kernel.at(float(t)) for t in times]
        fine=Kernel(lam,cut['L']+16);checks=[]
        for j in [0,4,10,20,30]:
            ref=fine.at(float(times[j]));base=rows[j]
            errors={key:abs(base[key]-ref[key]) for key in ['uncrossed','crossed_reduced','otoc','commutator_squared']}
            assert max(errors.values())<1e-10
            checks.append(dict(time=float(times[j]),reference=ref,absolute_errors=errors,
                otoc_relative_error=errors['otoc']/max(abs(ref['otoc']),1e-300)))
        assert max(abs(r['uncrossed_imag']) for r in rows)<1e-10
        assert min(r['commutator_squared'] for r in rows)>-1e-10
        assert abs(rows[0]['otoc']-(2+kernel.Q))<1e-10
        small=kernel.at(.0001)
        coefficient=2/lam*(1-kernel.r)**2*(5+2*kernel.r+3*kernel.r**2)
        residual=abs(small['commutator_squared']/.0001**2/coefficient-1)
        assert residual<1e-5
        groups.append(dict(lambda_value=lam,cutoff=cut,records=rows,L_plus_16_checks=checks,
            small_time=dict(time=.0001,value=small['commutator_squared'],analytic_t2_coefficient=coefficient,relative_residual=residual)))
        axes[0].plot(times,[r['otoc'] for r in rows],label=rf'$\lambda={lam}$')
        axes[1].plot(times,[r['commutator_squared'] for r in rows])
        print('main',lam,'L',cut['L'],'done',round(time.monotonic()-started,2),flush=True)
    axes[0].set_ylabel(r'$F_{\rm OTO}(t)$');axes[1].set_ylabel(r'$\langle [X(t),X]^\dagger[X(t),X]\rangle$')
    for ax in axes:ax.set_xlabel(r'$\mathcal{J}_s t$');ax.grid(alpha=.2)
    axes[0].legend(fontsize=8);axes[1].axhline(2,color='k',ls='--',lw=.7)
    fig.tight_layout();fig.savefig(OUT/'otoc-and-commutator.png',dpi=180);plt.close(fig)
    independent=[]
    for lam in [.25,.1]:
        tt=[0.,.2,.5,1.];coarse=uncrossed_integral(lam,tt,96);fine=uncrossed_integral(lam,tt,192)
        g=next(g for g in groups if g['lambda_value']==lam);kk=Kernel(lam,g['cutoff']['L'])
        ref=[kk.at(t)['uncrossed'] for t in tt]
        errors=[abs(v-r['value']) for v,r in zip(ref,fine['records'])]
        assert max(errors)<1e-10
        independent.append(dict(lambda_value=lam,coarse=coarse,fine=fine,jacobi=ref,absolute_errors=errors))
    mp.mp.dps=60;scaling=[];zs=[.025,.05,.1]
    reference=[float(mp.hyperu(4,1,1/mp.mpf(str(z)))/mp.mpf(str(z))**4) for z in zs]
    for lam in [.02,.01,.005]:
        tt=[log(2*z/lam)/4 for z in zs];cut=cutoff(lam,max(tt));kk=Kernel(lam,cut['L'])
        vals=[kk.at(t) for t in tt];fine=Kernel(lam,cut['L']+16)
        refs=[fine.at(t) for t in tt]
        errors=[abs(v['crossed_reduced']-w['crossed_reduced']) for v,w in zip(vals,refs)]
        assert max(errors)<1e-10
        scaling.append(dict(lambda_value=lam,z=zs,times=tt,cutoff=cut,records=vals,
            asymptotic_reference=reference,crossed_difference_from_asymptotic=[v['crossed_reduced']-r for v,r in zip(vals,reference)],
            crossed_L_plus_16_errors=errors))
        print('scaling',lam,'L',cut['L'],'crossed',[round(v['crossed_reduced'],6) for v in vals],flush=True)
    data=dict(status='passed',units='Js=1; beta=mu=0; DS at fixed lambda before time/small-lambda limits',
        source='BINT1811.02584v2 (3.9),(4.9),(4.10); q=e^-4lambda, r=q^2, Q=q^4 restored',
        groups=groups,independent_uncrossed=independent,scrambling_scaling=scaling,
        error_scope='Analytic bounds cover infinite Jacobi and chord-index truncations only; floating roundoff and quadrature assessed empirically. Scaling discrepancies are not truncation errors or certified asymptotic remainder bounds.',
        environment=dict(python=platform.python_version(),numpy=np.__version__,matplotlib=matplotlib.__version__,mpmath=mp.__version__),
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['code/singlet_otoc.py','code/run_singlet_otoc.py']})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print('saved',round(time.monotonic()-started,2),flush=True)

if __name__=='__main__':main()
