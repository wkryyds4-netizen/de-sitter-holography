"""Physical leading two-branch packets; never labels them exact finite-M data."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib,json,sys,platform
import numpy as np
import scipy
import mpmath
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from logtime_limit import pair_coefficients
from moving_threshold_curve import prefactor,fourier,cj

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/stage3/two-branch'
RUNS=[(24,40),(32,40),(24,64)]
ELLS=[8.,12.,16.,24.,32.,40.]

def load_spectra(n,m):
    dp=f'results/stage3/density-endpoints/run-w{n}-chi{m}.json'
    mp=f'results/stage3/moving-curve/run-w{n}-chi{m}-R10.json'
    d=json.loads((ROOT/dp).read_text());e=json.loads((ROOT/mp).read_text())
    def dec(z):return np.array(z)[:,0]+1j*np.array(z)[:,1]
    data={name:dict(chi=np.array(d['spectrum_nodes']['chi']),w=np.array(d['spectrum_nodes']['weights']),
        jq=dec(model['JQ']),xinf=model['plateau']['Xinf'],edge=model['plateau']['Dedge']) for name,model in d['models'].items()}
    data['gamma0']=dict(chi=np.array(e['spectrum']['chi']),w=np.array(e['spectrum']['weights']),
        jq=dec(e['spectrum']['JQ']),xinf=e['plateau']['Xinf'],edge=e['plateau']['Dedge'])
    return data,[dp,mp]

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    dest=OUT/'diagnostics.json'
    if dest.exists():raise FileExistsError(dest)
    checks=[];identities=[]
    def check(name,error,tol):
        checks.append(dict(name=name,absolute_error=float(error),tolerance=tol,passed=bool(error<tol)))
    c=mpmath.mp.clone();c.dps=60;nu=c.mpf('.6');sg=c.one
    offsets=[['.2','.5','.4','0'],['-.3','.2','.5','0'],['1.2','-.1','.4','0'],
             ['.2','.600001','.4','0'],['.2','.599999','.4','0'],['-.8','-.1','-.8','0']]
    for i,off in enumerate(offsets):
        x,y,z,d=map(c.mpf,off);ch=x-y+z-d
        vp,vm=pair_coefficients(c,off,nu);labs=[x,y,z,d]
        H=c.exp(-c.pi*(y-d))/(1024*c.pi**5)*c.fprod(c.gamma(c.mpf('.25')+c.j*(labs[j]-labs[(j+1)%4]+ep*nu)/2) for j in range(4) for ep in [-1,1])
        Q=x*x+(y-x)**2+(y-z)**2+z*z
        def fm(w,ch):
            return c.exp(-c.pi*w-w*w)*c.fprod(c.gamma(c.mpf('.25')+c.j*ch/4+eta*c.j*(w+ep*nu)/2) for eta in [-1,1] for ep in [-1,1])
        direct=H*vp/c.pi**3*c.exp(-Q/2)
        reduced=2**(-2*c.j*ch)*c.cosh(c.pi*ch/2)*c.exp(-ch*ch/2)*fm(x-ch/2,-ch)*fm(z-ch/2,-ch)/(512*c.pi**7)
        conj=c.conj(H*vm/c.pi**3*c.exp(-Q/2))
        check(f'plus-original-ratio-{i}',abs(direct-reduced)/abs(direct),1e-48)
        check(f'conjugacy-same-domain-{i}',abs(direct-conj)/abs(direct),1e-48)
        identities.append(dict(offsets=off,plus=[str(c.re(direct)),str(c.im(direct))],relative_error=str(abs(direct-reduced)/abs(direct))))
    allruns={};inputs=[]
    for n,m in RUNS:
        data,paths=load_spectra(n,m);inputs+=paths;rows=[]
        for name,v in data.items():
            chi,w,jq,xinf=(v[a] for a in ['chi','w','jq','xinf'])
            amp=prefactor(chi)*jq
            for ell in ELLS:
                s=ell/2;k=2*ell-s;t=ell-s
                xx,dd,scale=fourier(chi,w,jq,xinf,[-s,-k])
                xx2,dd2,_=fourier(chi,w,jq,xinf,[-s,-k],2.)
                # Combined amplitude has no pole, hence needs no PV contact term.
                direct=np.sum(w*(np.exp(-1j*s*chi)*amp+np.exp(1j*k*chi)*amp.conjugate()))
                mirror=xx[0]+xx[1].conjugate()-xinf
                neg=np.sum(w*(np.exp(-1j*k*chi)*amp+np.exp(1j*s*chi)*amp.conjugate()))
                check(f'combined-versus-mirror-{n}-{m}-{name}-{ell}',abs(direct-mirror),3e-12)
                check(f'time-reversal-{n}-{m}-{name}-{ell}',abs(neg-direct.conjugate()),3e-14)
                check(f'PV-subtraction-{n}-{m}-{name}-{ell}',max(abs(dd-dd2)),2e-13)
                slow=2*dd[0].real;fast=2*dd[1].real
                rows.append(dict(model=name,ell=ell,log_M=ell-np.log(2),physical_t=t,T=-s,k=k,
                    leading_cross=cj(mirror),combined_integral=cj(direct),
                    slow_complex_difference=cj(dd[0]),fast_complex_difference=cj(dd[1].conjugate()),
                    slow_real_correction=float(slow),fast_real_correction=float(fast),
                    total_real_correction=float(slow+fast),
                    ratio=float(fast/slow) if abs(slow)>1e-20 else None,
                    cancellation_scale=float(scale),scope='two-branch leading kernel with actual projected packets'))
        allruns[f'w{n}-chi{m}']=rows
        print(n,m,'completed',len(rows),'leading-packet pairs',flush=True)
    base=allruns['w24-chi40'];high=allruns['w24-chi64'];trans=allruns['w32-chi40']
    comparisons=[]
    for a,b,d in zip(base,high,trans):
        for key in ['slow_real_correction','fast_real_correction','total_real_correction']:
            delta=max(abs(b[key]-a[key]),abs(b[key]-d[key]))
            relative=delta/max(abs(b[key]),1e-300)
            # This is an explicit empirical resolution check, not a certified error.
            tol=2e-3 if b['model']!='gap1' else 2e-3
            check(f'actual-small-real-refinement-{b["model"]}-{b["ell"]}-{key}',relative,tol)
            comparisons.append(dict(model=b['model'],ell=b['ell'],quantity=key,absolute_change=delta,relative_change=relative))
    fig,axes=plt.subplots(1,2,figsize=(10.6,3.9),layout='constrained')
    for name,n,color in [('gamma0',3,'#275d8c'),('gamma1',5,'#b35d35'),('gamma2',5,'#498464')]:
        row=[r for r in high if r['model']==name and r['ell']>=24];x=np.array([r['ell'] for r in row]);rat=np.array([r['ratio'] for r in row])
        axes[0].plot(x,rat*3**n,'o-',label=name.replace('gamma','p='),color=color)
    axes[0].axhline(1,color='black',ls='--',lw=1)
    axes[0].set(xlabel='ell = log(2M)',ylabel='(fast / slow) divided by (1/3)^n',title='Actual leading-kernel integrals, s = ell/2')
    axes[0].set_xticks([24,32,40])
    a=np.linspace(.02,.95,150)
    for n,color in [(3,'#275d8c'),(5,'#b35d35')]:
        axes[1].semilogy(a,(a/(2-a))**n,color=color,label=f'power n={n}')
    axes[1].set(xlabel='a = s / ell',ylabel='Asymptotic fast / slow',title='Fixed fraction: relative contribution persists')
    for ax in axes:ax.legend(frameon=False);ax.grid(alpha=.18)
    fig.savefig(OUT/'branch-window.pdf');fig.savefig(OUT/'branch-window.png',dpi=180);plt.close(fig)
    paths=sorted(set(inputs+['code/two_branch_window.py','code/logtime_limit.py','code/crossing_assembly.py',
        'code/moving_threshold_curve.py','code/preparation_asymptotics.py','manuscript/twobranch.tex',
        'sources/kolchmeyer-liu-2411.08090v1.html']))
    result=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='passed' if all(v['passed'] for v in checks) else 'failed',
        checks=checks,identities=identities,runs=allruns,refinement=comparisons,
        parameters=dict(nu=.6,sigma=1.,beta=1.,ells=ELLS,s_over_ell=.5,orders=RUNS,coefficient_dps=60),
        versions=dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__,mpmath=mpmath.__version__,matplotlib=matplotlib.__version__),
        limitations=['leading two-branch kernel, not exact finite-M 6j packet integration','reuses frozen physical spectra; no new transverse integrations',
            'small-real refinement is empirical, not a total quadrature error bound','uniform absolute theorem is an author self-reviewed draft',
            'no explicit growing-box constants, computable relative window, finite-G dynamics or priority claim'],
        hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(result['status'],len(checks),'checks; max small-real change',max(v['relative_change'] for v in comparisons),flush=True)
    for v in checks:
        if not v['passed']:print('FAILED',v,flush=True)
    if result['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
