"""Run and validate the finite-G conformal eikonal baseline."""
from __future__ import annotations
import csv
import hashlib
import json
import platform
import sys
from datetime import datetime,timezone
from pathlib import Path
import mpmath as mp
from conformal_shock import (SOURCE,alpha,shock_profile,angular_average,
    shock_moments,j_boundary,j_complex,j_integral,j_pv_real,disconnected,
    conformal_ratio,imaginary_ratio_integral,perturbative_ratio)

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'/'stage4'
DPS=35


def textnum(x):
    return mp.nstr(x,30)


def save_csv(path,rows):
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)


def main():
    mp.mp.dps=DPS
    OUT.mkdir(parents=True,exist_ok=True)
    checks=[]
    def check(name,ok,error=None,tolerance=None):
        d={'name':name,'passed':bool(ok)}
        if error is not None:d['error']=textnum(error)
        if tolerance is not None:d['tolerance']=str(tolerance)
        checks.append(d)
    m=mp.mpf(1000)
    # Kernel checks are independent of the angular profile implementation.
    for c in [mp.mpc('.2','.1'),mp.mpc('-.2','.05'),mp.mpc('-.7','-.1'),mp.mpc('2','.3')]:
        error=abs(j_complex(c)-j_integral(c))
        check('rational kernel '+str(c),error<mp.mpf('1e-25'),error,'1e-25')
    for c in ['-0.1','-0.5','-2']:
        error=abs(mp.re(j_boundary(c,1))-j_pv_real(c))
        check('subtracted PV kernel '+c,error<mp.mpf('1e-25'),error,'1e-25')
        for side in [-1,1]:
            eps=mp.mpf('1e-10')
            error=abs(j_complex(mp.mpc(c,side*eps))-j_boundary(c,side))
            check('kernel cut side '+c+','+str(side),error<mp.mpf('1e-7'),error,'1e-7')
    check('zero-scattering kernel',j_boundary(0,1)+j_boundary(0,-1)-1==1)
    check('empty-dS normalization',abs(disconnected(m,0)-1/(64*mp.pi**2))<mp.mpf('1e-30'))
    moments=[]
    for deficit_str in ['0','0.001','0.05','0.25','0.6']:
        deficit=mp.mpf(deficit_str)
        h1,h2=shock_moments(m,deficit)
        n1=angular_average(lambda p:shock_profile(p,m,deficit),deficit)
        n2=angular_average(lambda p:shock_profile(p,m,deficit)**2,deficit)
        for n,exact,numerical in [(1,h1,n1),(2,h2,n2)]:
            error=abs(exact-numerical)
            check('angular moment '+str(n)+',deficit='+deficit_str,error<mp.mpf('1e-30'),error,'1e-30')
        moments.append({'deficit':deficit_str,'h1_exact':textnum(h1),'h1_integral':textnum(n1),
                        'h2_exact':textnum(h2),'h2_integral':textnum(n2)})
    data=[]
    for deficit_str in ['0','0.05','0.25']:
        deficit=mp.mpf(deficit_str)
        for u in [-4,-3,-2,-1,0,1,2]:
            t=mp.log(m)+u
            r12=conformal_ratio(m,deficit,t,1)
            r14=conformal_ratio(m,deficit,t,-1)
            imag=imaginary_ratio_integral(m,deficit,t)
            check('conjugate contours deficit='+deficit_str+',u='+str(u),abs(r14-mp.conj(r12))<mp.mpf('1e-27'))
            error=abs(mp.im(r12)-imag)
            check('residue imaginary integral deficit='+deficit_str+',u='+str(u),error<mp.mpf('1e-25'),error,'1e-25')
            f0=disconnected(m,deficit)
            data.append({'m':textnum(m),'deficit':deficit_str,'G':textnum(deficit/(4*m)),
                         'alpha_over_pi':textnum(alpha(deficit)/mp.pi),'u':u,'t':textnum(t),
                         'R12_real':textnum(mp.re(r12)),'R12_imag':textnum(mp.im(r12)),
                         'R14_real':textnum(mp.re(r14)),'R14_imag':textnum(mp.im(r14)),
                         'TOC_ratio_retained_order':'1','F0':textnum(f0),
                         'F12_real':textnum(mp.re(f0*r12)),'F12_imag':textnum(mp.im(f0*r12)),
                         'max_enhancement_c':textnum(mp.exp(t)/(8*m*(mp.sin(mp.pi*deficit)/(mp.pi*deficit) if deficit else 1)))})
        print('completed deficit',deficit_str,flush=True)
    # Regulated computations keep both contours; no real-part averaging.
    regulators=[]
    for u in [0,2]:
        t=mp.log(m)+u
        for contour in [-1,1]:
            boundary=conformal_ratio(m,'0.25',t,contour)
            errors=[]
            for eps in ['1e-2','1e-4','1e-6']:
                value=conformal_ratio(m,'0.25',t,contour,eps)
                error=abs(value-boundary);errors.append(error)
                regulators.append({'u':u,'t':textnum(t),'contour':contour,'epsilon':eps,
                  'real':textnum(mp.re(value)),'imag':textnum(mp.im(value)),
                  'boundary_real':textnum(mp.re(boundary)),'boundary_imag':textnum(mp.im(boundary)),
                  'absolute_error':textnum(error)})
            check('regulated convergence u='+str(u)+',contour='+str(contour),
                  all(mp.mpf('.008')<errors[i+1]/errors[i]<mp.mpf('.012') for i in range(2)))
    # Precision increases and cyclic rotation/time reversal are checked separately.
    precision=[]
    for u in [0,2]:
        with mp.workdps(45):
            value45=conformal_ratio(mp.mpf(1000),mp.mpf('.25'),mp.log(1000)+u,1)
        row=next(r for r in data if r['deficit']=='0.25' and r['u']==u)
        value35=mp.mpc(row['R12_real'],row['R12_imag'])
        error=abs(value45-value35)
        precision.append({'u':u,'error_35_vs_45_dps':textnum(error)})
        check('precision increase u='+str(u),error<mp.mpf('1e-24'),error,'1e-24')
        negative=conformal_ratio(m,'.25',-(mp.log(m)+u),1)
        check('cyclic exchange by time reversal u='+str(u),abs(negative-mp.conj(value35))<mp.mpf('1e-24'))
    # Compare the known coefficients only in the small-enhancement region.
    for t_str in ['0','1','2']:
        t=mp.mpf(t_str)
        value=conformal_ratio(m,'.25',t,1)
        approx=perturbative_ratio(m,'.25',t)
        hmax=1/(2*m*(mp.sin(mp.pi/4)/(mp.pi/4)))
        third_scale=mp.mpf(6)*hmax**3*(mp.exp(3*t)+mp.exp(-3*t))/64
        error=abs(value-approx)
        check('small-enhancement second order t='+t_str,error<2*third_scale,error,2*third_scale)
    # A genuinely continuous G->0 test at fixed m and fixed t.
    t=mp.log(m)
    pure=conformal_ratio(m,0,t,1)
    recoil_limit=[]
    for deficit_str in ['1e-2','1e-3','1e-4']:
        value=conformal_ratio(m,deficit_str,t,1)
        error=abs(value-pure)
        recoil_limit.append({'deficit':deficit_str,'absolute_error_to_G0':textnum(error)})
    errs=[mp.mpf(r['absolute_error_to_G0']) for r in recoil_limit]
    check('fixed-m pure-recoil G limit',all(mp.mpf('.08')<errs[i+1]/errs[i]<mp.mpf('.12') for i in range(2)))
    for phi in ['0','0.7','1.8','3']:
        error=abs(shock_profile(phi,m,0)-mp.cos(mp.mpf(phi))/(2*m))
        check('G0 profile phi='+phi,error<mp.mpf('1e-30'),error,'1e-30')
    # Fixed deficit, increasing mass at fixed t must remove scattering.
    deviations=[abs(conformal_ratio(mass,'.25',1,1)-1) for mass in [1000,2000,4000]]
    check('fixed-time no-scattering mass limit',all(deviations[i+1]<deviations[i] for i in range(2)))
    save_csv(OUT/'correlators.csv',data)
    save_csv(OUT/'regulators.csv',regulators)
    save_csv(OUT/'shock_moments.csv',moments)
    hashes={name:hashlib.sha256((ROOT/'code'/name).read_bytes()).hexdigest()
            for name in ['conformal_shock.py','run_conformal_shock.py']}
    study={'generated_utc':datetime.now(timezone.utc).isoformat(),'source':SOURCE,
        'source_equations':['G.2','2.12','2.25','2.26','2.33','2.35','2.39','2.40','2.42'],
        'scope':'Two conformal probes, finite-G eikonal approximation; no finite-clock state or selected physical contour.',
        'normalization':'F0=1/(16 alpha^2); disconnected/TOC baseline only at retained eikonal order.',
        'precision_dps':DPS,'nu':'i/2','nu_prime':'i/2','m':'1000',
        'correlators':data,'regulators':regulators,'moments':moments,
        'precision_checks':precision,'pure_recoil_limit':recoil_limit,
        'code_sha256':hashes,'environment':{'python':sys.version,'mpmath':mp.__version__,'platform':platform.platform()}}
    failed=sum(not c['passed'] for c in checks)
    validation={'status':'failed' if failed else 'passed','checks':len(checks),'failed':failed,
                'check_details':checks,'code_sha256':hashes,'scope':study['scope']}
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    study['environment']['matplotlib']=matplotlib.__version__
    fig,axs=plt.subplots(1,2,figsize=(10.4,4.2),constrained_layout=True)
    for deficit_str in ['0','0.05','0.25']:
        rows=[r for r in data if r['deficit']==deficit_str]
        xs=[r['u'] for r in rows]
        line=axs[0].plot(xs,[float(r['R12_real']) for r in rows],'o-',label='4Gm='+deficit_str)[0]
        color=line.get_color()
        axs[1].plot(xs,[float(r['R12_imag']) for r in rows],'o-',color=color,label='12: 4Gm='+deficit_str)
        axs[1].plot(xs,[float(r['R14_imag']) for r in rows],'--',color=color,label='14: 4Gm='+deficit_str)
    axs[0].axhline(1,color='gray',ls=':',label='disconnected / TOC baseline')
    axs[0].set(xlabel=r'$u=t-\log m$',ylabel=r'$\mathrm{Re}(\mathcal{F}/\mathcal{F}_0)$',title='Two conformal probes, m=1000')
    axs[1].set(xlabel=r'$u=t-\log m$',ylabel=r'$\mathrm{Im}(\mathcal{F}/\mathcal{F}_0)$',title='Both lateral contours retained')
    for ax in axs:ax.grid(alpha=.25);ax.legend(fontsize=7)
    fig.savefig(OUT/'conformal-shock.png',dpi=180);plt.close(fig)
    (OUT/'study.json').write_text(json.dumps(study,indent=2)+'\n',encoding='utf-8')
    (OUT/'validation.json').write_text(json.dumps(validation,indent=2)+'\n',encoding='utf-8')
    print(f'{len(checks)-failed} passed; {failed} failed',flush=True)
    for c in checks:
        if not c['passed']:print(c)
    if failed:raise SystemExit(1)


if __name__=='__main__':main()
