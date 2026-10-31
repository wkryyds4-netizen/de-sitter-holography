"""Endpoint coefficients and actual full-curve diagnostics; not certification."""
from datetime import datetime,timezone
import json,math,platform,sys
import numpy as np
import scipy,mpmath
from scipy.special import digamma
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from density_endpoints import (ROOT,MODELS,jq_family,jq_halfplane,plateaux,
    coefficients,corner_scaled,prefactor,f,fourier,cj)
from run_density_curve import OUT,hashes

def main():
    names=['run-w24-chi40.json','run-w32-chi40.json','run-w24-chi64.json']
    runs=[json.loads((OUT/n).read_text()) for n in names];checks=[]
    def check(name,value,tol):checks.append(dict(name=name,value=float(value),tolerance=float(tol),passed=bool(value<tol)))
    c=mpmath.mp.clone();c.dps=50;nu=c.mpf('.6')
    F=abs(c.gamma(c.mpf('.25')+c.j*nu/2))**4
    H=2*c.re(c.digamma(c.mpf('.25')+c.j*nu/2))+c.euler+2*c.log(2)
    check('mp-F',abs(coefficients(1,1)['F']-float(F)),2e-12)
    check('mp-H',abs(coefficients(1,1)['H']-float(H)),2e-13)
    geometry=[]
    for chi in [-1.1,-1.,-.9,-.1,0.,.1,1.]:
        vals=jq_family(chi,32)
        rad=jq_family(chi,32,12.)
        for model in MODELS:
            name=model['name'];a=vals[name];b=jq_halfplane(chi,model,32)
            rel=abs(a-b)/max(abs(a),1.)
            check(name+'-region-'+str(chi),rel,2e-11)
            check(name+'-R12-'+str(chi),abs(a-rad[name])/max(abs(a),1.),2e-12)
            geometry.append(dict(model=name,chi=chi,slice=cj(a),halfplane=cj(b),R12=cj(rad[name]),relative_difference=rel))
    corner=[]
    for p in range(4):
        base=float(F)**2/math.factorial(p+3)
        expected=coefficients(p,1.)['d']-2j*digamma(.25+.3j).real
        for h in [.02,.01,.005,.001,.00025]:
            low=corner_scaled(h,p,1.,n=16);high=corner_scaled(h,p,1.,n=24)
            ratio=high/base;slope=(ratio-1)/h
            check(f'corner-rule-p{p}-h{h}',abs(low-high)/base,2e-12)
            if h==.00025:check(f'corner-first-derivative-p{p}',abs(slope-expected),.002)
            corner.append(dict(p=p,h=h,scaled_integral=cj(high),ratio=cj(ratio),difference_quotient=cj(slope),expected_slope=cj(expected)))
    gapbase=f(-.5,-1)**2/6
    gapcorner=[]
    for h in [.02,.01,.005,.001,.00025]:
        a=corner_scaled(h,0,1.,gap=1.,n=16);b=corner_scaled(h,0,1.,gap=1.,n=24)
        check('gap-corner-rule-'+str(h),abs(a-b)/abs(gapbase),2e-12)
        gapcorner.append(dict(h=h,ratio=cj(b/gapbase)))
    check('gap-corner-coefficient',abs(complex(*gapcorner[-1]['ratio'])-1),.003)
    refinements=[];summary={}
    for model in MODELS:
        name=model['name'];p=plateaux(model,40,12.)
        check(name+'-independent-zero',abs(jq_family(0.,32)[name]/(256*np.pi**6)-p['Xinf']),3e-14)
        check(name+'-positive-edge',p['identity_error'],3e-14)
        for key in ['N','Xinf','Dedge']:
            check(name+'-plateau-'+key,abs(runs[0]['models'][name]['plateau'][key]-p[key]),3e-14)
        errors=[]
        for rows in zip(*(r['models'][name]['rows'] for r in runs)):
            T=rows[0]['T']
            z=[complex(*r['X']) for r in rows]
            # Compare the directly computed small differences, not X after
            # its plateau-sized addition has rounded away low bits.
            corr=[complex(*r['actual_complex_correction']) for r in rows]
            absx=[abs(z[1]-z[0]),abs(z[2]-z[0])]
            real_change=max(abs(corr[1].real-corr[0].real),abs(corr[2].real-corr[0].real))
            rel=real_change/max(abs(corr[2].real),1e-300)
            for i,kind in enumerate(['transverse','fourier']):check(f'{name}-{kind}-T{T}',absx[i],3e-14)
            check(f'{name}-subtraction-T{T}',rows[-1]['subtraction_absolute_difference'],3e-14)
            if abs(T)<=64:check(f'{name}-small-real-T{T}',rel,.002)
            errors.append(dict(model=name,T=T,transverse_absolute_difference=absx[0],fourier_absolute_difference=absx[1],
                direct_small_real_absolute_change=real_change,relative_small_real_change=rel,
                within_point_two_percent_empirical_change=bool(rel<.002)))
        refinements.extend(errors)
        zz=np.array([complex(*v) for v in runs[-1]['models'][name]['dense_curve']['X']])
        check(name+'-sampled-CS',max(0.,float(np.max(abs(zz))-p['N'])),3e-13)
        check(name+'-sampled-D-positive',max(0.,float(-np.min(2*(p['N']-zz.real)))),3e-13)
        summary[name]=dict(plateau=p,max_small_real_relative_change_to64=max(v['relative_small_real_change'] for v in errors if abs(v['T'])<=64),
            max_small_real_relative_change_all=max(v['relative_small_real_change'] for v in errors))
    # For the plotted gapped oscillations, recompute only the cheap Fourier
    # sums from each frozen spectrum. This is still an actual integral.
    s=np.linspace(64,128,513);gapcurves=[];C=prefactor(-1)*f(-.5,-1)**2
    for run in runs:
        nodes=run['spectrum_nodes'];mm=run['models']['gap1'];jq=np.array([complex(*v) for v in mm['JQ']])
        x,d,scale=fourier(np.array(nodes['chi']),np.array(nodes['weights']),jq,mm['plateau']['Xinf'],-s)
        gapcurves.append(2*d.real*s**4/(2*abs(C)))
    gap_refinement=max(np.max(abs(gapcurves[1]-gapcurves[0])),np.max(abs(gapcurves[2]-gapcurves[0])))
    check('gap-oscillation-amplitude-refinement',gap_refinement,1e-3)
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for name,color,label in [('gamma1','#275d8c','Linear endpoint: p = 1'),('gamma2','#bd6530','Quadratic endpoint: p = 2')]:
        rows=sorted([r for r in runs[-1]['models'][name]['rows'] if -64<=r['T']<=-16],key=lambda r:abs(r['T']))
        axes[0].plot([abs(r['T']) for r in rows],[r['early_real_tail_ratio'] for r in rows],'-o',ms=3,color=color,label=label)
    axes[0].axhline(1,color='gray',ls=':');axes[0].set(xlabel='|T|',ylabel='Signed real correction / leading term',title='Both leading real corrections are negative',xscale='log')
    axes[0].set_xticks([16,24,32,48,64],['16','24','32','48','64']);axes[0].xaxis.set_minor_formatter(NullFormatter())
    axes[0].legend(frameon=False,fontsize=8)
    axes[1].plot(s,gapcurves[-1],color='#275d8c',label='Actual full integral',lw=1.6)
    axes[1].plot(s,-np.cos(np.angle(C)+s),color='#bd6530',ls='--',lw=1,label='Leading endpoint term')
    axes[1].set(xlabel='|T|',ylabel='(D - D_edge) |T|^4 / (2 |C_g|)',title='Support gap g = 1: oscillating correction')
    axes[1].legend(frameon=False,fontsize=8)
    for a in axes:a.grid(alpha=.18)
    fig.savefig(OUT/'endpoint-curves.pdf');fig.savefig(OUT/'endpoint-curves.png',dpi=180);plt.close(fig)
    paths=['code/validate_density_endpoints.py','code/run_density_curve.py','code/density_endpoints.py',
        'code/moving_threshold_curve.py','code/preparation_asymptotics.py']+['results/stage3/density-endpoints/'+n for n in names]
    output=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='passed' if all(v['passed'] for v in checks) else 'failed',
        checks=checks,geometry_rows=geometry,corner_rows=corner,gap_corner_rows=gapcorner,
        references_50dps=dict(F=str(F),H=str(H)),summary=summary,refinement_rows=refinements,gap_C=cj(C),
        gap_plot=dict(s=s.tolist(),scaled_actual=gapcurves[-1].tolist(),scaled_leading=(-np.cos(np.angle(C)+s)).tolist(),
            maximum_rule_change=float(gap_refinement)),
        limits=['author self-reviewed proof drafts, not independently certified','total quadrature errors not certified',
            'Gamma p=1,2 at |T|=96,128 retain diminished small-real precision; primary quantitative comparison uses |T|<=64',
            'gap residual after leading coefficient still visible at |T|=128; no fitted time shift',
            'no finite-M window, finite-G result, or full holographic dual'],
        versions=dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__,mpmath=mpmath.__version__,matplotlib=matplotlib.__version__),
        hashes=hashes(paths))
    (OUT/'diagnostics.json').write_text(json.dumps(output,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=output['status'],checks=len(checks),failed=[v for v in checks if not v['passed']],summary=summary,gap_plot_refinement=float(gap_refinement)),indent=2))
    if output['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
