"""Check actual Q curves, keep complex data and small-signal diagnostics."""
from datetime import datetime,timezone
import hashlib,json,sys,platform
import numpy as np
import scipy
from scipy.special import gamma,erfc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from moving_threshold_curve import ROOT,f,jq_slice,jq_ordered,plateau,cj
from run_moving_curve import OUT,hashes

def main():
    names=['run-w24-chi40-R10.json','run-w32-chi40-R10.json','run-w24-chi64-R10.json']
    runs=[json.loads((OUT/n).read_text()) for n in names]
    checks=[]
    def check(name,value,tol):
        checks.append(dict(name=name,value=float(value),tolerance=float(tol),passed=bool(value<tol)))
    geometry=[]
    for c in [-3.,-1.,-.1,0.,.1,1.,3.]:
        direct=jq_slice(c,32,10);ordered=jq_ordered(c,32,10)
        rad8=jq_slice(c,32,8);rad12=jq_slice(c,32,12)
        rel=abs(direct-ordered)/max(abs(direct),1.)
        check('independent-half-plane-chi'+str(c),rel,2e-11)
        check('radius-10-vs-12-chi'+str(c),abs(direct-rad12)/max(abs(direct),1.),2e-12)
        geometry.append(dict(chi=c,slice=cj(direct),ordered=cj(ordered),R8=cj(rad8),R12=cj(rad12),
            relative_representation_difference=rel))
    p=plateau(40,12)
    check('zero-spectrum-vs-independent-plateau',abs(jq_slice(0.,32,10)/(256*np.pi**6)-p['Xinf']),2e-14)
    check('positive-edge-vs-N-X',p['difference_identity_error'],2e-14)
    for key in ['N','Xinf','Dedge']:
        check('plateau-refinement-'+key,abs(runs[0]['plateau'][key]-p[key]),2e-14)
    spectrum0=np.array([complex(*v) for v in runs[0]['spectrum']['JQ']])
    spectrum1=np.array([complex(*v) for v in runs[1]['spectrum']['JQ']])
    spectrum_rel=np.max(abs(spectrum0-spectrum1)/np.maximum(abs(spectrum1),1.))
    check('transverse-complete-spectrum-refinement',spectrum_rel,2e-11)
    errors=[]
    for base,trans,outer in zip(*(r['rows'] for r in runs)):
        T=base['T'];z=[complex(*r['X']) for r in [base,trans,outer]]
        transverse=abs(z[1]-z[0]);fourier=abs(z[2]-z[0])
        check('transverse-X-T'+str(T),transverse,2e-14)
        check('fourier-X-T'+str(T),fourier,2e-14)
        check('subtraction-X-T'+str(T),outer['subtraction_absolute_difference'],2e-14)
        signal=(abs(complex(*outer['Xinf_minus_X']).real) if T<0 else abs(z[2].real))
        relative=max(transverse,fourier)/max(signal,1e-300)
        check('small-signal-relative-T'+str(T),relative,1e-3)
        errors.append(dict(T=T,transverse_absolute_difference=transverse,
            fourier_absolute_difference=fourier,relative_to_real_plateau_difference_or_late_X=relative))
    dense=runs[-1]['dense_curve'];t=np.array(dense['T']);z=np.array([complex(*v) for v in dense['X']])
    check('sampled-Cauchy-Schwarz',max(0.,float(np.max(abs(z))-p['N'])),2e-13)
    check('sampled-positive-commutator',max(0.,float(-np.min(2*(p['N']-z.real)))),2e-13)
    # Bounds for separate tails, evaluated in floating point, not intervals.
    mass=np.exp(np.pi**2/4)*np.sqrt(np.pi);Jbound=gamma(.25)**4*mass
    transverse_tails={}
    for radius in [8.,10.,12.]:
        tail=gamma(.25)**4*mass/2*(erfc(radius+np.pi/2)+erfc(radius-np.pi/2))
        transverse_tails[str(radius)]=dict(J_scalar_tail=float(tail),JQ_box_tail=float(2*Jbound*tail))
    outer_bound=Jbound**2/(512*np.pi**7)*np.sqrt(np.pi/24/np.tanh(6*np.pi))*np.sqrt(2*np.pi)*erfc(12/np.sqrt(2))
    pole_bound=p['Xinf']/(24*np.pi)*np.sqrt(np.pi)*erfc(12)
    fig,ax=plt.subplots(1,2,figsize=(9.4,3.8),layout='constrained')
    ax[0].plot(t,2*(p['N']-z.real),color='#275d8c',label='D_Q(T)')
    ax[0].axhline(p['Dedge'],color='#bd6530',ls='--',label='Static boundary floor')
    ax[0].axhline(2*p['N'],color='#7c7c7c',ls=':',label='Late plateau')
    ax[0].set(xlabel='T = t - log(2M)',ylabel='Commutator norm',title='Moving hard threshold')
    ax[0].legend(frameon=False,fontsize=8)
    for side,key,color in [(-1,'early_cubic_ratio','#275d8c'),(1,'late_cubic_ratio','#bd6530')]:
        selected=sorted([r for r in runs[-1]['rows'] if side*r['T']>=16],key=lambda r:abs(r['T']))
        ax[1].plot([abs(r['T']) for r in selected],[r[key] for r in selected],'-o',ms=3,
            color=color,label='Early: (D - D_edge) |T|^3 / (2B)' if side<0 else 'Late: Re(X) T^3 / B')
    ax[1].axhline(1.,color='#7c7c7c',ls=':');ax[1].set(xlabel='|T|',ylabel='Ratio to leading cubic tail',title='Actual integrals, no time shift',xscale='log')
    ax[1].set_xticks([16,32,64,128],['16','32','64','128'])
    ax[1].xaxis.set_minor_formatter(NullFormatter())
    ax[1].legend(frameon=False,fontsize=7)
    for a in ax:a.grid(alpha=.18)
    fig.savefig(OUT/'moving-curve.pdf');fig.savefig(OUT/'moving-curve.png',dpi=180);plt.close(fig)
    paths=['code/validate_moving_curve.py','code/run_moving_curve.py','code/moving_threshold_curve.py',
        'code/preparation_asymptotics.py']+['results/stage3/moving-curve/'+n for n in names]
    out=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='passed' if all(c['passed'] for c in checks) else 'failed',
        checks=checks,geometry_rows=geometry,independent_plateau=p,curve_refinement_rows=errors,
        complete_spectrum_max_relative_change=float(spectrum_rel),
        tail_expressions_evaluated=dict(transverse=transverse_tails,outer_A_at_12=float(outer_bound),
            outer_subtracted_pole_at_12=float(pole_bound),scope='separate analytic bounds, evaluated without outward rounding; no total PV quadrature error bound'),
        versions=dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__),
        limitations=['empirical quadrature only; no certified total error','complex X retained; D uses real part',
            'ordered M limit; no finite-M or finite-G claim','sampled positivity is not a proof of global monotonicity',
            'same f implementation used by both region representations; this is not independent special-function certification'],
        hashes=hashes(paths))
    (OUT/'diagnostics.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=out['status'],checks=len(checks),failed=[c for c in checks if not c['passed']],
        spectrum_relative_change=float(spectrum_rel),max_transverse_error=max(r['transverse_absolute_difference'] for r in errors),
        max_fourier_error=max(r['fourier_absolute_difference'] for r in errors),
        max_small_signal_relative_change=max(r['relative_to_real_plateau_difference_or_late_X'] for r in errors),tails=out['tail_expressions_evaluated']),indent=2))
    if out['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
