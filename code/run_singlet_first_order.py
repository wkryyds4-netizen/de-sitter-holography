"""Finite-h comparisons to the sourced leading connected four-time function."""
from pathlib import Path
import hashlib,json,platform,time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from singlet_otoc import cutoff
from singlet_general_fourpoint import GeneralKernel
from singlet_first_order import connected_coefficient,alternating_commutator,adjacent_cumulant
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-first-order'


def main():
    path=OUT/'finite-lambda.json'
    if path.exists():raise FileExistsError(path)
    old=json.loads((ROOT/'results/susskind-singlet-otoc/otoc-numerics.json').read_text())
    tt=[.1,.2,.35,.5];configs=[[.3,.1,-.1,0],[.2,-.05,.35,.1],[.15,.15,-.15,-.15]]
    groups=[];validations=[];start=time.monotonic()
    for lam in [.04,.02,.01,.005]:
        cut=cutoff(lam,.6);k=GeneralKernel(lam,cut['L']);h=4*lam
        og=next((v for v in old['groups'] if v['lambda_value']==lam),None)
        alternating=[]
        for t in tt:
            if og:
                row=min(og['records'],key=lambda r:abs(r['time']-t));assert abs(row['time']-t)<1e-12
            else:row=k.at(t)
            cc=float(alternating_commutator(2*t));kc=float(adjacent_cumulant(2*t))
            actualc=row['commutator_squared'];actualk=row['f22']-1-2*row['g2']**2
            alternating.append(dict(time=t,raw=row,reused_previous=bool(og),commutator_over_h=actualc/h,
                commutator_first=cc,commutator_coefficient_error=actualc/h-cc,
                adjacent_over_h=actualk/h,adjacent_first=kc,adjacent_coefficient_error=actualk/h-kc))
        general=[]
        for times in configs:
            val=k.four(times);leading=connected_coefficient(2*np.array(times));pred=[leading['direct'],leading['nested'],leading['crossed']]
            actual=(np.array(val['channels'])-np.array(val['pair_products']))/h
            val.update(first_order_channels=pred,connected_channels_over_h=actual.tolist(),
                channel_coefficient_errors=(actual-pred).tolist(),first_order_total=leading['total'],connected_over_h=val['connected']/h)
            general.append(val)
        # A selected L+16 refinement for every lambda, including nonalternating times.
        fine=GeneralKernel(lam,cut['L']+16);ref=fine.four(configs[0]);base=general[0]
        errors=np.abs(np.array(ref['channels'])-np.array(base['channels']))
        assert max(errors)<1e-10
        validations.append(dict(kind='L+16',lambda_value=lam,channel_absolute_errors=errors.tolist(),
            reference=ref,connected_absolute_difference=abs(ref['connected']-base['connected'])))
        if lam==.04:
            g=k.four([.2,0,.2,0]);r=k.at(.2)
            assert abs(g['total']-r['otoc'])<1e-11
            validations.append(dict(kind='alternating reduction',absolute_difference=abs(g['total']-r['otoc'])))
            ts=np.array(configs[0]);base=general[0]
            for name,newts in [('cyclic',np.roll(ts,1)),('reversal',ts[::-1]),('shift',ts+.23)]:
                new=k.four(newts);err=abs(new['total']-base['total']);assert err<1e-11
                validations.append(dict(kind=name,absolute_difference=err))
        assert max(abs(v) for r in general for v in r['channel_imag'])<1e-10
        groups.append(dict(lambda_value=lam,h=h,cutoff=cut,alternating=alternating,general=general))
        print('finished',lam,'L',cut['L'],'seconds',round(time.monotonic()-start,1),flush=True)
    extrap=[]
    for i,t in enumerate(tt):
        a=groups[-1]['alternating'][i];b=groups[-2]['alternating'][i]
        extrap.append(dict(time=t,commutator_richardson=2*a['commutator_over_h']-b['commutator_over_h'],
            commutator_first=a['commutator_first'],adjacent_richardson=2*a['adjacent_over_h']-b['adjacent_over_h'],adjacent_first=a['adjacent_first']))
    fig,axs=plt.subplots(1,2,figsize=(9,3.7));fine_t=np.linspace(0,.55,150)
    axs[0].plot(fine_t,alternating_commutator(2*fine_t),'k-',label='leading coefficient')
    axs[1].plot(fine_t,adjacent_cumulant(2*fine_t),'k-')
    for group in groups:
        lam=group['lambda_value'];rows=group['alternating']
        axs[0].plot(tt,[r['commutator_over_h'] for r in rows],'o--',ms=3,label=fr'$\lambda={lam}$')
        axs[1].plot(tt,[r['adjacent_over_h'] for r in rows],'o--',ms=3)
    axs[0].set(xlabel=r'$\mathcal{J}_s t$',ylabel=r'$\mathcal{C}(t)/(4\lambda)$')
    axs[1].set(xlabel=r'$\mathcal{J}_s t$',ylabel=r'$\kappa_{22}(t)/(4\lambda)$')
    axs[0].legend(fontsize=8);axs[1].axhline(0,color='gray',lw=.6)
    fig.tight_layout();fig.savefig(OUT/'fixed-time-first-order.png',dpi=180);plt.close(fig)
    inputs=['code/singlet_first_order.py','code/singlet_general_fourpoint.py','code/run_singlet_first_order.py','code/singlet_otoc.py','results/susskind-singlet-otoc/otoc-numerics.json']
    data=dict(status='passed',status_scope='Algorithmic checks passed; asymptotic discrepancies are measured, not asserted zero or rigorously bounded',
        units='Js=1, h=4 lambda, u_j=2t_j, beta=mu=0, fixed-time small lambda after DS',
        groups=groups,validations=validations,richardson=extrap,
        environment=dict(python=platform.python_version(),numpy=np.__version__,matplotlib=matplotlib.__version__),
        error_scope='Previous cutoff parameters reused conservatively. L+16 and symmetry checks assess truncation and floating arithmetic empirically; the general four-time total bound is not certified here. No uniform h remainder certification.',
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    path.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(status='passed',groups=len(groups),richardson=extrap,seconds=time.monotonic()-start)))


if __name__=='__main__':main()
