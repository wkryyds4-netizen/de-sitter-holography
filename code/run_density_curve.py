"""Full integrals for endpoint-vanishing and gapped packets, not fitted tails."""
from datetime import datetime,timezone
import argparse,json,hashlib,platform,sys,time
import numpy as np
import scipy
from density_endpoints import ROOT,MODELS,jq_family,plateaux,tail,coefficients,f,prefactor,fourier,cj,grid

OUT=ROOT/'results/stage3/density-endpoints'
TIMES=[-128,-96,-64,-48,-32,-24,-16,-12,-8,-6,-4,-2,0,2,4,6,8,12,16,24,32,48,64,96,128]
BREAKS=sorted(set([-12.,-10.,10.,12.]+np.arange(-8.,8.25,.5).tolist()))
def hashes(paths):return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}

def run(n,m):
    OUT.mkdir(parents=True,exist_ok=True);dest=OUT/f'run-w{n}-chi{m}.json'
    if dest.exists():raise FileExistsError(dest)
    start=time.perf_counter();chi,w=grid(BREAKS,m);spectra={v['name']:[] for v in MODELS}
    for i,c in enumerate(chi):
        values=jq_family(c,n)
        for name,z in values.items():spectra[name].append(z)
        if (i+1)%200==0:print(n,m,i+1,'/',len(chi),'seconds',round(time.perf_counter()-start,1),flush=True)
    dense=np.linspace(-16,16,257);models={}
    for model in MODELS:
        name=model['name'];p=plateaux(model,n);jq=np.array(spectra[name])
        x,d,scale=fourier(chi,w,jq,p['Xinf'],TIMES)
        x2,d2,_=fourier(chi,w,jq,p['Xinf'],TIMES,2.)
        xd,dd,_=fourier(chi,w,jq,p['Xinf'],dense)
        rows=[]
        for i,T in enumerate(TIMES):
            correction=-d[i] if T<0 else x[i]
            predicted=tail(model,T) if T else 0j
            if model['gap']>0:
                C=model['beta']*prefactor(-model['gap'])*f(-model['gap']/2,-model['gap'])**2
                demod=correction*T**4*np.exp(1j*model['gap']*T)/C if T else 0j
                leading=None
            else:
                cc=coefficients(model['p'],model['beta']);porder=model['p']
                leading=(2*(-1)**(porder//2)*cc['S']/(porder+3) if porder%2==0
                         else 2*(-1)**((porder+1)//2)*cc['S']*cc['d'])
                power=porder+3 if porder%2==0 else porder+4
                leading_ratio=2*d[i].real*abs(T)**power/leading if T<0 else None
                demod=None
            rows.append(dict(T=T,X=cj(x[i]),Xinf_minus_X=cj(d[i]),D=float(2*(p['N']-x[i].real)),
                D_minus_edge=float(2*d[i].real),subtraction_absolute_difference=float(abs(x[i]-x2[i])),
                predicted_complex_correction=cj(predicted),actual_complex_correction=cj(correction),
                early_real_tail_ratio=leading_ratio if leading is not None and T<0 else None,
                gap_demodulated_ratio=cj(demod) if demod is not None else None,
                cancellation_scale_over_real_correction=float(scale/max(abs(correction.real),1e-300))))
        models[name]=dict(model=model,plateau=p,rows=rows,JQ=[cj(z) for z in jq],
            dense_curve=dict(T=dense.tolist(),X=[cj(z) for z in xd],Xinf_minus_X=[cj(z) for z in dd]),
            regular_integrand_absolute_scale=scale)
    paths=['code/run_density_curve.py','code/density_endpoints.py','code/moving_threshold_curve.py',
        'code/preparation_asymptotics.py','manuscript/preparation.tex','manuscript/movingcurve.tex']
    out=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),kind='actual hard-threshold limiting integrals; empirical quadrature',
        parameters=dict(nu=.6,sigma=1.,transverse_order=n,fourier_order=m,radius=10.,chi_breaks=BREAKS,
            time_origin='t_M=log(2M)+T; M first',models=MODELS),
        spectrum_nodes=dict(chi=chi.tolist(),weights=w.tolist()),models=models,elapsed_seconds=time.perf_counter()-start,
        versions=dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__),hashes=hashes(paths))
    dest.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(dest,flush=True)
    for name,model in models.items():
        print(name,model['plateau'],flush=True)
        for row in model['rows']:
            if row['T'] in [-128,-64,-32,-16,-8]:
                print(row['T'],'D-edge',row['D_minus_edge'],'ratio',row['early_real_tail_ratio'],'demod',row['gap_demodulated_ratio'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,required=True);p.add_argument('--m',type=int,required=True)
    a=p.parse_args();run(a.n,a.m)
