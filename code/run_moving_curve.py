"""Reproducible, separately refined Q transverse/Fourier quadrature."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,hashlib,json,platform,sys,time
import numpy as np
import scipy
from scipy.special import gamma,erfc
from moving_threshold_curve import (ROOT,NU,SIGMA,BETA,f,grid,cj,jq_slice,
    jq_ordered,plateau,fourier)

OUT=ROOT/'results/stage3/moving-curve'
TIMES=[-128,-96,-64,-48,-32,-24,-16,-12,-8,-6,-4,-2,0,2,4,6,8,12,16,24,32,48,64,96,128]
BREAKS=sorted(set([-12.,-10.,10.,12.]+np.arange(-8.,8.25,.5).tolist()))

def hashes(paths):
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}

def run(n,m,radius=10.):
    OUT.mkdir(parents=True,exist_ok=True)
    dest=OUT/f'run-w{n}-chi{m}-R{radius:g}.json'
    if dest.exists():
        raise FileExistsError(f'Preserve saved run: {dest}')
    t=time.perf_counter();chi,w=grid(BREAKS,m);values=[]
    for i,c in enumerate(chi):
        values.append(jq_slice(c,n,radius))
        if (i+1)%200==0:
            print('transverse',n,'Fourier',m,i+1,'/',len(chi),'seconds',round(time.perf_counter()-t,1),flush=True)
    jq=np.array(values);p=plateau(n,radius)
    x,d,scale=fourier(chi,w,jq,p['Xinf'],TIMES)
    x2,d2,scale2=fourier(chi,w,jq,p['Xinf'],TIMES,2.)
    dense=np.linspace(-16,16,257)
    xd,dd,_=fourier(chi,w,jq,p['Xinf'],dense)
    twoB=BETA*f(0,0).real**2/(768*np.pi**7)
    rows=[]
    for i,T in enumerate(TIMES):
        rows.append(dict(T=T,X=cj(x[i]),Xinf_minus_X=cj(d[i]),
            D=float(2*(p['N']-x[i].real)),D_minus_edge=float(2*d[i].real),
            alternate_subtraction_X=cj(x2[i]),subtraction_absolute_difference=float(abs(x[i]-x2[i])),
            early_cubic_ratio=float(2*d[i].real*abs(T)**3/twoB) if T<0 else None,
            late_cubic_ratio=float(2*x[i].real*T**3/twoB) if T>0 else None,
            early_cancellation_scale_over_signal=float(scale/max(abs(d[i].real),1e-300))))
    paths=['code/run_moving_curve.py','code/moving_threshold_curve.py',
        'code/preparation_asymptotics.py','manuscript/preparation.tex','manuscript/logtime.tex',
        'sources/kolchmeyer-liu-2411.08090v1.html']
    out=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),kind='actual limiting Q curve; empirical quadrature',
        parameters=dict(nu=NU,sigma=SIGMA,beta=BETA,transverse_order=n,fourier_order=m,
            transverse_radius=radius,fourier_breaks=BREAKS,subtraction_widths=[1.,2.],
            time_origin='t_M=log(2M)+T; M->infinity first'),
        plateau=p,two_B=twoB,regular_integral_absolute_scale=scale,rows=rows,
        spectrum=dict(chi=chi.tolist(),weights=w.tolist(),JQ=[cj(z) for z in jq]),
        dense_curve=dict(T=dense.tolist(),X=[cj(z) for z in xd],Xinf_minus_X=[cj(z) for z in dd]),
        elapsed_seconds=time.perf_counter()-t,
        versions=dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__),
        hashes=hashes(paths))
    dest.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(dest,flush=True)
    print(json.dumps({'plateau':p,'rows':[v for v in rows if v['T'] in [-128,-64,-32,-16,-8,0,16,64,128]]},indent=2),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--n',type=int,required=True)
    parser.add_argument('--m',type=int,required=True);parser.add_argument('--radius',type=float,default=10.)
    a=parser.parse_args();run(a.n,a.m,a.radius)
