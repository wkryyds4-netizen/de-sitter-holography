"""Parallel refinement of the conditional exact crossed four-point task."""
from pathlib import Path
from datetime import datetime,timezone
from concurrent.futures import ProcessPoolExecutor
import argparse,json,hashlib,time,platform,os
import numpy as np,scipy
from run_crossing_conditional import rule,positive_norms,cj,sha
from quantized_observer import mass
from observer_crossing import CrossingEvaluator
ROOT=Path(__file__).resolve().parents[1]
_EV=None
def init(dps):
    global _EV
    _EV=CrossingEvaluator(dps,8)
def point(task):
    a,b,c,w,sigma=task;d=.8;nu=.6
    labels=tuple(format(float(x),".17g") for x in (a,b,c,d))
    r=_EV.evaluate(str(nu),*labels)
    raw=complex(r.hh_cross/_EV.rho_hh(str(d)))
    ma,mb,mc,md=map(float,map(mass,(a,b,c,d)))
    damp=np.exp(-.5*sigma*sigma*((md-mc)**2+(mc-mb)**2+(mb-ma)**2+(ma-md)**2))
    z=raw*damp*w
    values=[z*np.exp(-1j*t*(md-mc+mb-ma)) for t in (0.,.5,1.,2.)]
    return values,abs(z),float(r.cancellation_ratio),abs(raw.imag),labels,raw
def run(n,upper,sigma,workers,dps):
    aa,wa=rule(n,.2,upper);bb,wb=rule(n,.2,upper)
    tasks=[]
    for a,aw in zip(aa,wa):
        cc,wc=rule(n,a,upper)
        for b,bw in zip(bb,wb):
            for c,cw in zip(cc,wc):tasks.append((a,b,c,2*aw*bw*cw,sigma))
    start=time.perf_counter();vals=[];absolute=0.;cancel=0.;im=0.;worst=None
    with ProcessPoolExecutor(max_workers=workers,initializer=init,initargs=(dps,)) as pool:
        for j,(v,av,cr,imag,labels,raw) in enumerate(pool.map(point,tasks,chunksize=4),1):
            vals.append(v);absolute+=av;im=max(im,imag)
            if cr>cancel:cancel=cr;worst=dict(labels=labels,conditional=cj(raw),cancellation=cr)
            if j%128==0 or j==len(tasks):print(json.dumps(dict(order=n,points=j,total=len(tasks),seconds=time.perf_counter()-start)),flush=True)
    return dict(order=n,upper=upper,lower=.2,nu=.6,d=.8,sigma_proper=sigma,dps=dps,guard_digits=8,
        points=len(tasks),workers=workers,values=[cj(z) for z in np.sum(vals,axis=0)],
        times_left_to_right=[[t,0,t,0] for t in (0.,.5,1.,2.)],
        quadrature_absolute_sum=absolute,maximum_cancellation=cancel,maximum_raw_imaginary=im,
        worst_point=worst,elapsed_seconds=time.perf_counter()-start,
        positive_norms=positive_norms(.6,.8,.2,upper,sigma))
def main():
    p=argparse.ArgumentParser();p.add_argument("--orders",type=int,nargs="+",default=[10,12])
    p.add_argument("--upper",type=float,default=3.);p.add_argument("--sigma",type=float,default=1.)
    p.add_argument("--workers",type=int,default=4);p.add_argument("--dps",type=int,default=15)
    p.add_argument("--label",default="conditional-crossing-refined")
    a=p.parse_args();out=ROOT/"results/stage3/fourpoint"/(a.label+".json")
    if out.exists():raise FileExistsError("Use a new label.")
    s=dict(generated_utc=datetime.now(timezone.utc).isoformat(),status="running",
        scope="Conditional fixed d=.8 exact finite-mass crossed task, finite cutoff, no normalized external packet.",
        numerical_limits="Gauss refinement and positive-tail estimates, no interval certification.",
        source="https://arxiv.org/html/2411.08090v1",
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        code_sha256={n:sha(ROOT/"code"/n) for n in ["run_crossing_refined.py","run_crossing_conditional.py","observer_crossing.py","quantized_observer.py"]},
        input_sha256={"sources/kolchmeyer-liu-2411.08090v1.html":sha(ROOT/"sources/kolchmeyer-liu-2411.08090v1.html")},
        rows=[])
    def save():out.write_text(json.dumps(s,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    save()
    try:
        for n in a.orders:
            row=run(n,a.upper,a.sigma,a.workers,a.dps);s["rows"].append(row);save()
            print(json.dumps(dict(completed_order=n,values=row["values"])),flush=True)
        s["status"]="completed";save()
    except BaseException as exc:
        s["status"]="failed";s["exception"]=repr(exc);save();raise
if __name__=="__main__":main()
