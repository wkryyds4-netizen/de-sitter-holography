"""Conditional crossed dS2 task; finite spectral cutoff and exact 6j kernel.

nu=.6, external spectral label d=.8, per-insertion Q at s=.2.
Each Gaussian has proper-time width sigma. Times are left-to-right (t,0,t,0).
A/C reflection reduces the zero/OTO task to a<c; no singular nodes sampled.
Fixed d is a conditional density, not a normalizable delta-mass state.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,json,hashlib,time,platform
import numpy as np
import scipy
from scipy.special import roots_legendre
from scipy.integrate import quad
from quantized_observer import mass,right_kernel
from observer_crossing import CrossingEvaluator

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/stage3/fourpoint"
def cj(z):return dict(real=float(z.real),imag=float(z.imag))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rule(n,lo,hi):
    x,w=roots_legendre(n)
    return lo+(hi-lo)*(x+1)/2,w*(hi-lo)/2
def positive_norms(nu,d,lo,upper,sigma,cap=12.):
    def k(a,b):return float(right_kernel(nu,a,b))*np.exp(-(sigma*(mass(a)-mass(b)))**2)
    def h(a,low,high):return quad(lambda b:k(b,a),low,high,epsabs=2e-13,epsrel=2e-10,limit=150)[0]
    ncut=quad(lambda a:k(a,d)*h(a,lo,upper),lo,upper,epsabs=2e-13,epsrel=2e-10,limit=150)[0]
    # Direct positive tails avoid subtraction of nearly equal norms.
    outer_tail=quad(lambda a:k(a,d)*h(a,lo,cap),upper,cap,epsabs=2e-14,epsrel=2e-9,limit=150)[0]
    inner_tail=quad(lambda a:k(a,d)*h(a,upper,cap),lo,upper,epsabs=2e-14,epsrel=2e-9,limit=150)[0]
    eps=outer_tail+inner_tail
    n=ncut+eps
    bound=np.sqrt(eps)*(np.sqrt(n)+np.sqrt(ncut))
    return dict(N_cut=ncut,N_to_cap=n,epsilon_to_cap=eps,outer_tail=outer_tail,
                inner_tail=inner_tail,cap=cap,cross_tail_bound_estimate=bound,
                scope="Positive quadrature estimate to finite cap; not a certified infinite-tail bound.")
def evaluate_order(order,upper,sigma,dps=15,guard=8):
    nu=.6;d=.8;lo=.2
    ev=CrossingEvaluator(dps,guard)
    ctx=ev.ctx
    aa,wa=rule(order,lo,upper);bb,wb=rule(order,lo,upper)
    times=(0.,.5,1.,2.)
    sums=[0j]*len(times);absolute=0.;max_cancel=0.;max_imag=0.;max_seconds=0.
    samples=[];start=time.perf_counter();done=0
    for ai,(a,aw) in enumerate(zip(aa,wa)):
        cc,wc=rule(order,a,upper)
        for b,bw in zip(bb,wb):
            for c,cw in zip(cc,wc):
                labels=tuple(format(float(x),".17g") for x in (a,b,c,d))
                r=ev.evaluate(str(nu),*labels)
                raw=complex(r.hh_cross/ev.rho_hh(str(d)))
                ma,mb,mc,md=map(float,map(mass,(a,b,c,d)))
                damp=np.exp(-.5*sigma*sigma*((md-mc)**2+(mc-mb)**2+(mb-ma)**2+(ma-md)**2))
                weighted=raw*damp*(2*aw*bw*cw)
                for j,t in enumerate(times):
                    sums[j]+=weighted*np.exp(-1j*t*(md-mc+mb-ma))
                absolute+=abs(weighted)
                max_cancel=max(max_cancel,float(r.cancellation_ratio))
                max_imag=max(max_imag,abs(raw.imag))
                max_seconds=max(max_seconds,r.elapsed_seconds)
                if len(samples)<4 or float(r.cancellation_ratio)>max(s["cancellation"] for s in samples):
                    samples.append(dict(labels=labels,conditional=cj(raw),cancellation=float(r.cancellation_ratio)))
                    samples=sorted(samples,key=lambda s:s["cancellation"],reverse=True)[:6]
                done+=1
        print(json.dumps(dict(order=order,upper=upper,outer_node=ai+1,points=done,
                             elapsed_seconds=time.perf_counter()-start)),flush=True)
    return dict(order=order,upper=upper,lower=lo,nu=nu,d=d,sigma_proper=sigma,
                dps=dps,guard_digits=guard,points=done,
                times_left_to_right=[[t,0,t,0] for t in times],
                values=[cj(z) for z in sums],quadrature_absolute_sum=absolute,
                maximum_cancellation=max_cancel,maximum_raw_imaginary=max_imag,
                slowest_point_seconds=max_seconds,diagnostic_samples=samples,
                elapsed_seconds=time.perf_counter()-start,
                positive_norms=positive_norms(nu,d,lo,upper,sigma))
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--orders",type=int,nargs="+",default=[4,6,8])
    p.add_argument("--upper",type=float,default=3.)
    p.add_argument("--sigma",type=float,default=1.)
    p.add_argument("--label",default="conditional-crossing-pilot")
    p.add_argument("--dps",type=int,default=15)
    args=p.parse_args()
    path=OUT/(args.label+".json")
    if path.exists():raise FileExistsError("Use a new label; results are not overwritten.")
    study=dict(generated_utc=datetime.now(timezone.utc).isoformat(),
      source="https://arxiv.org/html/2411.08090v1",equations=["6.26","A.14","B.21-B.24"],
      scope="Conditional fixed external label; exact finite-mass kernel with finite spectral cutoff; no packet average.",
      approximation="Numerical quadrature, no interval certification; G->0 and free matter are source model assumptions.",
      versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
      code_sha256={n:sha(ROOT/"code"/n) for n in ["run_crossing_conditional.py","observer_crossing.py","quantized_observer.py"]},
      input_sha256={"sources/kolchmeyer-liu-2411.08090v1.html":sha(ROOT/"sources/kolchmeyer-liu-2411.08090v1.html")},
      rows=[],status="running")
    OUT.mkdir(parents=True,exist_ok=True)
    def save():path.write_text(json.dumps(study,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    save()
    try:
        for n in args.orders:
            row=evaluate_order(n,args.upper,args.sigma,args.dps)
            study["rows"].append(row);save()
            print(json.dumps(dict(completed_order=n,values=row["values"])),flush=True)
        study["status"]="completed";save()
    except BaseException as exc:
        study["status"]="failed";study["exception"]=repr(exc);save();raise
if __name__=="__main__":main()
