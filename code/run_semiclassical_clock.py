"""Compare finite quantized observers with the dS2 common-clock limit."""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,platform
import numpy as np,scipy
from scipy.integrate import quad
from quantized_observer import DoublePacket,mass,right_kernel
from observer_fourpoint import FourPointQuadrature
from semiclassical_clock import stationary_transition,packet_survival,ClockLimit,energy_tail_moments

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/stage3/fourpoint"
def serial(v):
    if isinstance(v,np.generic):return serial(v.item())
    if isinstance(v,complex):return {"real":v.real,"imag":v.imag}
    if isinstance(v,dict):return {k:serial(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [serial(x) for x in v]
    return v

def main():
    checks=[]
    def check(n,b,**kw):checks.append(dict(name=n,passed=bool(b),**serial(kw)))
    def close(n,x,y,tol):
        e=float(abs(x-y));check(n,e<=tol,error=e,tolerance=tol)
    packet=DoublePacket();nu=1.;sigma=.5
    for q in (0.,2.,2.5,3.7,5.,8.5,10.,12.):
        ref=sum(quad(lambda p:float(packet.probability(p)),max(q,a),b,epsabs=1e-13)[0]
                for a,b in packet.intervals if b>q)
        close(f"survival vs independent packet integral {q}",float(packet_survival(packet,q)),ref,2e-12)
    for d in (.2,1.,2.):
        close(f"limiting detailed balance {d}",stationary_transition(d,nu),
              np.exp(-2*np.pi*d)*stationary_transition(-d,nu),1e-14)
    limit=ClockLimit(packet,nu,sigma)
    fine=ClockLimit(packet,nu,sigma,order=40,p_order=64,extra=18)
    times=[[0.,0.,0.,0.],[1.,0.,1.,0.],[2.,0.,2.,0.]]
    lv=[limit.values(t) for t in times]
    fv=[fine.values(t) for t in times]
    for i in range(len(times)):
        for k in ("C2","AABB","ABBA"):close(f"limit quadrature {i,k}",lv[i][k],fv[i][k],2e-11)
    tails=energy_tail_moments(packet,nu,sigma,order=32,radius=14)
    tails_fine=energy_tail_moments(packet,nu,sigma,order=48,radius=14)
    tails_extended=energy_tail_moments(packet,nu,sigma,order=48,radius=18)
    for k in ("C2","AABB","ABBA"):
        close(f"energy-tail refinement {k}",tails[k],tails_fine[k],2e-9)
        close(f"energy-tail cutoff extension {k}",tails_extended[k],tails_fine[k],2e-9)
        close(f"independent common-clock survival formula {k}",tails_fine[k],fv[0][k],2e-9)
    mass_rows=[]
    for m0 in (2.,5.,20.,80.):
        calc=FourPointQuadrature(packet,m0,nu,sigma,True,p_order=48,mass_order=28,tail_extra=14)
        values=[calc.values(t) for t in times]
        c2=calc.two_point(0)
        errors=[]
        for i,row in enumerate(values):
            row["C2"]=c2
            errors.append({k:float(abs(row[k]-fv[i][k])) for k in ("C2","AABB","ABBA")})
        mass_rows.append(dict(m0=m0,values=values,absolute_differences_from_limit=errors))
        print(m0,errors[0],flush=True)
    for k in ("C2","AABB","ABBA"):
        errors=[row["absolute_differences_from_limit"][0][k] for row in mass_rows]
        check("observed zero-time convergence "+k,all(x>y for x,y in zip(errors,errors[1:])),errors=errors,
              scope="empirical sampled monotonicity, not a uniform rate theorem")
    kernel_rows=[]
    for mi in (5.,20.,80.):
        bi=np.sqrt(mi*mi-.25)
        for delta in (-1.,-.2,.7,1.5):
            mf=mi+delta;af=np.sqrt(mf*mf-.25)
            value=float(right_kernel(nu,af,bi)*mf/af)
            ref=float(stationary_transition(delta,nu))
            kernel_rows.append(dict(initial_mass=mi,delta=delta,finite=value,limit=ref,error=value-ref))
    for delta in (-1.,-.2,.7,1.5):
        es=[abs(r["error"]) for r in kernel_rows if r["delta"]==delta]
        check(f"pointwise mass-kernel convergence {delta}",all(x>y for x,y in zip(es,es[1:])),errors=es)
    scripts=["semiclassical_clock.py","run_semiclassical_clock.py","observer_fourpoint.py","quantized_observer.py"]
    meta=dict(generated_utc=datetime.now(timezone.utc).isoformat(),
              code_sha256={n:hashlib.sha256((ROOT/"code"/n).read_bytes()).hexdigest() for n in scripts},
              versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
              source="https://arxiv.org/html/2411.08090v1",equations=["6.18","6.21","6.23","6.24"],
              scope="Fixed proper-time Gaussian, fixed compact flat-clock packet, two noncrossing flavor tasks in dS2. Does not establish finite-time-recoil ABAB or cross-dimensional equivalence.",
              numerical_limits="Quadrature refinement and truncation comparisons, not interval certification.")
    study=dict(**meta,nu=nu,sigma_proper=sigma,packet={"c1":3,"c2":9,"width":1},
               times_left_to_right=times,limit_values=fv,energy_tail=tails_fine,energy_tail_extended=tails_extended,
               finite_mass_rows=mass_rows,kernel_rows=kernel_rows)
    failed=sum(not x["passed"] for x in checks)
    validation=dict(**meta,status="passed" if failed==0 else "failed",checks=len(checks),failed=failed,check_details=checks)
    OUT.mkdir(parents=True,exist_ok=True)
    for n,d in [("clock-limit-study.json",study),("clock-limit-validation.json",validation)]:
        (OUT/n).write_text(json.dumps(serial(d),ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(dict(checks=len(checks),failed=failed,limit_zero=serial(fv[0]),tails=tails_fine)),flush=True)
    if failed:raise SystemExit(1)

if __name__=="__main__":main()
