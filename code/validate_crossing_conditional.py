"""Independent diagnostics for finite-cutoff conditional crossed integrals."""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,platform
import numpy as np
from scipy.integrate import quad
from observer_crossing import CrossingEvaluator
from observer_crossing_regular import RegularCrossingEvaluator
from quantized_observer import mass,right_kernel
from run_crossing_conditional import positive_norms,cj,sha
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/stage3/fourpoint"
TOL=2e-6
def z(row):return complex(row["real"],row["imag"])
def main():
    filenames=["conditional-crossing-pilot.json","conditional-crossing-refined.json","conditional-crossing-cartesian.json"]
    data={n:json.loads((OUT/n).read_text(encoding="utf-8-sig")) for n in filenames}
    for n,d in data.items():
        if d["status"]!="completed":raise RuntimeError("Not terminal: "+n)
    pilot=data[filenames[0]]["rows"];tri=data[filenames[1]]["rows"];cart=data[filenames[2]]["rows"]
    checks=[];diagnostics={}
    def check(name,passed,**info):checks.append(dict(name=name,passed=bool(passed),**info))
    def close(name,a,b,tol):
        e=abs(a-b);check(name,e<=tol,error=float(e),tolerance=tol)
    fine=cart[-1];tv=[z(v) for v in tri[-1]["values"]];cv=[z(v) for v in fine["values"]]
    for i in range(4):
        close("Cartesian order refinement "+str(i),cv[i],z(cart[-2]["values"][i]),TOL)
        close("Independent triangle versus Cartesian "+str(i),tv[i],cv[i],TOL)
    norms=fine["positive_norms"]
    nc=norms["N_cut"];ne=norms["N_to_cap"]
    for i,val in enumerate(cv):
        check("Finite-cutoff Cauchy Schwarz "+str(i),abs(val)<=nc+TOL,absolute_value=abs(val),N_cut=nc,tolerance=TOL)
    close("Zero-time reality",cv[0].imag,0.,1e-16)
    # Conditional two-point from adaptive s integral, independent of cross grids.
    nu,d,lo,upper,sigma=fine["nu"],fine["d"],fine["lower"],fine["upper"],fine["sigma_proper"]
    def k(a):return float(right_kernel(nu,a,d))*np.exp(-(sigma*(mass(a)-mass(d)))**2)
    hcut,herr=quad(k,lo,upper,epsabs=1e-13,epsrel=1e-12)
    dcut=hcut*hcut;m4=dcut+nc+cv[0].real
    commutator_norm=2*(nc-cv[0].real)
    check("Equal-time two-species commutator norm nonnegative",commutator_norm>=-2*TOL,norm_squared=commutator_norm,tolerance=2*TOL)
    check("Same-species four-order positive form bound",dcut-TOL<=m4<=dcut+2*nc+TOL,D_cut=dcut,M4_cut=m4,upper=dcut+2*nc,tolerance=TOL)
    extended=positive_norms(nu,d,lo,upper,sigma,cap=16.)
    for key in ("N_cut","N_to_cap","epsilon_to_cap","cross_tail_bound_estimate"):
        close("Positive norm cap 12 to 16 "+key,norms[key],extended[key],2e-12)
    check("Positive tails nonnegative",norms["outer_tail"]>=0 and norms["inner_tail"]>=0,outer_tail=norms["outer_tail"],inner_tail=norms["inner_tail"])
    # Independent higher precision checks at difficult actual quadrature nodes.
    ev=CrossingEvaluator(35,15);sample_rows=[]
    samples=pilot[-1]["diagnostic_samples"][:2]+[tri[-1]["worst_point"]]
    for i,s in enumerate(samples):
        labels=s["labels"];r=ev.evaluate(str(nu),*labels)
        ref=complex(r.hh_cross/ev.rho_hh(labels[-1]))
        old=z(s["conditional"])
        close("Actual quadrature node precision "+str(i),old,ref,2e-14*max(1e-8,abs(ref)))
        # A/C reflection independently evaluates the opposite triangle.
        mirror=ev.evaluate(str(nu),labels[2],labels[1],labels[0],labels[3])
        mr=complex(mirror.hh_cross/ev.rho_hh(labels[-1]))
        close("Opposite triangle reflection "+str(i),ref,mr,2e-25*max(1e-8,abs(ref)))
        sample_rows.append(dict(labels=labels,old=cj(old),high_precision=cj(ref),mirror=cj(mr),requested_dps=35,guard_digits=15))
    regular=RegularCrossingEvaluator(35,15)
    s=fine["worst_point"];labels=s["labels"]
    r=regular.evaluate(str(nu),*labels)
    with regular.ctx.workdps(r.working_dps):
        ref=complex(r.hh_cross/regular.rho_hh(labels[-1]))
    close("Cartesian difficult node precision",z(s["conditional"]),ref,2e-14*max(1e-8,abs(ref)))
    sample_rows.append(dict(labels=labels,old=s["conditional"],high_precision=cj(ref),method=r.method,requested_dps=35,guard_digits=15))
    diagnostics=dict(parameters=dict(nu=nu,d=d,lower=lo,upper=upper,sigma_proper=sigma),
        times_left_to_right=fine["times_left_to_right"],values=fine["values"],
        conditional_C2_cut=hcut,conditional_C2_adaptive_error=herr,D_cut=dcut,N_cut=nc,M4_cut=m4,
        crossing_to_nested_ratio=cv[0].real/nc,commutator_norm_squared_cut=commutator_norm,
        commutator_form_by_time=[dict(t=fine["times_left_to_right"][i][0],norm_squared=2*(nc-val.real)) for i,val in enumerate(cv)],positive_norms=norms,positive_norms_extended=extended,
        point_precision=sample_rows,
        order_history=[dict(method=method,order=r["order"],values=r["values"]) for method,rows in [("triangle",pilot+tri),("Cartesian",cart)] for r in rows],
        coarse_grid_differences=[dict(order=r["order"],max_difference_from_final=max(abs(z(x)-y) for x,y in zip(r["values"],cv))) for r in pilot],
        finite_domain_tolerance=TOL,
        scope="Conditional fiber d, not a normalized external wavepacket. Finite-domain quadrature tolerance is empirical. Tail is a positive finite-cap estimate inserted into a proved inequality; infinite-tail certification was not performed.")
    failed=sum(not c["passed"] for c in checks)
    output=dict(generated_utc=datetime.now(timezone.utc).isoformat(),status="passed" if not failed else "failed",
        checks=len(checks),failed=failed,check_details=checks,diagnostics=diagnostics,
        code_sha256={n:sha(ROOT/"code"/n) for n in ["validate_crossing_conditional.py","observer_crossing.py","observer_crossing_regular.py","run_crossing_conditional.py","quantized_observer.py"]},
        input_sha256={("results/stage3/fourpoint/"+n):sha(OUT/n) for n in filenames})
    target=OUT/"conditional-validation.json"
    if target.exists():raise FileExistsError("Preserve the prior validation before rerunning.")
    target.write_text(json.dumps(output,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(dict(status=output["status"],checks=len(checks),failed=failed,values=fine["values"],D_cut=dcut,N_cut=nc,M4_cut=m4)),flush=True)
    if failed:raise SystemExit(1)
if __name__=="__main__":main()
