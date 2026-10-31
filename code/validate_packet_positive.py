"""Refine the cheap full-packet positive baseline and compare independent coordinates."""
from pathlib import Path
from datetime import datetime,timezone
import json,platform
import numpy as np
from scipy.integrate import quad
from packet_positive import positive_packet,positive_at_nodes,average_positive
from packet_quadrature import gauss_packet,inner_rule,packet_density
from quantized_observer import right_kernel,mass
from run_crossing_conditional import positive_norms
from crossing_point_cache import sha
ROOT=Path(__file__).resolve().parents[1]


def main():
    checks=[]
    def check(name,passed,**info):checks.append(dict(name=name,passed=bool(passed),**info))
    def close(name,a,b,tol):check(name,abs(a-b)<=tol,error=float(abs(a-b)),tolerance=tol)
    settings=[(24,64,12.),(32,64,12.),(32,96,12.),(48,96,16.)]
    rows=[dict(outer_order=o,inner_order=i,cap=c,values=positive_packet(o,i,c)) for o,i,c in settings]
    fine=rows[-1]["values"]
    fields=["C2_cut","D_cut","N_cut","epsilon_to_cap","cross_tail_bound_estimate"]
    for row in rows[:-1]:
        for name in fields:close(f"Positive grid {row['outer_order']}/{row['inner_order']}/{row['cap']} {name}",row["values"][name],fine[name],5e-12)
    for d in (.7,.8,.9):
        vector=positive_at_nodes([d],96,12.)
        adaptive=positive_norms(.6,d,.2,3.,1.,cap=12.)
        for name in ("N_cut","epsilon_to_cap","inner_tail","outer_tail"):
            close(f"Independent adaptive at d={d} {name}",vector[name][0],adaptive[name],5e-12)
        def physical_kernel(m):
            s=np.sqrt(m*m-.25)
            return float(right_kernel(.6,s,d))*np.exp(-(m-float(mass(d)))**2)*m/s
        physical,error=quad(physical_kernel,float(mass(.2)),float(mass(3.)),epsabs=1e-13,epsrel=1e-12)
        close(f"Independent physical mass C2 at d={d}",vector["C2_cut"][0],physical,5e-13)
    ds,weights=inner_rule(48,.7,.9);weights*=packet_density(ds)
    c2s=positive_at_nodes(ds,96,16.)["C2_cut"]
    variance=float(weights @ (c2s-fine["C2_cut"])**2)
    close("Stable external conditional variance",variance,fine["conditional_two_point_variance"],5e-14)
    check("Finite packet normalization",abs(fine["external_quadrature_mass"]-1)<1e-13)
    check("Positive packet direct tails",fine["inner_tail"]>0 and fine["outer_tail"]>0)
    check("Common external mass prevents factorization",variance>1e-8,variance=variance)
    low_rules=[]
    for n in (2,3):
        rule=gauss_packet(n)
        values=average_positive(list(map(float,rule["nodes"])),list(map(float,rule["weights"])),96,16.)
        low_rules.append(dict(order=n,rule=rule,values=values))
    failed=sum(not c["passed"] for c in checks)
    result=dict(status="passed" if not failed else "failed",checks=len(checks),failed=failed,check_details=checks,
        generated_utc=datetime.now(timezone.utc).isoformat(),parameters=dict(nu=.6,packet=[.7,.9],lower=.2,upper=3.,sigma_proper=1.),
        rows=rows,weighted_gauss_baselines=low_rules,final=fine,
        scope="Actual normalized finite packet; independent positive finite-cap baselines. Not an infinite-tail certification.",
        code_sha256={n:sha(ROOT/"code"/n) for n in ["validate_packet_positive.py","packet_positive.py","packet_quadrature.py","run_crossing_conditional.py","quantized_observer.py"]})
    output=ROOT/"results/stage3/packet-crossing/positive-validation.json"
    if output.exists():raise FileExistsError("Preserve previous validation")
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(dict(status=result["status"],checks=len(checks),failed=failed,final=fine)),flush=True)
    if failed:raise SystemExit(1)
if __name__=="__main__":main()