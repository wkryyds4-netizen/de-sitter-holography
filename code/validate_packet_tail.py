"""Check algebraic constants and evaluate the proved conservative far-tail envelope.

The proof is in stage3-packet-tail-envelope.md. mpmath evaluation of its sharp
expression is not interval arithmetic. Coarse bounds use exact integer or
rational inequalities together with the analytic inequalities proved there.
"""
from pathlib import Path
from datetime import datetime,timezone
from fractions import Fraction
import json,mpmath,math
from packet_positive import positive_packet
from crossing_point_cache import sha
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/stage3/packet-crossing"


def envelope(r,d,sigma=1,dps=70):
    c=mpmath.mp.clone();c.dps=dps
    s=c.mpf(sigma);d=c.mpf(d);r=c.mpf(r)
    C=c.gamma(c.mpf(1)/4)**8/(16*c.pi**5)
    M=c.sqrt(d*d+c.mpf(1)/4);K=c.sqrt(r*r+c.mpf(1)/4)
    mu=M+c.pi/s**2;Q=2*c.pi*M+c.pi**2/s**2
    if K<=mu:raise ValueError("This choice of Chernoff parameter requires K>mu")
    pref=c.pi*C*C/s**2*c.cosh(c.pi*d)**2*c.exp(Q)
    bx=K*K+3/(4*s*s)
    by=K*K/2+(K+mu)**2/8+3/(4*s*s)
    total=pref*(bx*c.exp(-s*s*(K-mu)**2)+by*c.exp(-s*s*(K-mu)**2/2))
    norm=pref*(mu*mu+3/(4*s*s))
    return dict(B_R=c.nstr(total,45),full_norm_envelope=c.nstr(norm,45),mu=c.nstr(mu,45),K=c.nstr(K,45),C0=c.nstr(C,45))


def main():
    checks=[]
    def check(name,passed,**info):checks.append(dict(name=name,passed=bool(passed),**info))
    mu=Fraction(26,25)+Fraction(22,7)
    check("C0 coarse rational bound",Fraction(4**8,16*3**5)<17)
    check("Q rational bound",2*Fraction(22,7)*Fraction(26,25)+Fraction(22,7)**2<17)
    check("Chernoff exponent at R20",(20-mu)**2/2>125)
    check("Elementary e lower bound",sum(Fraction(1,math.factorial(k)) for k in range(5))>Fraction(27,10))
    check("Exact integer exponential comparison",795328*10**140<27**102)
    check("Exact norm prefactor bound",1156*3**23*Fraction(103,4)<3*10**15)
    check("Squared conservative pairing bound",4*Fraction(3*10**15,10**38)<Fraction(4,10**22))
    context=mpmath.mp.clone();context.dps=60
    evaluations=[]
    for r in (16,20):
        previous=0.
        for d in ("0.7","0.8","0.9"):
            row=envelope(r,d);row.update(spectral_cutoff=r,d=d,sigma=1)
            fine=envelope(r,d,dps=90)
            a,b=context.mpf(row["B_R"]),context.mpf(fine["B_R"])
            check(f"70 to 90 digit envelope {r}/{d}",abs(a-b)<=abs(b)*context.mpf("1e-40"))
            check(f"Sample monotonicity {r}/{d}",float(a)>previous)
            previous=float(a);evaluations.append(row)
    cap20=positive_packet(32,96,20.)
    refined=positive_packet(48,128,20.)
    old=json.loads((OUT/"positive-validation.json").read_text())["final"]
    for field in ("C2_cut","D_cut","N_cut","epsilon_to_cap","cross_tail_bound_estimate"):
        error=abs(cap20[field]-refined[field])
        check("Cap20 positive quadrature refinement "+field,error<=5e-12,error=error,tolerance=5e-12)
        error=abs(refined[field]-old[field])
        check("Cap16 to cap20 positive baseline "+field,error<=5e-12,error=error,tolerance=5e-12)
    failed=sum(not c["passed"] for c in checks)
    result=dict(status="passed" if not failed else "failed",checks=len(checks),failed=failed,check_details=checks,
        generated_utc=datetime.now(timezone.utc).isoformat(),sharp_expression_evaluations=evaluations,
        strict_coarse_bounds=dict(sigma=1,external_d_max=.9,spectral_cutoff=20,
            epsilon_beyond_cap_upper="1e-38",full_N_upper="3e15",all_real_time_cross_pairing_tail_upper="2e-11"),
        cap20_positive_estimate=cap20,cap20_refined_positive_estimate=refined,
        scope="The analytic far-tail bound is strict; positive finite-cap quadrature and its differences are empirical. No strict total error for the Lambda=3 computation is claimed.",
        code_sha256={n:sha(ROOT/"code"/n) for n in ["validate_packet_tail.py","packet_positive.py","packet_quadrature.py","quantized_observer.py"]},
        input_sha256={"results/stage3/packet-crossing/positive-validation.json":sha(OUT/"positive-validation.json"),
            "notes/stage3-packet-tail-envelope.md":sha(ROOT/"notes/stage3-packet-tail-envelope.md")})
    target=OUT/"tail-validation.json"
    if target.exists():raise FileExistsError("Keep earlier validation")
    target.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(dict(status=result["status"],checks=len(checks),failed=failed,bounds=result["strict_coarse_bounds"],positive=refined)),flush=True)
    if failed:raise SystemExit(1)
if __name__=="__main__":main()