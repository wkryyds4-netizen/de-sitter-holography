"""Independent moment, cache-identity and phase-symmetry checks before the run."""
from pathlib import Path
from datetime import datetime,timezone
import json,math,copy
import numpy as np
from scipy.integrate import quad
from packet_quadrature import gauss_packet,packet_density,cartesian_tasks,rephase,inner_rule
from crossing_point_cache import configuration,point_identity,point_key,point_evaluate,PointCache,sha
from observer_crossing import CrossingEvaluator
from quantized_observer import mass
ROOT=Path(__file__).resolve().parents[1]


def main():
    checks=[]
    def check(name,passed,**details):checks.append(dict(name=name,passed=bool(passed),**details))
    def close(name,a,b,tol):check(name,abs(a-b)<=tol,error=float(abs(a-b)),tolerance=tol)
    pi=math.pi
    for n in (2,3):
        rule=gauss_packet(n)
        ds=np.array(list(map(float,rule["nodes"])));ws=np.array(list(map(float,rule["weights"])))
        xs=(ds-.8)/.1
        check("Positive normalized rule "+str(n),np.all(ws>0) and np.all(ds>.7) and np.all(ds<.9))
        for k in range(2*n):
            direct=quad(lambda x:x**k*(.5+2*math.cos(pi*x)/3+math.cos(2*pi*x)/6),-1,1,epsabs=1e-14)[0]
            close(f"Independent centered moment n={n} k={k}",np.dot(ws,xs**k),direct,3e-15)
        if n==2:close("Analytic G2 displacement",abs(ds[0]-.8),.1*math.sqrt(1/3-5/(2*pi*pi)),2e-16)
        else:
            mu2=1/3-5/(2*pi*pi);mu4=1/5-5/pi**2+63/(2*pi**4)
            close("Analytic G3 displacement",abs(ds[0]-.8),.1*math.sqrt(mu4/mu2),3e-16)
    close("Independent packet normalization",quad(lambda d:float(packet_density(d)),.7,.9,epsabs=1e-14)[0],1.,2e-15)
    cfg=configuration();identity=point_identity(cfg,"0.6",["0.61","0.71","0.83","0.8"])
    record=point_evaluate(identity)
    old=CrossingEvaluator(35,15);ref=complex(old.conditional_cross("0.6",*identity["labels"]))
    raw=complex(float(record["conditional"]["real"]),float(record["conditional"]["imag"]))
    close("Cached raw point versus old high precision formula",raw,ref,2e-16)
    with PointCache(":memory:") as cache:
        check("Initial miss",cache.get(identity) is None)
        cache.put(record);check("Raw record round trip",cache.get(identity)==record)
        for field in ("dps","guard_digits","source_sha256"):
            other=copy.deepcopy(identity)
            other["configuration"][field]=str(other["configuration"][field])+"_changed"
            check("Cache invalidation "+field,cache.get(other) is None and point_key(other)!=record["key"])
        cache.db.execute("UPDATE points SET payload=payload||' ' WHERE key=?",(record["key"],))
        try:cache.get(identity);detected=False
        except ValueError:detected=True
        check("Corrupted payload rejected",detected)
    # A separately defined symmetric raw kernel checks weights and phases.
    # Direct full cube uses all nodes, without triangular multiplicities.
    rule=gauss_packet(3);tasks=cartesian_tasks(3,rule,cfg);records={}
    for task in tasks:
        vals=list(map(float,task["identity"]["labels"]))
        value=1/(1+sum(vals))
        records[task["key"]]=dict(conditional=dict(real=str(value),imag="0"))
    times=(0.,.5,1.,2.,-2.)
    total,_,_=rephase(tasks,records,rule,times)
    nodes,weights=inner_rule(3);direct=[]
    for t in times:
        answer=0j
        for d,pw in zip(map(float,rule["nodes"]),map(float,rule["weights"])):
            for a,aw in zip(nodes,weights):
                for b,bw in zip(nodes,weights):
                    for c,cw in zip(nodes,weights):
                        md,ma,mb,mc=map(float,map(mass,(d,a,b,c)))
                        differences=np.array([md-mc,mc-mb,mb-ma,ma-md])
                        answer+=pw*aw*bw*cw/(1+a+b+c+d)*np.exp(-np.dot(differences,differences)/2)*np.exp(-1j*np.dot(differences,[t,0,t,0]))
        direct.append(answer)
    for i,(a,b) in enumerate(zip(total,direct)):close("Full cube versus folded phase "+str(i),a,b,3e-14)
    close("Time reversal synthetic kernel",total[-1],total[-2].conjugate(),1e-14)
    failed=sum(not c["passed"] for c in checks)
    result=dict(status="passed" if failed==0 else "failed",checks=len(checks),failed=failed,check_details=checks,
        generated_utc=datetime.now(timezone.utc).isoformat(),actual_point=record,
        code_sha256={n:sha(ROOT/"code"/n) for n in ["validate_packet_setup.py","packet_quadrature.py","crossing_point_cache.py","observer_crossing_regular.py","observer_crossing.py","quantized_observer.py"]},
        input_sha256={"sources/kolchmeyer-liu-2411.08090v1.html":sha(ROOT/"sources/kolchmeyer-liu-2411.08090v1.html")})
    out=ROOT/"results/stage3/packet-crossing/setup-validation.json";out.parent.mkdir(exist_ok=True)
    if out.exists():raise FileExistsError("Keep the previous validation")
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(dict(status=result["status"],checks=len(checks),failed=failed)),flush=True)
    if failed:raise SystemExit(1)
if __name__=="__main__":main()