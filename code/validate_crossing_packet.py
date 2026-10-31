"""Validate the actual finite packet, including frozen raw-point replay.

Tolerances were set before the 16th-order results. Refinement differences
are empirical and are never labeled rigorous quadrature/infinite-tail bounds.
"""
from pathlib import Path
from datetime import datetime,timezone
import json,gzip,math,sqlite3
import numpy as np
from crossing_point_cache import sha,configuration,point_key,digest_text
from packet_quadrature import cartesian_tasks,rephase
from observer_crossing_regular import RegularCrossingEvaluator
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/stage3/packet-crossing"
INNER_TOL=2e-6
OUTER_TOL=2e-7


def z(row):return complex(row["real"],row["imag"])
def cj(value):return dict(real=float(value.real),imag=float(value.imag))


def main():
    files=["packet-crossing-study.json","setup-validation.json","positive-validation.json","tail-validation.json"]
    study,setup,positive,tail=[json.loads((OUT/n).read_text(encoding="utf-8-sig")) for n in files]
    if study["status"]!="completed":raise RuntimeError("Wait for the live run to complete")
    if setup["status"]!="passed" or positive["status"]!="passed" or tail["status"]!="passed":raise RuntimeError("Prerequisites not passed")
    archive_path=ROOT/study["frozen_point_archive"]
    if sha(archive_path)!=study["input_sha256"][study["frozen_point_archive"]]:raise ValueError("Frozen raw archive hash mismatch")
    with gzip.open(archive_path,"rt",encoding="utf-8") as f:archive=json.load(f)
    records={r["key"]:r for r in archive["records"]}
    checks=[]
    def check(name,passed,**info):checks.append(dict(name=name,passed=bool(passed),**info))
    def close(name,a,b,tol):check(name,abs(a-b)<=tol,error=float(abs(a-b)),tolerance=tol)
    check("Raw point identities unique",len(records)==len(archive["records"])==study["unique_points"])
    check("Current kernel/cache identity unchanged",study["point_configuration"]==configuration() and archive["configuration"]==configuration())
    check("Every raw identity and key agree",all(point_key(r["identity"])==k for k,r in records.items()))
    check("All recorded raw values finite",all(math.isfinite(float(r["conditional"][v])) for r in records.values() for v in ("real","imag")))
    cache_db=sqlite3.connect((OUT/"point-cache.sqlite").resolve().as_uri()+"?mode=ro",uri=True)
    cache_mismatches=[]
    for key,record in records.items():
        stored=cache_db.execute("SELECT payload,payload_sha256 FROM points WHERE key=?",(key,)).fetchone()
        if stored is None or digest_text(stored[0])!=stored[1] or json.loads(stored[0])!=record:cache_mismatches.append(key)
    cache_db.close()
    check("Persistent cache agrees with every frozen raw point",not cache_mismatches,checked=len(records),mismatch_count=len(cache_mismatches))
    rows={(r["outer_order"],r["inner_order"]):r for r in study["rows"]}
    check("All planned rule pairs present",set(rows)=={(2,12),(3,12),(2,16),(3,16)})
    if set(rows)!={(2,12),(3,12),(2,16),(3,16)}:raise ValueError("Incomplete grid plan")
    all_keys=set();replayed={};task_by_pair={}
    for pair,row in rows.items():
        tasks=cartesian_tasks(pair[1],row["outer_rule"],study["point_configuration"])
        task_by_pair[pair]=tasks;keys={t["key"] for t in tasks};all_keys.update(keys)
        check("Raw coverage "+str(pair),keys==set(row["point_keys"]) and keys<=set(records))
        total,conditional,_=rephase(tasks,records,row["outer_rule"])
        replayed[pair]=total
        for i,(value,old) in enumerate(zip(total,row["values"])):close(f"Frozen raw replay {pair} tindex={i}",value,z(old),2e-16)
    check("No missing or extraneous exported points",all_keys==set(records))
    node_errors={};weighted_inner={}
    for outer in (2,3):
        low,high=rows[(outer,12)],rows[(outer,16)]
        errors=[]
        for j,(lo,hi) in enumerate(zip(low["conditional_by_node"],high["conditional_by_node"])):
            check(f"Same external node for internal refinement {outer}/{j}",lo["d"]==hi["d"] and lo["weight"]==hi["weight"])
            err=[abs(z(b)-z(a)) for a,b in zip(lo["values"],hi["values"])]
            errors.append(err)
            for i,e in enumerate(err):check(f"Internal conditional refinement {outer}/{j}/{i}",e<=INNER_TOL,error=e,tolerance=INNER_TOL)
        ws=np.array(list(map(float,high["outer_rule"]["weights"])))
        weighted_inner[outer]=(ws @ np.array(errors)).tolist();node_errors[outer]=errors
        for i in range(4):close(f"Internal full-packet refinement {outer}/{i}",replayed[(outer,12)][i],replayed[(outer,16)][i],INNER_TOL)
    outer_changes={}
    for inner in (12,16):
        errors=[]
        for i in range(4):
            a,b=replayed[(2,inner)][i],replayed[(3,inner)][i]
            close(f"External G2 versus G3 at inner={inner} tindex={i}",a,b,OUTER_TOL);errors.append(abs(a-b))
        outer_changes[inner]=errors
    noise=(np.array(weighted_inner[2])+np.array(weighted_inner[3])).tolist()
    fine=rows[(3,16)];values=replayed[(3,16)];baseline=positive["final"]
    nc,dc=baseline["N_cut"],baseline["D_cut"]
    for i,value in enumerate(values):check("Packet Cauchy Schwarz "+str(i),abs(value)<=nc+INNER_TOL,absolute_value=abs(value),N_cut=nc,tolerance=INNER_TOL)
    close("Zero time packet reality",values[0].imag,0.,1e-16)
    m4=dc+nc+values[0].real
    check("Same-field fourth-order form",dc-INNER_TOL<=m4<=dc+2*nc+INNER_TOL,D_cut=dc,M4=m4,upper=dc+2*nc,tolerance=INNER_TOL)
    comm=[2*(nc-x.real) for x in values]
    for i,v in enumerate(comm):check("Two-step difference form "+str(i),v>=-2*INNER_TOL and v<=4*nc+2*INNER_TOL,norm_squared=v,tolerance=2*INNER_TOL)
    reverse,_,_=rephase(task_by_pair[(3,16)],records,fine["outer_rule"],times=(0.,-.5,-1.,-2.))
    for i,(a,b) in enumerate(zip(reverse,values)):close("Time reflection from frozen nodes "+str(i),a,b.conjugate(),1e-14)
    old=json.loads((ROOT/"results/stage3/fourpoint/conditional-crossing-cartesian.json").read_text())
    for order in (12,16):
        previous=next(r for r in old["rows"] if r["order"]==order)
        current=rows[(3,order)]["conditional_by_node"][1]
        for i,(a,b) in enumerate(zip(previous["values"],current["values"])):close(f"Center node regression {order}/{i}",z(a),z(b),1e-12)
    # High precision at actual final nodes: largest cancellation and largest
    # positive/negative damped weighted contribution. No inferred zero values.
    tasks=task_by_pair[(3,16)]
    def weighted_score(task):
        from quantized_observer import mass
        labels=list(map(float,task["identity"]["labels"]))
        ma,mb,mc,md=map(float,map(mass,labels));ds=[md-mc,mc-mb,mb-ma,ma-md]
        return float(records[task["key"]]["conditional"]["real"])*math.exp(-sum(x*x for x in ds)/2)*task["inner_weight"]*float(fine["outer_rule"]["weights"][task["outer_index"]])
    chosen=[fine["worst_point_key"],max(tasks,key=weighted_score)["key"],min(tasks,key=weighted_score)["key"]]
    high_precision=[]
    ev=RegularCrossingEvaluator(35,15)
    for i,key in enumerate(dict.fromkeys(chosen)):
        record=records[key];identity=record["identity"]
        result=ev.evaluate(identity["nu"],*identity["labels"])
        with ev.ctx.workdps(result.working_dps):ref=complex(result.hh_cross/ev.rho_hh(identity["labels"][-1]))
        raw=complex(float(record["conditional"]["real"]),float(record["conditional"]["imag"]))
        close("Actual packet-node precision "+str(i),raw,ref,5e-14*max(1e-8,abs(ref)))
        high_precision.append(dict(key=key,identity=identity,old=cj(raw),high=cj(ref),working_dps=result.working_dps))
    failed=sum(not c["passed"] for c in checks)
    diagnostics=dict(values=[cj(v) for v in values],times=[0.,.5,1.,2.],positive_baseline=baseline,
        M4_cut=m4,commutator_form_by_time=comm,far_cap_analytic_bounds=tail["strict_coarse_bounds"],
        node_internal_changes=node_errors,weighted_internal_changes=weighted_inner,
        external_changes_at_fixed_inner=outer_changes,combined_internal_empirical_noise=noise,
        outer_changes_below_internal_noise=[e<=n for e,n in zip(outer_changes[16],noise)],
        inner_tolerance=INNER_TOL,outer_tolerance=OUTER_TOL,high_precision=high_precision,
        scope="All errors are empirical diagnostics of this normalized finite packet and finite spectral cutoff. Far-cap analytic bounds are recorded separately; finite-cap positive quadrature and the total Lambda=3 error remain empirical. Full CG physical-domain identification remains separate.")
    inputs={("results/stage3/packet-crossing/"+n):sha(OUT/n) for n in files}
    inputs[study["frozen_point_archive"]]=sha(archive_path)
    inputs["results/stage3/fourpoint/conditional-crossing-cartesian.json"]=sha(ROOT/"results/stage3/fourpoint/conditional-crossing-cartesian.json")
    result=dict(status="passed" if not failed else "failed",checks=len(checks),failed=failed,check_details=checks,diagnostics=diagnostics,
        generated_utc=datetime.now(timezone.utc).isoformat(),
        code_sha256={n:sha(ROOT/"code"/n) for n in ["validate_crossing_packet.py","crossing_point_cache.py","packet_quadrature.py","observer_crossing_regular.py","observer_crossing.py","quantized_observer.py"]},input_sha256=inputs)
    target=OUT/"validation.json"
    if target.exists():raise FileExistsError("Keep earlier validation before rerun")
    target.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(dict(status=result["status"],checks=len(checks),failed=failed,values=diagnostics["values"],commutator_form=comm)),flush=True)
    if failed:raise SystemExit(1)
if __name__=="__main__":main()