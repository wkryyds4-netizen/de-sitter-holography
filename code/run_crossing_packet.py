"""Finite normalized external packet, exact finite-mass crossed kernel.

Every table row specifies weighted outer order and Cartesian inner order.
Raw points are cached and exported with identities/precision/provenance;
time dependence is summed afterwards. Output completion is not validation.
"""
from pathlib import Path
from datetime import datetime,timezone
from concurrent.futures import ProcessPoolExecutor,as_completed
import argparse,json,gzip,time,platform,re
import numpy as np,scipy,mpmath
from crossing_point_cache import PointCache,configuration,point_evaluate,sha
from packet_quadrature import gauss_packet,cartesian_tasks,rephase
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/stage3/packet-crossing"


def cj(z):return dict(real=float(z.real),imag=float(z.imag))


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--pairs",nargs="+",default=["2:12","3:12","2:16","3:16"])
    p.add_argument("--workers",type=int,default=16)
    p.add_argument("--dps",type=int,default=15)
    p.add_argument("--label",default="packet-crossing-study")
    a=p.parse_args()
    if not re.fullmatch(r"[a-zA-Z0-9_-]+",a.label):raise ValueError("Plain label required")
    pairs=[tuple(map(int,x.split(":"))) for x in a.pairs]
    if any(len(x)!=2 or not 1<=x[0]<=12 or not 2<=x[1]<=100 for x in pairs):raise ValueError("Invalid order pair")
    OUT.mkdir(parents=True,exist_ok=True)
    output=OUT/(a.label+".json");archive=OUT/(a.label+"-points.json.gz")
    if output.exists() or archive.exists():raise FileExistsError("Preserve old runs; use a new label")
    cfg=configuration(a.dps,8)
    study=dict(schema="packet-crossing-study-v1",status="running",generated_utc=datetime.now(timezone.utc).isoformat(),
        parameters=dict(nu="0.6",packet_lo="0.7",packet_hi="0.9",internal_lo=.2,internal_hi=3.,sigma_proper=1.,
            times_left_to_right=[[t,0,t,0] for t in (0.,.5,1.,2.)]),
        scope="dS2 L=1 G->0 free principal series, normalized finite sin^4 external packet, every mass segment projected",
        numerical_limits="Empirical order refinements; no interval certification or infinite-tail certification",
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,mpmath=mpmath.__version__),
        point_configuration=cfg,workers=a.workers,requested_pairs=pairs,rows=[],
        code_sha256={n:sha(ROOT/"code"/n) for n in ["run_crossing_packet.py","packet_quadrature.py",
            "crossing_point_cache.py","observer_crossing_regular.py","observer_crossing.py","quantized_observer.py"]},
        input_sha256={"sources/kolchmeyer-liu-2411.08090v1.html":sha(ROOT/"sources/kolchmeyer-liu-2411.08090v1.html")})
    def save():output.write_text(json.dumps(study,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    save();all_records={};global_start=time.perf_counter()
    try:
        with PointCache(OUT/"point-cache.sqlite") as cache:
            for outer,inner in pairs:
                started=time.perf_counter();rule=gauss_packet(outer)
                tasks=cartesian_tasks(inner,rule,cfg)
                unique={t["key"]:t["identity"] for t in tasks}
                missing=[];records={}
                for key,identity in unique.items():
                    record=all_records.get(key)
                    if record is None:record=cache.get(identity)
                    if record is None:missing.append(identity)
                    else:records[key]=record
                cached=len(records)
                # Large labels were the slow tail of the earlier ordered map.
                # Put those jobs first and collect completed futures out of order;
                # phase summation later uses deterministic quadrature order.
                missing.sort(key=lambda ident:sum(map(float,ident["labels"])),reverse=True)
                print(json.dumps(dict(outer=outer,inner=inner,points=len(unique),cached=cached,new=len(missing))),flush=True)
                if missing:
                    with ProcessPoolExecutor(max_workers=a.workers) as pool:
                        futures=[pool.submit(point_evaluate,identity) for identity in missing]
                        for done,future in enumerate(as_completed(futures),1):
                            record=future.result();cache.put(record);records[record["key"]]=record
                            if done%128==0 or done==len(missing):
                                print(json.dumps(dict(outer=outer,inner=inner,evaluated=done,new_total=len(missing),
                                    seconds=time.perf_counter()-started)),flush=True)
                cache.commit();all_records.update(records)
                total,conditional,absolute=rephase(tasks,records,rule)
                worst=max(records.values(),key=lambda r:float(r["cancellation_ratio"]))
                row=dict(outer_order=outer,inner_order=inner,outer_rule=rule,points=len(tasks),
                    cache_hits=cached,new_points=len(missing),values=[cj(z) for z in total],
                    conditional_by_node=[dict(d=d,weight=w,values=[cj(z) for z in vals],absolute_sum=ab)
                        for d,w,vals,ab in zip(rule["nodes"],rule["weights"],conditional,absolute)],
                    point_keys=sorted(records),worst_point_key=worst["key"],maximum_cancellation=worst["cancellation_ratio"],
                    maximum_working_dps=max(r["working_dps"] for r in records.values()),
                    elapsed_seconds=time.perf_counter()-started)
                study["rows"].append(row);save()
                print(json.dumps(dict(completed_pair=[outer,inner],values=row["values"])),flush=True)
        with gzip.open(archive,"xt",encoding="utf-8") as stream:
            json.dump(dict(schema="frozen-crossing-points-v1",configuration=cfg,
                records=[all_records[k] for k in sorted(all_records)]),stream,sort_keys=True,allow_nan=False)
        study["input_sha256"][archive.relative_to(ROOT).as_posix()]=sha(archive)
        study["frozen_point_archive"]=archive.relative_to(ROOT).as_posix()
        study["unique_points"]=len(all_records)
        study["elapsed_seconds"]=time.perf_counter()-global_start
        study["status"]="completed";save()
        print(json.dumps(dict(status="completed",rows=len(study["rows"]),unique_points=len(all_records),seconds=study["elapsed_seconds"])),flush=True)
    except BaseException as exc:
        study["status"]="failed";study["exception"]=repr(exc);save();raise


if __name__=="__main__":main()