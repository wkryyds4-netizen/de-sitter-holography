"""Exact-surface and heavy-mass diagnostics, preserving both precision levels."""
import hashlib
import json
import platform
import time
from pathlib import Path
import mpmath as mp
from ds3_racah_collision import physical_resolved,collision_data,omega_collision,omega_collision_derivative
from global_racah import omega_quad

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/stage4/ds3-collision'


def clean(v):
    if isinstance(v,mp.mpc):return [str(v.real),str(v.imag)]
    if isinstance(v,mp.mpf):return str(v)
    if isinstance(v,dict):return {k:clean(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [clean(x) for x in v]
    return v


def run():
    BASE.mkdir(parents=True,exist_ok=True);dest=BASE/'study.json'
    if dest.exists():raise RuntimeError('Keep prior run')
    data=dict(status='running',versions=dict(python=platform.python_version(),mpmath=mp.__version__),rows=[],
              parameters=dict(p='.6',q='.6',u='.3',v='-.2',masses=[2,4,8,16],offsets=['0','.1','.4']),
              error_scope='Pointwise precision/series refinement only; not a distributional limit or certified truncation bound')
    start=time.time()
    for M in data['parameters']['masses']:
        for offset in data['parameters']['offsets']:
            vals=[];details=[];runs=[]
            for dps,mult,K in [(50,8,20),(80,16,28)]:
                with mp.workdps(dps):
                    c=mp.mpf(M);b=c+mp.mpf('.3');a=c+mp.mpf('.1');d=a+c-b+mp.mpf(offset)
                    N=max(128,mult*M);t=time.time()
                    val,detail=physical_resolved(b,d,[a,mp.mpf('.6'),c,mp.mpf('.6')],N,K)
                    vals.append(val);details.append(clean(detail));runs.append(dict(dps=dps,N=N,K=K,seconds=time.time()-t,value=clean(val)))
            with mp.workdps(80):
                err=abs(vals[0]-vals[1]);rel=err/abs(vals[1])
                row=dict(M=M,offset=offset,runs=runs,absolute_difference=str(err),relative_difference=str(rel),details=details)
            data['rows'].append(row);data['seconds']=time.time()-start
            dest.write_text(json.dumps(data,indent=2),encoding='utf-8')
            print(M,offset,'U',mp.nstr(vals[-1],12),'relative',mp.nstr(rel,5),'seconds',round(data['seconds'],1),flush=True)
    paths=['code/run_ds3_collision.py','code/ds3_racah_collision.py','code/global_racah.py',
           'code/ds3_crossing_dictionary.py','code/ds3_local_vertex.py']
    data.update(status='completed',hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(data,indent=2),encoding='utf-8')


if __name__=='__main__':run()
