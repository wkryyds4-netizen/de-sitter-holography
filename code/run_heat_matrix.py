from pathlib import Path
import json,hashlib,platform,time
import numpy as np
import scipy
from heat_matrix import fft_family,hermite_family

def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/heat-matrix'
    base.mkdir(exist_ok=True)
    if (base/'matrix-runs.jsonl').exists() or (base/'matrix-study.json').exists():raise RuntimeError('Do not overwrite')
    start=time.time();rows=[]
    cases=[('fft',n,L,False) for n,L in [(512,20),(1024,20),(2048,20),(2048,28)]]
    cases+=[('hermite',n,q,False) for n,q in [(64,256),(128,512),(256,768)]]
    cases+=[('fft',1024,20,True),('hermite',64,512,True)]
    with (base/'matrix-runs.jsonl').open('x') as raw:
        for method,n,param,flat in cases:
            eps=[0,.01,.0001,1e-8,1e-10] if method=='fft' else [0,.01,1e-10]
            run=fft_family(n,param,1,eps,flat) if method=='fft' else hermite_family(n,param,1,eps,flat)
            for row in run:rows.append(row);raw.write(json.dumps(row)+'\n')
            raw.flush();print(method,n,param,flat,'completed',flush=True)
    data=dict(status='unresolved continuum identity; quadrature not certified',rows=rows,seconds=time.time()-start,
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        scope='Actual discrete matrix elements of D_epsilon, not certified continuous values. Heat error is separately bounded analytically.',
        hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in
            ['code/heat_matrix.py','code/run_heat_matrix.py','code/unitary_flip.py']})
    (base/'matrix-study.json').write_text(json.dumps(data,indent=2));print(len(rows),'records',flush=True)
if __name__=='__main__':main()
