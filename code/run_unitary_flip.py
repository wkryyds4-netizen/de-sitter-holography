"""Retain all runs, including failure to reproduce the non-flat pentagon."""
from pathlib import Path
import json,hashlib,platform,time
import numpy,scipy
from unitary_flip import fft_run,hermite_run
root=Path(__file__).resolve().parents[1];base=root/'results/stage4/unitary-flip'
base.mkdir(exist_ok=True,parents=True)
if (base/'runs.jsonl').exists():raise RuntimeError('Do not overwrite runs')
rows=[]
with (base/'runs.jsonl').open('x') as f:
    for method,sizes in [('fft',[(256,16),(512,24),(1024,32)]),('hermite',[(64,None),(128,None),(256,None),(512,None)])]:
        for n,L in sizes:
            for h,flat in [(.5,False),(1.,False),(1.,True)]:
                t=time.perf_counter();r=fft_run(n,L,h,flat) if method=='fft' else hermite_run(n,h,flat)
                r['seconds']=time.perf_counter()-t;rows.append(r);f.write(json.dumps(r)+'\n');f.flush()
                print(method,n,h,flat,r['pentagon_residual'],flush=True)
inputs=['code/unitary_flip.py','code/run_unitary_flip.py']
result=dict(status='completed; non-flat pentagon reproduction unresolved',rows=rows,
    versions=dict(python=platform.python_version(),numpy=numpy.__version__,scipy=scipy.__version__),
    hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in inputs})
(base/'study.json').write_text(json.dumps(result,indent=2))
