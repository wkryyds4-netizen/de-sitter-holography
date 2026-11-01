from pathlib import Path
import json,hashlib,time,platform
import numpy as np,scipy
from flip_domain import project_run
root=Path(__file__).resolve().parents[1];base=root/'results/stage4/flip-domain'
base.mkdir(parents=True,exist_ok=True)
if (base/'runs.jsonl').exists():raise RuntimeError('Do not overwrite')
cases=[(16,128,1.,False),(32,256,1.,False),(64,256,1.,False),(64,512,1.,False),
    (128,512,1.,False),(128,768,1.,False),(128,1024,1.,False),
    (64,512,.5,False),(128,768,.5,False),(32,256,1.,True),(64,512,1.,True)]
rows=[]
with (base/'runs.jsonl').open('x') as log:
    for i,(n,q,h,flat) in enumerate(cases):
        t=time.perf_counter();r,arrays=project_run(n,q,h,flat)
        path=base/f'coefficients-{i:02d}.npz'
        if path.exists():raise RuntimeError('Do not overwrite coefficients')
        np.savez_compressed(path,**arrays)
        r.update(seconds=time.perf_counter()-t,coefficients=path.relative_to(root).as_posix(),
            coefficients_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        rows.append(r);log.write(json.dumps(r)+'\n');log.flush();print(i,r['residual'],r['quadrature_dependent_margin'],flush=True)
inputs=['code/flip_domain.py','code/run_flip_domain.py','code/unitary_flip.py']
(base/'study.json').write_text(json.dumps(dict(status='completed; continuum certification pending',rows=rows,
    versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
    hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in inputs}),indent=2))
