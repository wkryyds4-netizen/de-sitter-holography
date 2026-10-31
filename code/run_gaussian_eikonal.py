"""Save fresh smeared-eikonal quadratures without overwriting evidence."""
import hashlib,json,platform,time
from pathlib import Path
import numpy as np
import scipy
from gaussian_eikonal import smeared,enhancement,moment

def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/gaussian-eikonal'
    base.mkdir(parents=True,exist_ok=True);dest=base/'study.json';raw=base/'runs.jsonl'
    if dest.exists() or raw.exists():raise RuntimeError('Preserve prior evidence')
    rows=[];start=time.time()
    cases=[(g,s,u) for g in [0,.25] for s in [.4,.8] for u in [0,2]]
    for gamma,sigma,u in cases:
        t=np.log(1000)+u
        coarse=smeared(1000,gamma,t,sigma,20,24)
        fine=smeared(1000,gamma,t,sigma,32,40)
        boost=smeared(1000,gamma,u,sigma,32,40,boosted=True)
        naive=smeared(1000,gamma,t,sigma,32,40,center_only=True)
        row=dict(m=1000,gamma=gamma,sigma=sigma,u=u,t=t,coarse=coarse,fine=fine,boost=boost,
                 center_only=naive,difference=abs(complex(*coarse['value'])-complex(*fine['value'])))
        rows.append(row)
        with raw.open('a',encoding='utf-8') as f:f.write(json.dumps(row)+'\n');f.flush()
        print(gamma,sigma,u,fine['value'],'delta',row['difference'],flush=True)
    controls=[]
    for t in [0,-np.log(1000)]:controls.append(dict(t=t,result=smeared(1000,.25,t,.8,32,40)))
    moments=[dict(sigma=s,A=[moment(r,s) for r in [1,2,3]],Q=[enhancement(n,s) for n in [1,2]]) for s in [.1,.2,.4,.8]]
    paths=['code/gaussian_eikonal.py','code/run_gaussian_eikonal.py','code/conformal_shock.py',
           'sources/chen-stanford-tang-yang-2607.14042v1.html']
    result=dict(status='computed; convergence diagnostics follow',rows=rows,controls=controls,moments=moments,
        seconds=time.time()-start,dependencies=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        scope='Double precision quadrature of a specified four-time smeared eikonal model. Not a total error or full gravity result.',
        hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8')

if __name__=='__main__':main()
