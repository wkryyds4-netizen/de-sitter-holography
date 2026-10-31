"""One-shot recorded principal eikonal experiment; refuses existing outputs."""
import hashlib,json,platform,time
from pathlib import Path
import numpy as np, scipy, mpmath as mp
from principal_eikonal import spectral_moment,sobol_integral

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/stage4/principal-eikonal'

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    if (OUT/'runs.jsonl').exists():raise RuntimeError('Do not overwrite runs')
    start=time.time(); moments=[]
    for sigma in [.4,.8]:
        for nu in [.6,.9]:
            for n in range(4):
                value=spectral_moment(n,nu,sigma,50)
                moments.append(dict(n=n,nu=nu,sigma=sigma,value=mp.nstr(value,48)))
    (OUT/'moments.json').write_text(json.dumps(moments,indent=2))
    lookup={(x['n'],x['nu'],x['sigma']):float(x['value']) for x in moments}
    rows=[]
    with (OUT/'runs.jsonl').open('x') as file:
        for sigma,gamma,u in [(.4,0.,0.),(.4,.25,0.),(.8,.25,0.),(.8,.25,2.)]:
            base=lookup[0,.6,sigma]*lookup[0,.9,sigma]
            for power in [18,20]:
                for seed in range(8):
                    row=sobol_integral(.6,.9,sigma,np.log(1000)+u,1000,gamma,power,7100+seed,base)
                    row.update(nu=.6,nup=.9,sigma=sigma,gamma=gamma,m=1000,u=u,base=base)
                    rows.append(row);file.write(json.dumps(row)+'\n');file.flush()
                selected=rows[-8:]; a=np.array([r['value'] for r in selected]);
                print(json.dumps(dict(sigma=sigma,gamma=gamma,u=u,power=power,
                      mean=a.mean(axis=0).tolist(),empirical_stderr=(a.std(axis=0,ddof=1)/np.sqrt(8)).tolist())),flush=True)
    inputs=['code/principal_eikonal.py','code/run_principal_eikonal.py','sources/chen-stanford-tang-yang-2607.14042v1.html']
    result=dict(status='completed; empirical quadrature errors only',elapsed_seconds=time.time()-start,
                versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__),
                moment_dps=50,rows=rows,moments=moments,
                hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    (OUT/'study.json').write_text(json.dumps(result,indent=2))

if __name__=='__main__':main()
