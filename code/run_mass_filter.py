"""Preserve paired full-conditional simulations; not a mass-integrated trace."""
import json,hashlib,platform
from pathlib import Path
import numpy as np, scipy, mpmath as mp
from principal_eikonal import spectral_moment
from mass_filter import paired_integral
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/stage4/mass-filter'

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    if (OUT/'runs.jsonl').exists():raise RuntimeError('Do not overwrite')
    rows=[];summaries=[];moments=[]
    with (OUT/'runs.jsonl').open('x') as file:
        for delta,u in [(1/np.sqrt(1000),0.),(.4,0.),(.8,0.),(.8,2.)]:
            b=[spectral_moment(0,v,delta,50) for v in [.6,.9]];base=float(b[0]*b[1])
            moments.append(dict(delta=delta,B0=[mp.nstr(v,48) for v in b]))
            for power in [19,21]:
                group=[]
                for seed in range(8):
                    row=paired_integral(delta,u,power,15100+seed,base)
                    rows.append(row);group.append(row);file.write(json.dumps(row)+'\n');file.flush()
                s=dict(delta=delta,u=u,power=power)
                for key in ['correlated','independent','difference']:
                    x=np.array([r[key] for r in group]);s[key]=x.mean(axis=0).tolist();s[key+'_se']=(x.std(axis=0,ddof=1)/np.sqrt(8)).tolist()
                summaries.append(s);print(json.dumps(s),flush=True)
    inputs=['code/mass_filter.py','code/run_mass_filter.py','code/principal_eikonal.py']
    result=dict(status='completed; conditional eikonal only',rows=rows,summaries=summaries,moments=moments,
                versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__),
                hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    (OUT/'study.json').write_text(json.dumps(result,indent=2))
if __name__=='__main__':main()
