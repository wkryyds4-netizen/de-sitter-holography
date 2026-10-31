"""A separate late-time scan; retains all seeds and signed estimates."""
import hashlib,json
from pathlib import Path
import numpy as np
from principal_eikonal import sobol_integral
from check_principal_eikonal import fixed_rapidity
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/stage4/principal-eikonal'

def main():
    if (OUT/'late-runs.jsonl').exists():raise RuntimeError('Already run')
    moments=json.loads((OUT/'moments.json').read_text())
    lookup={(r['n'],r['nu'],r['sigma']):float(r['value']) for r in moments}
    fixed=[]
    for lam in [-10.,-100.,-1000.]:
        values=[fixed_rapidity(lam,.6,.9,0,0,n) for n in [48,72]]
        fixed.append(dict(lam=lam,values=[[v.real,v.imag] for v in values],difference=float(abs(values[1]-values[0]))))
    rows=[];summaries=[]
    with (OUT/'late-runs.jsonl').open('x') as file:
        for u in [4.,6.,8.]:
            for rate,power,offset in [(1.,20,11100),(.7,21,13100)]:
                group=[];base=lookup[0,.6,.4]*lookup[0,.9,.4]
                for seed in range(8):
                    row=sobol_integral(.6,.9,.4,np.log(1000)+u,1000,.25,power,offset+seed,base,rate)
                    row.update(u=u,sigma=.4,gamma=.25,m=1000,nu=.6,nup=.9,base=base)
                    rows.append(row);group.append(row);file.write(json.dumps(row)+'\n');file.flush()
                a=np.array([r['raw_value'] for r in group]);b=np.array([r['value'] for r in group])
                summary=dict(u=u,rate=rate,power=power,raw_mean=a.mean(axis=0).tolist(),
                             raw_empirical_stderr=(a.std(axis=0,ddof=1)/np.sqrt(8)).tolist(),
                             cv_mean=b.mean(axis=0).tolist(),cv_empirical_stderr=(b.std(axis=0,ddof=1)/np.sqrt(8)).tolist(),
                             absolute_raw_mean=float(np.mean([r['absolute_raw'] for r in group])))
                summaries.append(summary);print(json.dumps(summary),flush=True)
    inputs=['code/principal_eikonal.py','code/check_principal_eikonal.py','code/late_principal_eikonal.py','results/stage4/principal-eikonal/moments.json']
    result=dict(status='completed; exploratory late-time sign scan',fixed=fixed,rows=rows,summaries=summaries,
                hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    (OUT/'late.json').write_text(json.dumps(result,indent=2))

if __name__=='__main__':main()
