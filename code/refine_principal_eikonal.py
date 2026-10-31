"""Opposite momentum rotations and changed-importance independent scrambles."""
import hashlib,json
from pathlib import Path
import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss
from principal_eikonal import sobol_integral
from check_principal_eikonal import fixed_rapidity
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/stage4/principal-eikonal'

def rotated(lam,order):
    z,w=leggauss(order);r=-18+22*(z+1)/2;w*=11;x=np.exp(r)
    with mp.workdps(45):
        phase=mp.exp(-1j*mp.pi/4)
        f=np.array([complex(mp.besselk(mp.mpc(0,.6),phase*v)**2) for v in x])
        g=np.array([complex(mp.besselk(mp.mpc(0,.9),phase*v)**2) for v in x])
    return -np.sum((w*x*x*f)[:,None]*(w*x*x*g)[None,:]*np.exp(1j*lam*x[:,None]*x[None,:]))

def main():
    if (OUT/'refinement-runs.jsonl').exists():raise RuntimeError('Do not overwrite')
    checks=[];rotation=[]
    for lam in [.2,-.2,-1.]:
        values=[rotated(lam,n) for n in [180,280]]
        reference=fixed_rapidity(lam,.6,.9,0,0,56)
        row=dict(lam=lam,rotated=[[v.real,v.imag] for v in values],rapidity=[reference.real,reference.imag])
        rotation.append(row)
        for name,error,tol in [('opposite rotation',abs(values[-1]-reference),2e-9),('rotation refinement',abs(values[1]-values[0]),2e-8)]:
            checks.append(dict(name=name+str(lam),error=float(error),tolerance=tol,passed=bool(error<tol)))
        print(json.dumps(row),flush=True)
    old=json.loads((OUT/'study.json').read_text()); rows=[];summaries=[]
    with (OUT/'refinement-runs.jsonl').open('x') as file:
        for sigma,gamma,u in [(.4,0.,0.),(.4,.25,0.),(.8,.25,0.),(.8,.25,2.)]:
            previous=[r for r in old['rows'] if (r['sigma'],r['gamma'],r['u'],r['power'])==(sigma,gamma,u,20)]
            base=previous[0]['base']
            new=[]
            for seed in range(8):
                row=sobol_integral(.6,.9,sigma,np.log(1000)+u,1000,gamma,21,9100+seed,base,rate=.7)
                row.update(sigma=sigma,gamma=gamma,u=u,m=1000,nu=.6,nup=.9,base=base)
                rows.append(row);new.append(row);file.write(json.dumps(row)+'\n');file.flush()
            a=np.array([r['value'] for r in previous]);b=np.array([r['value'] for r in new])
            sem=np.sqrt((a.var(axis=0,ddof=1)+b.var(axis=0,ddof=1))/8)
            diff=np.abs(a.mean(axis=0)-b.mean(axis=0))
            summary=dict(sigma=sigma,gamma=gamma,u=u,mean=b.mean(axis=0).tolist(),
                         empirical_stderr=(b.std(axis=0,ddof=1)/np.sqrt(8)).tolist(),difference=diff.tolist(),
                         combined_stderr=sem.tolist(),raw_mean=np.mean([r['raw_value'] for r in new],axis=0).tolist(),
                         baseline_ratio_mean=float(np.mean([r['baseline_ratio'] for r in new])),
                         absolute_correction_mean=np.mean([r['absolute_correction'] for r in new],axis=0).tolist())
            summaries.append(summary);print(json.dumps(summary),flush=True)
            # Explicit empirical consistency criterion; no certified confidence interval.
            for k in range(2):checks.append(dict(name='changed importance '+str((sigma,gamma,u,k)),error=float(diff[k]),tolerance=float(6*sem[k]+1e-6),passed=bool(diff[k]<6*sem[k]+1e-6)))
    inputs=['code/principal_eikonal.py','code/check_principal_eikonal.py','code/refine_principal_eikonal.py','results/stage4/principal-eikonal/study.json']
    out=dict(status='passed' if all(c['passed'] for c in checks) else 'failed',checks=checks,rotations=rotation,rows=rows,summaries=summaries,
             scope='Complex original momentum controls and changed importance/seeds; empirical consistency, not a total error certificate.',
             hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    (OUT/'refinement.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(dict(status=out['status'],checks=len(checks),failures=[c for c in checks if not c['passed']])))

if __name__=='__main__':main()
