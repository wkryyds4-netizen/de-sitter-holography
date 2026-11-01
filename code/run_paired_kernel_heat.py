from pathlib import Path
import json,hashlib,platform,time
import numpy as np
import scipy
from paired_kernel_heat import normal_nodes,phase_loss,heat_loss_from_mean,vacuum_heat_error

def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/paired-kernel'
    base.mkdir(exist_ok=True)
    if (base/'heat-runs.jsonl').exists() or (base/'heat-study.json').exists():raise RuntimeError('Do not overwrite')
    rows=[];start=time.time()
    with (base/'heat-runs.jsonl').open('x') as raw:
        for seed in [1729,2718,3141,5772]:
            normals=normal_nodes(16,seed)
            for h in [None,.5,1.,2.]:
                for delta in [.01,.001,.0001,.00001,.000001]:
                    t=2*delta;loss,tails=phase_loss(normals,t,h)
                    for power in [12,14,16]:
                        v=loss[:2**power];mean=float(np.mean(v));heat=heat_loss_from_mean(t,mean)
                        row=dict(seed=seed,h=h,delta=delta,t=t,power=power,n=len(v),
                            mean_phase_loss=mean,mean_square=float(np.mean(v*v)),maximum=float(np.max(v)),
                            count_above_point1=int(np.sum(v>.1)),estimated_one_minus_C=float(heat),
                            estimated_weak_error_functional=float(2*vacuum_heat_error(delta)+np.sqrt(heat)),
                            product_tails=tails)
                        if h is None:row['exact_one_minus_C']=float(2*np.pi*t/(1+2*np.pi*t))
                        rows.append(row);raw.write(json.dumps(row)+'\n');raw.flush()
            print('completed seed',seed,flush=True)
    hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in
            ['code/paired_kernel_heat.py','code/run_paired_kernel_heat.py','code/unitary_flip.py']}
    data=dict(status='empirical diagnostic completed',rows=rows,seconds=time.time()-start,
        state='normalized oscillator Gaussian; v=exp(i*x*y/pi)*F_h(x,y)*f',
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        hashes=hashes,scope='QMC estimates of a true heat-tail upper-bound functional. Samples and replicate spread do not certify an upper bound. No pentagon word evaluated.')
    (base/'heat-study.json').write_text(json.dumps(data,indent=2));print(len(rows),'records',flush=True)
if __name__=='__main__':main()
