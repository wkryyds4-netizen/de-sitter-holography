"""High precision comparisons at influential, ill-conditioned and generic points."""
import hashlib
import json
import time
from pathlib import Path
import mpmath as mp
import numpy as np
from ds3_racah_fast import physical
from ds3_crossing_dictionary import physical_scalar

ROOT=Path(__file__).resolve().parents[1]


def run():
    mp.mp.dps=45
    base=ROOT/'results/stage4/ds3-crossed-packet'
    pilot=json.loads((base/'pilot-n4-L8.json').read_text())
    samples=[dict(kind='influential',**{k:r[k] for k in ['a','b','c','d']}) for r in pilot['worst_weighted_nodes'][:6]]
    with np.load(base/'pilot-n4-L8.npz') as f:
        i,j=np.unravel_index(np.argmax(f['condition']),f['condition'].shape)
        ia,ib,id=np.unravel_index(j,(len(f['x']),)*3)
        samples.append(dict(kind='max_condition',a=float(f['x'][ia]),b=float(f['x'][ib]),d=float(f['x'][id]),c=float(f['c'][i])))
    rng=np.random.default_rng(7361)
    for _ in range(6):
        a,b,d=rng.uniform(.2,12,3);c=rng.uniform(.4,.8)
        samples.append(dict(kind='generic',a=float(a),b=float(b),d=float(d),c=float(c)))
    rows=[]
    for r in samples:
        a,b,c,d=[r[k] for k in ['a','b','c','d']]
        v,cond=physical(b,d,[a,.9,c,.6],128,12)
        v2,_=physical(b,d,[a,.9,c,.6],256,16)
        start=time.time()
        # Preserve actual binary64 input values in the high-precision reference.
        reference=physical_scalar(mp.mpf(b),mp.mpf(d),list(map(mp.mpf,[a,.9,c,.6])))
        err=abs(mp.mpc(complex(v))-reference)
        row=dict(**r,value=[float(v.real),float(v.imag)],refined=[float(v2.real),float(v2.imag)],
                 reference=[str(reference.real),str(reference.imag)],absolute_error=str(err),
                 relative_error=str(err/abs(reference)),condition=float(cond),
                 refinement_difference=float(abs(v-v2)),seconds=time.time()-start)
        rows.append(row);print(r['kind'],float(err/abs(reference)),flush=True)
        data=dict(dps=45,seed=7361,rows=rows,status='running',
                  note='Selected comparisons; generic high-condition relative errors are retained, not all required to pass.')
        (base/'accelerator-validation.json').write_text(json.dumps(data,indent=2),encoding='utf-8')
    paths=['code/ds3_racah_fast.py','code/validate_ds3_racah_fast.py',
           'code/ds3_crossing_dictionary.py','code/global_racah.py',
           'results/stage4/ds3-crossed-packet/pilot-n4-L8.json',
           'results/stage4/ds3-crossed-packet/pilot-n4-L8.npz']
    data.update(status='completed',hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})
    (base/'accelerator-validation.json').write_text(json.dumps(data,indent=2),encoding='utf-8')


if __name__=='__main__':run()
