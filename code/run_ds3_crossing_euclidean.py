"""Reproducible Euclidean quadratures, checkpointed after each grid.

The original 8,16,24,40 power-4 run is retained in euclidean.json.
This independent refinement uses a smoother radial map and puncture cuts.
"""
import hashlib
import json
import platform
import time
from pathlib import Path
import mpmath as mp
from ds3_crossing_dictionary import euclidean_overlap
from global_racah import scalar_6j
from validate_ds3_crossing_dictionary import encode

ROOT=Path(__file__).resolve().parents[1]


def run():
    mp.mp.dps=35
    out=ROOT/'results/stage4/ds3-crossing-dictionary/euclidean-refined.json'
    if out.exists(): raise RuntimeError('Preserve existing run; use a new output path.')
    cases=[(['.4','.7','.93','.6'],'.65','.57',[(24,8,'0'),(40,8,'0'),
             (32,4,'0.001'),(32,4,'0.0001')]),
           (['1.1','.43','.82','1.3'],'.37','1.17',[(24,8,'0'),(40,8,'0')])]
    paths=['code/ds3_crossing_dictionary.py','code/global_racah.py',
           'code/run_ds3_crossing_euclidean.py']
    data=dict(dps=35,python=platform.python_version(),mpmath=mp.__version__,
              hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
              rows=[],status='running',error_kind='empirical, not interval certified')
    for labels,bs,ds,grids in cases:
        e=list(map(mp.mpf,labels));b,d=mp.mpf(bs),mp.mpf(ds)
        exact=scalar_6j(b,d,e)['value']
        for n,p,lo in grids:
            start=time.time();r=euclidean_overlap(b,d,e,n,p,lo)
            r.update(external=e,b=b,d=d,analytic=exact,
                     relative_error=abs(r['value']/exact-1),seconds=time.time()-start)
            data['rows'].append(r)
            out.write_text(json.dumps(encode(data),indent=2),encoding='utf-8')
            print(labels,n,p,lo,mp.nstr(r['relative_error'],12),r['seconds'],flush=True)
    data['status']='completed'
    out.write_text(json.dumps(encode(data),indent=2),encoding='utf-8')


if __name__=='__main__':run()
