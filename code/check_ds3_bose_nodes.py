"""Preserved high-precision checks at the largest weighted condition proxies."""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from ds3_crossing_dictionary import physical_scalar

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/stage4/ds3-bose'

def main():
    dest=BASE/'weighted-node-checks.json'
    if dest.exists():raise RuntimeError('Preserve prior run')
    summary=json.loads((BASE/'summary.json').read_text())
    mp.mp.dps=55;rows=[]
    for r in summary['results'][-1]['worst_weighted_nodes']:
        a,b,c,d=[mp.mpf(r[k]) for k in ['a','b','c','d']]
        ref=physical_scalar(b,d,[a,mp.mpf(.6),c,mp.mpf(.6)])
        err=abs(mp.mpc(*r['u'])-ref)
        rows.append(dict(**r,reference=[str(ref.real),str(ref.imag)],absolute_error=str(err),
                         relative_error=str(err/abs(ref)),weighted_error=str(mp.mpf(r['weight'])*err)))
        print(r['a'],r['b'],r['d'],str(err/abs(ref)),flush=True)
    paths=['code/check_ds3_bose_nodes.py','code/ds3_crossing_dictionary.py','code/global_racah.py',
           'results/stage4/ds3-bose/summary.json']
    result=dict(status='completed',dps=55,rows=rows,
                scope='Selected largest weighted proxies only; not a rigorous total quadrature error estimate',
                hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8')

if __name__=='__main__':main()
