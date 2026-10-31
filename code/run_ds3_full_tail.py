"""Checkpointed direct-sum full-U study; never replace saved scientific output."""
import hashlib
import json
import platform
import time
from pathlib import Path
import mpmath as mp
import sympy as sp
from ds3_full_tail_bound import physical
from run_ds3_collision import clean


def main():
    root=Path(__file__).resolve().parents[1]
    dest=root/'results/stage4/ds3-full-tail'
    dest.mkdir(exist_ok=True)
    raw=dest/'runs.jsonl'
    if raw.exists():
        raise RuntimeError('Preserve previous checkpoint stream; do not restart blindly')
    refs=json.loads((root/'results/stage4/ds3-collision/validation.json').read_text())
    references={r['M']:r['reference'] for r in refs['heavy']}
    study=json.loads((root/'results/stage4/ds3-collision/study.json').read_text())
    references[2]=next(r['runs'][-1]['value'] for r in study['rows'] if r['M']==2 and r['offset']=='0')
    rows=[]
    for M in [2,8,32]:
        for refinement,dps in [(1,110),(2,110),(2,150)]:
            mp.mp.dps=dps;t0=time.monotonic()
            b=str(sp.Rational(M)+sp.Rational(3,10));d=str(sp.Rational(M)-sp.Rational(1,5))
            external=[str(sp.Rational(M)+sp.Rational(1,10)),'3/5',str(M),'3/5']
            value,bound,details=physical(b,d,external,refinement)
            reference=mp.mpc(*references[M]);error=abs(value-reference)
            row=clean(dict(M=M,refinement=refinement,dps=dps,value=value,tail=bound,
                           reference=reference,reference_error=error,inside_tail=error<=bound,
                           seconds=time.monotonic()-t0,details=details))
            with raw.open('a',encoding='utf-8') as f:
                f.write(json.dumps(row)+'\n');f.flush()
            rows.append(row)
            print(M,refinement,dps,'error',mp.nstr(error,6),'tail',mp.nstr(bound,6),flush=True)
    paths=['code/ds3_full_tail_bound.py','code/run_ds3_full_tail.py','code/global_racah.py',
           'code/ds3_crossing_dictionary.py','code/ds3_local_vertex.py',
           'results/stage4/ds3-collision/validation.json','results/stage4/ds3-collision/study.json']
    data=dict(status='completed',rows=rows,versions=dict(python=platform.python_version(),mpmath=mp.__version__,sympy=sp.__version__),
              scope='Exact rational ratio certificates; evaluated exact-arithmetic tail bounds exclude rounding and reference error.',
              hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths})
    (dest/'study.json').write_text(json.dumps(data,indent=2),encoding='utf-8')


if __name__=='__main__':main()
