"""Checkpointed reflected full-kernel sums, with their exact-arithmetic tails."""
import hashlib,json,platform,time
from pathlib import Path
import mpmath as mp
import sympy as sp
from ds3_shadow_duality import physical
from run_ds3_collision import clean


def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/ds3-shadow-duality';base.mkdir(exist_ok=True)
    raw=base/'runs.jsonl';target=base/'study.json'
    if raw.exists() or target.exists():raise RuntimeError('Preserve prior study')
    old=json.loads((root/'results/stage4/ds3-collision/validation.json').read_text())
    refs={r['M']:r['reference'] for r in old['heavy']};rows=[]
    for M in [2,8,16,32,64,128]:
        for refinement,dps in [(1,70),(2,70),(2,110)]:
            mp.mp.dps=dps;t0=time.monotonic()
            value,tail,details=physical(str(sp.Rational(M)+sp.Rational(3,10)),str(sp.Rational(M)-sp.Rational(1,5)),[str(sp.Rational(M)+sp.Rational(1,10)),'3/5',str(M),'3/5'],refinement)
            ref=mp.mpc(*refs[M]) if M in refs else None
            error=abs(value-ref) if ref is not None else None
            row=clean(dict(M=M,dps=dps,refinement=refinement,value=value,tail=tail,reference=ref,reference_error=error,details=details,seconds=time.monotonic()-t0))
            with raw.open('a',encoding='utf-8') as f:f.write(json.dumps(row)+'\n');f.flush()
            rows.append(row);print(M,refinement,dps,'U',mp.nstr(value,12),'tail',mp.nstr(tail,5),'error',mp.nstr(error,5),flush=True)
    paths=['code/ds3_shadow_duality.py','code/run_ds3_shadow_duality.py','code/ds3_full_tail_bound.py','code/global_racah.py','code/ds3_crossing_dictionary.py','results/stage4/ds3-collision/validation.json']
    data=dict(status='completed',rows=rows,versions=dict(python=platform.python_version(),mpmath=mp.__version__,sympy=sp.__version__),
              scope='Reflected exact-arithmetic series-tail bounds; no interval rounding or full physical mass integral.',
              hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths})
    target.write_text(json.dumps(data,indent=2),encoding='utf-8')


if __name__=='__main__':main()
