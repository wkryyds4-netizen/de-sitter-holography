"""Checkpointed transverse full-kernel tests; never overwrites a prior study."""
import hashlib,json,platform,time
from pathlib import Path
import mpmath as mp
import sympy as sp
from ds3_transverse_phase import evaluate
from run_ds3_collision import clean


def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/ds3-transverse-phase';base.mkdir(exist_ok=True)
    target=base/'study.json';raw=base/'runs.jsonl'
    if target.exists() or raw.exists():raise RuntimeError('Preserve prior run')
    mp.mp.dps=90;rows=[]
    for mode in ['fixed','scaled']:
        for M in [16,64,256,1024]:
            for sign in [-1,1]:
                t0=time.monotonic()
                r=sign*(mp.mpf('.08') if mode=='fixed' else 1/mp.log(M))
                row=evaluate(M,'.1','-.1','.4','.4',str(r))
                row.update(mode=mode,lambda_value=sign if mode=='scaled' else None,
                           scaled_value=mp.pi*row['value']/mp.log(M),
                           sinc_target=mp.sin(1) if mode=='scaled' else None,seconds=time.monotonic()-t0)
                row=clean(row);rows.append(row)
                with raw.open('a',encoding='utf-8') as f:f.write(json.dumps(row)+'\n');f.flush()
                print(mode,M,sign,'finite error',row['finite_error'][:12],'bound',row['finite_bound'][:12],flush=True)
    paths=['code/ds3_transverse_phase.py','code/run_ds3_transverse_phase.py','code/ds3_shadow_duality.py',
           'code/ds3_full_tail_bound.py','code/ds3_controlled_omega.py','code/ds3_peak_normalization.py']
    data=dict(status='passed' if all(r['finite_passed'] for r in rows) else 'failed',rows=rows,dps=90,
              versions=dict(python=platform.python_version(),mpmath=mp.__version__,sympy=sp.__version__),
              scope='Punctured transfer window; explicit finite-M bound and pointwise fixed-lambda scaling. No distributional limit or rounding certification.',
              hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths})
    target.write_text(json.dumps(data,indent=2),encoding='utf-8');print(data['status'],len(rows))


if __name__=='__main__':main()
