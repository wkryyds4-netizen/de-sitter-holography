"""Checkpointed full-kernel tests across the transverse removable singularity."""
import hashlib,json,platform,time
from pathlib import Path
import mpmath as mp
import sympy as sp
from ds3_uniform_local import evaluate
from run_ds3_collision import clean


def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/ds3-uniform-local';base.mkdir(exist_ok=True)
    dest=base/'study.json';raw=base/'runs.jsonl'
    if dest.exists() or raw.exists():raise RuntimeError('Preserve prior run')
    mp.mp.dps=100;rows=[]
    for M in [16,64,256]:
        for r in ['-.2','-1e-12','-1e-24','0','1e-24','1e-12','.2']:
            t=time.monotonic();row=evaluate(M,'.1','-.1','.4','.4',r)
            row['seconds']=time.monotonic()-t;row=clean(row);rows.append(row)
            with raw.open('a',encoding='utf-8') as f:f.write(json.dumps(row)+'\n');f.flush()
            print(M,r,'error',row['finite_error'][:12],'uniform bound',row['uniform_finite_bound'][:12],flush=True)
    paths=['code/ds3_uniform_local.py','code/run_ds3_uniform_local.py','code/ds3_full_tail_bound.py',
           'code/ds3_shadow_duality.py','code/ds3_transverse_phase.py','code/ds3_crossing_dictionary.py']
    data=dict(status='passed' if all(r['passed'] for r in rows) else 'failed',rows=rows,dps=100,
              versions=dict(python=platform.python_version(),mpmath=mp.__version__,sympy=sp.__version__),
              scope='Uniform finite-M series remainder across zero; final sine expansion has an analytic asymptotic, not an evaluated universal constant. Excludes rounding.',
              hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(data,indent=2),encoding='utf-8');print(data['status'])


if __name__=='__main__':main()
