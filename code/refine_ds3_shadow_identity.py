"""Fresh high-precision checks for two retained failed old-reference comparisons."""
import hashlib,json
from pathlib import Path
import mpmath as mp
from ds3_shadow_duality import reflection
from global_racah import omega_hyper
from run_ds3_collision import clean


def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/ds3-shadow-duality';raw=base/'refinement.jsonl';dest=base/'refinement.json'
    if raw.exists() or dest.exists():raise RuntimeError('Preserve prior refinements')
    rows=[]
    for rstr in ['.1','.4']:
        for dps in [100,130]:
            mp.mp.dps=dps;M=mp.mpf(16);a=M+mp.mpf('.1');b=M+mp.mpf('.3');d=M-mp.mpf('.2')+mp.mpf(rstr)
            hs=[(1+1j*v)/2 for v in [a,mp.mpf('.6'),M,mp.mpf('.6')]];H=(1+1j*b)/2;hp=(1+1j*d)/2
            direct,condition=omega_hyper(H,1-hp,hs)
            other=omega_hyper(1-H,hp,hs)[0];pred=reflection(H,1-hp,hs)*mp.conj(other)
            error=abs(pred-direct)/abs(direct)
            row=clean(dict(M=16,offset=rstr,dps=dps,direct=direct,reflected=pred,relative_error=error,original_pair_condition=condition,tolerance='1e-75',passed=error<mp.mpf('1e-75')))
            with raw.open('a',encoding='utf-8') as f:f.write(json.dumps(row)+'\n');f.flush()
            rows.append(row);print(rstr,dps,mp.nstr(error,5),row['passed'],flush=True)
    paths=['code/refine_ds3_shadow_identity.py','code/ds3_shadow_duality.py','code/global_racah.py','results/stage4/ds3-shadow-duality/validation.json']
    data=dict(status='passed' if all(r['passed'] for r in rows) else 'failed',rows=rows,hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(data,indent=2),encoding='utf-8')


if __name__=='__main__':main()
