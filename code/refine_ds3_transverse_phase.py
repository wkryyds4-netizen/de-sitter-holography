"""New 120-digit, doubled-cutoff checks, preserving the original 90-digit data."""
import hashlib,json
from pathlib import Path
import mpmath as mp
from ds3_transverse_phase import evaluate
from run_ds3_collision import clean


def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/ds3-transverse-phase'
    dest=base/'refinement.json'
    if dest.exists():raise RuntimeError('Preserve prior refinement')
    study=json.loads((base/'study.json').read_text());mp.mp.dps=120;rows=[]
    for index in [0,1,14,15]:
        old=study['rows'][index];new=evaluate(*[old[k] for k in ['M','u','v','p','q','r']],refinement=4)
        difference=abs(new['value']-mp.mpc(*old['value']))
        budget=new['tail']+mp.mpf(old['tail'])+mp.mpf('1e-75')
        rows.append(clean(dict(index=index,difference=difference,relative=difference/abs(new['value']),
                               comparison_budget=budget,passed=difference<=budget,result=new)))
    paths=['code/refine_ds3_transverse_phase.py','code/ds3_transverse_phase.py','results/stage4/ds3-transverse-phase/study.json']
    data=dict(status='passed' if all(r['passed'] for r in rows) else 'failed',rows=rows,dps=120,
              scope='Four selected points, double cutoff; absolute 1e-75 rounding allowance is empirical, not interval certified.',
              hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(data,indent=2),encoding='utf-8');print(data['status'],[(r['index'],r['relative']) for r in rows])


if __name__=='__main__':main()
