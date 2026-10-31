"""Finish JSON serialization from saved quadratures; do not recompute them."""
import hashlib,json,platform
from pathlib import Path
import mpmath as mp
import numpy as np
from check_ds3_fock_vacuum import finite_fock
from ds3_fock_vacuum import raw_rate,block_bound
from ds3_two_matter import tilted
from run_ds3_collision import clean

if __name__=='__main__':
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/ds3-fock-vacuum';dest=base/'diagnostics-recovered.json'
    if dest.exists():raise RuntimeError('Preserve recovery')
    rows=[json.loads(x) for x in (base/'vacuum-runs.jsonl').read_text().splitlines()]
    assert len(rows)==20
    mp.mp.dps=60;rates=[]
    for a,b,p in [('.4','.7','.6'),('1.1','.3','1.2'),('3','4','.6'),('16.2','16','.8'),('64.4','64','1.2')]:
        a,b,p=map(mp.mpf,[a,b,p]);x=raw_rate(a,b,p);y=tilted(a,b,p);rel=abs(x-y)/y
        rates.append(clean(dict(a=a,b=b,p=p,raw=x,rate=y,relative_difference=rel,passed=rel<mp.mpf('1e-48'))))
    finite=finite_fock();passed=all(r['passed'] for r in rates+rows+finite['bounds']+finite['fourpoint']+finite['commutators'])
    paths=['code/ds3_fock_vacuum.py','code/check_ds3_fock_vacuum.py','code/recover_ds3_fock_vacuum.py',
        'code/ds3_two_matter.py','code/ds3_group_pairing.py','results/stage4/ds3-fock-vacuum/vacuum-runs.jsonl']
    data=dict(status='passed' if passed else 'failed',rates=rates,vacuum=rows,finite_fock=finite,
        working_dps=60,block_bound=clean(block_bound('.2')),
        dependencies=dict(python=platform.python_version(),mpmath=mp.__version__,numpy=np.__version__),
        initial_failure='Driver completed and saved all 20 quadrature pairs, then json.dumps rejected numpy.bool_ in finite_fock.bounds. No diagnostics.json was written. Original script and raw output retained.',
        recovery='Reuse all saved quadratures byte-for-byte; recompute five inexpensive rate identities and auxiliary finite matrices; encode numpy scalar via item(). No change to criterion.',
        scope='Vacuum return and finite Fock-path diagnostics; no fresh U evaluations or all-spin numerical certification. Quadrature is not interval certified.',
        hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths})
    def scalar(x):
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x).__name__)
    dest.write_text(json.dumps(data,indent=2,default=scalar),encoding='utf-8')
    print(data['status'],len(rates)+len(rows)+len(finite['bounds'])+len(finite['fourpoint'])+len(finite['commutators']))
    if not passed:raise SystemExit(1)
