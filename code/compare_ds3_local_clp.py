"""Supplementary local normalization comparison; does not rerun mode diagnostics."""
from pathlib import Path
from datetime import datetime,timezone
from itertools import product
import hashlib,json,platform
import mpmath as mp
from ds3_local_vertex import reduced_coefficient
from recoil_poles import direct_c_pair
ROOT=Path(__file__).resolve().parents[1]

def predicted_ratio(nus,cc):
    count=sum(cc);factor=mp.mpf(1)
    if count in [1,2]:
        selected=cc.index(True) if count==1 else cc.index(False)
        factor=mp.cosh(mp.pi*nus[selected])**2
    return factor/(64*mp.pi**4)*mp.fprod(n*mp.exp(-mp.pi*n)/mp.sinh(mp.pi*n) for n in nus)

def run():
    mp.mp.dps=60;rows=[]
    for ss in [[".7","1.1",".6"],["1.4",".9","1.2"],["3.2","2.7",".5"]]:
        nus=list(map(mp.mpf,ss))
        for cc in product([False,True],repeat=3):
            coeff=reduced_coefficient(nus,cc)
            clp=direct_c_pair(nus)
            ratio=abs(coeff)**2/clp
            expected=predicted_ratio(nus,cc)
            err=abs(ratio/expected-1)
            rows.append(dict(nus=ss,conjugates=cc,coefficient=mp.nstr(coeff,65),
                norm=mp.nstr(abs(coeff)**2,65),direct_C_pair=mp.nstr(clp,65),
                ratio=mp.nstr(ratio,65),expected=mp.nstr(expected,65),
                relative_error=mp.nstr(err,65),tolerance="1e-50",passed=bool(err<mp.mpf("1e-50"))))
    files=["code/compare_ds3_local_clp.py","code/ds3_local_vertex.py","code/recoil_poles.py",
           "notes/stage4-recoil-poles.md","notes/stage4-recoil-spectral.md"]
    out=dict(status="passed" if all(r["passed"] for r in rows) else "failed",
        created_at_utc=datetime.now(timezone.utc).isoformat(),python=platform.python_version(),
        mpmath=mp.__version__,dps=60,rows=rows,
        scope="Mass-dependent paired C0 differs from the specified local coefficient norm by explicit leg factors; not a physical two-point n or finite-kappa dictionary.",
        hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files})
    dest=ROOT/"results/stage4/ds3-local-vertex/clp-comparison.json"
    assert not dest.exists(),"Do not overwrite evidence."
    dest.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(dict(status=out["status"],checks=len(rows))))
if __name__=="__main__":run()
