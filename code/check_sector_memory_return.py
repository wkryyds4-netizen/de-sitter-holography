"""Fresh checks only; refuse to overwrite a successful scientific output."""
from datetime import datetime,timezone
from pathlib import Path
from collections import Counter
import hashlib,json,platform,time
import numpy as np
import scipy
import sector_memory_return as mr

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/susskind-sector-memory-return/interface-checks.json'


def formal_inverse(maximum):
    # Free words H^a F_b H^c; every coefficient is an exact integer.
    rows=[]
    def kpoly(n,left=0,right=0):
        return Counter({(a+left,b,n-2-a-b+right):1
                        for a in range(n-1) for b in range(n-1-a)})
    for n in range(maximum+1):
        out=Counter()
        for sign,poly in [(1,kpoly(n+2)),(-1,kpoly(n+1,1,0)),
                          (-1,kpoly(n+1,0,1)),(1,kpoly(n,1,1))]:
            for word,coefficient in poly.items():out[word]+=sign*coefficient
        clean={word:v for word,v in out.items() if v}
        assert clean=={(0,n,0):1}
        rows.append(dict(degree=n,expanded_words=len(out),survivors=1))
    return rows


def main():
    if OUT.exists():raise FileExistsError(OUT)
    start=time.time()
    formal=formal_inverse(20)
    finite=[]
    for args in [(3,1,2,3),(4,1,2,2),(4,2,3,1)]:
        finite.append(mr.inverse_rows(*args))
        print(json.dumps({'finite':args,'status':'passed'}),flush=True)
    dilation=[]
    times=[-1.1,-.2,0.,.4,1.1]
    for n,p,m,seed in [(3,1,2,731),(4,1,2,732),(4,2,3,733)]:
        dilation.append(mr.dilation_rows(n,p,m,seed,times))
    counts=dict(formal_degrees=len(formal),finite_rows=sum(len(x['rows']) for x in finite),
                nonzero_lower_chaos_rows=sum(y['lower_chaos_norm']>1e-10 for x in finite for y in x['rows']),
                dilation_times=sum(len(x['rows']) for x in dilation),
                schatten_bounds=2*3*sum(len(x['rows']) for x in dilation))
    assert counts['nonzero_lower_chaos_rows']>0
    inputs=['code/check_sector_memory_return.py','code/sector_memory_return.py','code/complex_syk_dictionary.py']
    result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),
                seconds=time.time()-start,counts=counts,formal=formal,finite=finite,dilation=dilation,
                parameters=dict(C_m=1,times=times,hermite='probabilists; full Gaussian chaos retained',
                                norm='normalized source-dimension Schatten; finite checks are not asymptotic proofs'),
                environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
                input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    OUT.parent.mkdir(parents=True,exist_ok=True)
    with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({'status':'passed','counts':counts,'sha256':hashlib.sha256(OUT.read_bytes()).hexdigest()}))


if __name__=='__main__':main()
