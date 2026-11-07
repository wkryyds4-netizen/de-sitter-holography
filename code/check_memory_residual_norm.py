"""Check the finite norm against saved full CAR Grams and contiguous identities."""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import hashlib,json,platform
from memory_residual_norm import pair_norm,pair_norm_raw,exact_norm
from memory_coherent_cross import exchanged
from cubic_linear_dynamics import crossing

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-residual-norm'


def main():
    OUT.mkdir(exist_ok=True);target=OUT/'finite-norm.json'
    if target.exists():raise FileExistsError(target)
    identities=pairs=0
    for N in range(4,22,2):
        m=N//2
        for a in range(m):
            for b in range(m):
                chi=crossing(N,m,a,b)
                assert crossing(N,m-1,a,b)==chi+(exchanged(N,a,b)-chi)/(m+1)
                identities+=1
                for p in range(max(1,a,b),m):
                    assert pair_norm(N,p,a,b)==pair_norm_raw(N,p,a,b)
                    pairs+=1
    old=ROOT/'results/susskind-memory-coherent-cross/car-cross.json'
    saved=json.loads(old.read_text());rows=[]
    for row in saved['rows']:
        N,p,m=row['N'],row['p'],row['m']
        C=comb(m,p)*comb(N-m+p,p)
        predicted=C*C*exact_norm(N,p)
        assert predicted==F(row['exact_gram'][3][3])
        rows.append(dict(N=N,p=p,m=m,exact_norm=str(predicted),C=C,
                         independent_car_coordinate_pairs=row['coordinate_pairs']))
    names=['code/memory_residual_norm.py','code/check_memory_residual_norm.py',
           'code/memory_coherent_cross.py','code/cubic_linear_dynamics.py',
           'results/susskind-memory-coherent-cross/car-cross.json']
    data=dict(status='passed',lower_charge_contiguous_identities=identities,
              rational_norm_identities=pairs,rows=rows,
              method='Exact Fraction expressions compared with previously computed complete same-J CAR Grams. Previous outputs are read, not regenerated.',
              environment=dict(python=platform.python_version()),
              input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print('passed',identities,'contiguous,',pairs,'norm identities; CAR',rows)


if __name__=='__main__':main()
