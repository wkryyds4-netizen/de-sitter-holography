"""All-coordinate CAR check of general-charge coherent variance."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,platform,time
import numpy as np
from memory_central_norm import pair,exact_gram
from memory_quadratic_contraction import Contraction
from memory_rank_blocks import exact_integer
from check_memory_coherent_car import dot

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-central-norm'


def main():
    target=OUT/'finite-norm.json'
    if target.exists():raise FileExistsError(target)
    start=time.time();identities=0
    for N in range(5,15):
        for m in range(2,N-1):
            for p in range(1,min(4,m,N-m+1)):
                for a in range(p+1):
                    for b in range(p+1):
                        pair(N,p,m,a,b);identities+=1
    rows=[]
    for N,p,m in [(5,1,2),(5,2,3),(6,1,4),(6,2,4)]:
        o=Contraction(N,p,m)
        def values(h,k):return exact_integer(np.array([h@h,o.A(h,k)]))
        diag=[values(h,k) for h,k in zip(o.hs,o.ks)]
        G=[[F(0) for _ in range(2)] for _ in range(2)]
        for i in range(len(o.rs)):
            for j in range(i,len(o.rs)):
                zz=2*diag[i] if i==j else values(o.hs[i]+o.hs[j],o.ks[i]+o.ks[j])-diag[i]-diag[j]
                factor=F(int(o.rs[i]*o.rs[j])*(1 if i==j else 2),8*o.d)
                for a in range(2):
                    for b in range(a,2):G[a][b]+=factor*dot(zz[a],zz[b])
        G[1][0]=G[0][1]
        prediction=[[o.C**2*v for v in row] for row in exact_gram(N,p,m)]
        assert G==prediction,(N,p,m,G,prediction)
        rows.append(dict(N=N,p=p,m=m,coordinate_pairs=len(o.rs)*(len(o.rs)+1)//2,
                         gram=list(map(lambda row:list(map(str,row)),G)),C=o.C))
        print('general CAR',N,p,m,'passed',flush=True)
    names=['code/check_memory_central_norm.py','code/memory_central_norm.py',
           'code/memory_quadratic_contraction.py','code/memory_rank_blocks.py',
           'code/check_memory_coherent_car.py','code/cubic_linear_dynamics.py',
           'code/singlet_wick_moments.py','code/complex_syk_dictionary.py']
    target.write_text(json.dumps(dict(status='passed',general_pair_checks=identities,rows=rows,
        exact_CAR_coordinate_pairs=sum(r['coordinate_pairs'] for r in rows),elapsed_seconds=time.time()-start,
        scope='General-charge all-coordinate same-disorder 2x2 Gram for W2H2 and W2A. Rational pair checks include contiguous relation and norm positivity. Not a proof of a general formula by sampling alone.',
        environment=dict(python=platform.python_version(),numpy=np.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')


if __name__=='__main__':main()
