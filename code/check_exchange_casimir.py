"""Complete independent CAR sums for the exchanged-product Casimir moment."""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import hashlib,json,platform,time
import numpy as np
import sympy as s
from singlet_wick_moments import gue_basis
from memory_rank_blocks import casimir_matrix,projectors,exact_integer,dot_integer
from cubic_linear_dynamics import crossing

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-residual-norm'


def main():
    target=OUT/'exchange-casimir.json'
    if target.exists():raise FileExistsError(target)
    start=time.time();rows=[]
    for N,m in [(4,1),(4,2),(5,2)]:
        _,vs,rs=gue_basis(N,m)
        ix=np.array([a for a in range(1<<N) if a.bit_count()==m]);D=len(ix)
        vs=vs[:,ix[:,None],ix];count=len(rs);cas=casimir_matrix(N,m)
        Ps=projectors(cas,N,min(m,N-m))
        dims=[comb(N,r)**2-(comb(N,r-1)**2 if r else 0) for r in range(len(Ps))]
        projected=[exact_integer(vs.reshape(count,-1)@num.T).reshape(vs.shape) for num,_ in Ps]
        ww=np.broadcast_to(rs[:,None,None]*rs[None,:,None],(count,count,D*D))
        for a in range(len(Ps)):
            for b in range(a,len(Ps)):
                prod=exact_integer(np.einsum('aij,bjk->abik',projected[a],projected[b],optimize=True)).reshape(count,count,-1)
                rev=exact_integer(np.einsum('bij,ajk->abik',projected[b],projected[a],optimize=True)).reshape(prod.shape)
                twice_x=exact_integer(prod@cas.T)
                den=4*dims[a]*dims[b]*(Ps[a][1]*Ps[b][1])**2
                total=F(D*dot_integer(rev,prod,ww),den)
                swapped=F(D*dot_integer(rev,twice_x,ww),2*den)
                ordered=F(D*dot_integer(prod,twice_x,ww),2*den)
                chi=crossing(N,m,a,b);xa=a*(N+1-a);xb=b*(N+1-b)
                cf=F(m*(N-m)*(N+1),N)
                assert total==chi
                assert swapped==(xa+xb)*chi,(N,m,a,b,swapped,chi)
                assert ordered==xa+xb-F(xa*xb,1)/cf
                rows.append(dict(N=N,m=m,a=a,b=b,ordered_gaussian_pairs=count*count,
                                 exchange=str(total),exchange_x=str(swapped),ordered_x=str(ordered)))
        print('exchange Casimir',N,m,'passed',flush=True)
    names=['code/check_exchange_casimir.py','code/memory_rank_blocks.py','code/cubic_linear_dynamics.py',
           'code/singlet_wick_moments.py','code/complex_syk_dictionary.py']
    data=dict(status='passed',rows=rows,identities=3*len(rows),
              ordered_coordinate_pairs=sum(r['ordered_gaussian_pairs'] for r in rows),elapsed_seconds=time.time()-start,
              scope='Exact general-charge products; complete Gaussian coordinate sums, no sampled input tensors.',
              environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__),
              input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(data,indent=2)+'\n')


if __name__=='__main__':main()
