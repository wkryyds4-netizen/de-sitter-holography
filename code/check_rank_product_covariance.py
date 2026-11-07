"""Direct CAR/projector products versus the independent Racah probabilities."""
from fractions import Fraction as F
from math import comb
import hashlib,json,platform
import numpy as np
from singlet_wick_moments import gue_basis
from memory_rank_blocks import ROOT,OUT,exact_integer,dot_integer
from rank_product_moments import product_distribution


def main():
    path=OUT/'product-covariance-checks.json'
    if path.exists():raise FileExistsError(path)
    saved=json.loads((OUT/'rank-blocks.json').read_text());N,m=5,3;D=comb(N,m)
    _,vs,rs=gue_basis(N,m);ix=np.array([a for a in range(1<<N) if a.bit_count()==m]);vs=vs[:,ix[:,None],ix]
    Ps=[(np.array(r['numerator'],dtype=np.int64),r['denominator']) for r in saved['output_projectors']]
    dims=[comb(N,j)**2-(comb(N,j-1)**2 if j else 0) for j in range(len(Ps))]
    tensors=[exact_integer((v@vs.reshape(len(vs),-1).T).T).reshape(vs.shape) for v,_ in Ps]
    weights=np.broadcast_to(rs[:,None,None]*rs[None,:,None],(len(rs),len(rs),D*D))
    rows=[]
    for a,b in [(1,2),(2,2)]:
        expected=product_distribution(N,m,a,b);da=Ps[a][1];db=Ps[b][1]
        product=exact_integer(np.einsum('aij,bjk->abik',tensors[a],tensors[b],optimize=True)).reshape(len(rs),len(rs),-1)
        reverse=exact_integer(np.einsum('bij,ajk->abik',tensors[b],tensors[a],optimize=True)).reshape(product.shape)
        probabilities=[];symmetric=[]
        for r,(pn,pd) in enumerate(Ps):
            z=exact_integer(np.einsum('abz,yz->aby',product,pn,optimize=True))
            zz=exact_integer(np.einsum('abz,yz->aby',reverse,pn,optimize=True))
            denominator=4*dims[a]*dims[b]*(da*db*pd)**2
            pr=F(D*dot_integer(z,z,weights),denominator)
            sy=F(D*dot_integer(z+zz,z+zz,weights),denominator)
            assert pr==expected[r] and 0<=sy<=4*pr
            probabilities.append(str(pr));symmetric.append(str(sy))
        rows.append(dict(a=a,b=b,probabilities=probabilities,symmetric_strength_over_CaCb=symmetric,
                         bound='symmetric_strength <= 4 * ordered_probability'))
        print('Independent rank product',a,b,probabilities,flush=True)
    names=['code/check_rank_product_covariance.py','code/rank_product_moments.py','code/memory_rank_blocks.py',
           'results/susskind-memory-rank-blocks/rank-blocks.json','code/singlet_wick_moments.py']
    path.write_text(json.dumps(dict(status='passed',N=N,m=m,rows=rows,
        independent_method='Full CAR-generated GUE coordinates projected by exact Casimir polynomials, all 10000 ordered Gaussian coordinate pairs, all output ranks. Fraction normalization.',
        environment=dict(python=platform.python_version(),numpy=np.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')


if __name__=='__main__':main()
