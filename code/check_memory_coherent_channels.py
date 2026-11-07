"""Check the pure-rank charged covariance as a complete rectangular map."""
from pathlib import Path
from functools import reduce
from math import gcd
from fractions import Fraction as F
import hashlib,json,platform,time
import numpy as np
import sympy as s
from memory_quadratic_contraction import Contraction
from memory_rank_blocks import casimir_matrix,projectors,exact_integer
from memory_coherent_cross import charged
from cubic_linear_dynamics import weights

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-coherent-cross'


def charged_projectors(o):
    d,e=o.d,o.dm;I=np.eye(d*e,dtype=np.int64)
    ops={}
    for i in range(o.N):
        for j in range(o.N):
            # Lower-sector E_ij = delta_ij - c_j c_i^dagger.
            upper=o.cs[i].conj().T@o.cs[j]
            lower=(np.eye(e) if i==j else np.zeros((e,e)))-o.cs[j]@o.cs[i].conj().T
            ops[i,j]=np.kron(upper,np.eye(e))-np.kron(np.eye(d),lower.T)
    cas=exact_integer(sum(ops[i,j]@ops[j,i] for i in range(o.N) for j in range(o.N))).real.astype(np.int64)
    roots=[o.N+2*t*(o.N-t) for t in range(min(o.m-1,o.N-o.m)+1)]
    result=[]
    for t,rt in enumerate(roots):
        num=I.copy();den=1
        for u,ru in enumerate(roots):
            if u!=t:num=num@(cas-ru*I);den*=rt-ru
        if den<0:num=-num;den=-den
        common=reduce(gcd,[den]+[abs(int(x)) for x in num.flat]);num//=common;den//=common
        assert np.array_equal(num@num,den*num)
        assert np.array_equal(cas@num,rt*num)
        result.append((num,den))
    dd=int(s.ilcm(*[den for _,den in result]))
    assert np.array_equal(sum(num*(dd//den) for num,den in result),dd*I)
    return result


def main():
    target=OUT/'charged-map-checks.json'
    if target.exists():raise FileExistsError(target)
    start=time.time();rows=[]
    for N,p,m in [(4,1,2),(5,1,2),(5,2,3),(6,2,3)]:
        o=Contraction(N,p,m)
        pp=projectors(casimir_matrix(N,m),N,min(m,N-m))
        pm=projectors(casimir_matrix(N,m-1),N,min(m-1,N-m+1))
        pt=charged_projectors(o);omega=weights(N,m,p)
        for a in range(p+1):
            numh,dh=pp[a];numk,dk=pm[a]
            h=exact_integer((o.hs.reshape(len(o.rs),-1)@numh.T).reshape(o.hs.shape))
            k=exact_integer((o.ks.reshape(len(o.rs),-1)@numk.T).reshape(o.ks.shape))
            raw=exact_integer(np.einsum('a,aij,akl->iljk',o.rs,h,k,optimize=True).reshape(o.d*o.dm,-1))
            coeff=[F(o.C)*omega[a]*F(m-p,m)*charged(N,m,a,t) for t in range(len(pt))]
            den=int(s.ilcm(*[(v/F(pd)).denominator for v,(_,pd) in zip(coeff,pt)]))
            expected=sum(int(v*den/F(pd))*pn for v,(pn,pd) in zip(coeff,pt))
            left=exact_integer(den*raw);right=exact_integer(2*dh*dk*expected)
            assert np.array_equal(left,right),(N,p,m,a)
            rows.append(dict(N=N,p=p,m=m,input_rank=a,map_entries=raw.size,eigenvalues=list(map(str,coeff))))
        print('complete charged maps',N,p,m,'passed',flush=True)
    names=['code/check_memory_coherent_channels.py','code/memory_coherent_cross.py',
           'code/memory_quadratic_contraction.py','code/memory_rank_blocks.py',
           'code/cubic_linear_dynamics.py','code/singlet_wick_moments.py','code/complex_syk_dictionary.py']
    result=dict(status='passed',rows=rows,complete_maps=len(rows),complex_entries=sum(r['map_entries'] for r in rows),
                elapsed_seconds=time.time()-start,scope='Complete rectangular covariance maps; exact CAR, neutral input and charged output Casimir projectors. No random test vectors.',
                environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__),
                input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
