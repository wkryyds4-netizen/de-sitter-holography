"""Independent pure-rank CAR check of coherent Casimir insertions."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,platform,time
import numpy as np
import sympy as s
from memory_quadratic_contraction import Contraction
from memory_rank_blocks import casimir_matrix,projectors,exact_integer,dot_integer
from memory_central_norm import exchanged
from cubic_linear_dynamics import crossing,weights

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-central-cross'


def main():
    target=OUT/'charged-casimir.json'
    if target.exists():raise FileExistsError(target)
    start=time.time();rows=[]
    for N,p,m in [(4,1,2),(5,1,2),(5,2,3),(6,1,4)]:
        o=Contraction(N,p,m);count=len(o.rs)
        cp=casimir_matrix(N,m);cm=casimir_matrix(N,m-1)
        pp=projectors(cp,N,min(m,N-m));pm=projectors(cm,N,min(m-1,N-m+1))
        hs=[exact_integer((o.hs.reshape(count,-1)@num.T).reshape(o.hs.shape)) for num,den in pp[:p+1]]
        ks=[exact_integer((o.ks.reshape(count,-1)@num.T).reshape(o.ks.shape)) for num,den in pm[:p+1]]
        omega=weights(N,m,p);g=F((m-p)**2,m)
        ct=F((m-1)*(N-m)*(N+1),N-1)
        ww=np.broadcast_to(o.rs[:,None,None]*o.rs[None,:,None],(count,count,o.d*o.d))
        for a in range(p+1):
            for b in range(a,p+1):
                prod=exact_integer(np.einsum('aij,bjk->abik',hs[a],hs[b],optimize=True))
                rev=exact_integer(np.einsum('bij,ajk->abik',hs[b],hs[a],optimize=True))
                low=exact_integer(np.einsum('aij,bjk->abik',ks[a],ks[b],optimize=True))
                removed=exact_integer(sum(np.einsum('ij,abjk,kl->abil',c.conj().T,low,c,optimize=True) for c in o.cs))
                prod=prod.reshape(count,count,-1);rev=rev.reshape(count,count,-1);removed=removed.reshape(count,count,-1)
                casprod=exact_integer(prod@cp.T);casrev=exact_integer(rev@cp.T)
                den=4*o.d*pp[a][1]*pp[b][1]*pm[a][1]*pm[b][1]
                actual=[F(dot_integer(v,removed,ww),den*(2 if i>=2 else 1))
                        for i,v in enumerate([prod,rev,casprod,casrev])]
                xa=a*(N+1-a);xb=b*(N+1-b);chi=crossing(N,m,a,b);xi=exchanged(N,m,a,b)
                scale=o.C**2*omega[a]*omega[b]*g
                expected=[scale,scale*xi,scale*(xa+xb-F(xa*xb,ct)),scale*((xa+xb)*xi-m*(xi-chi))]
                assert actual==expected,(N,p,m,a,b,actual,expected)
                rows.append(dict(N=N,p=p,m=m,a=a,b=b,ordered_coordinate_pairs=count*count,
                                 actual=list(map(str,actual)),expected=list(map(str,expected))))
        print('pure rank charged Casimir',N,p,m,'passed',flush=True)
    names=['code/check_central_charged_casimir.py','code/memory_central_norm.py',
           'code/memory_quadratic_contraction.py','code/memory_rank_blocks.py',
           'code/cubic_linear_dynamics.py','code/singlet_wick_moments.py','code/complex_syk_dictionary.py']
    target.write_text(json.dumps(dict(status='passed',rows=rows,exact_identities=4*len(rows),
        ordered_coordinate_pairs=sum(r['ordered_coordinate_pairs'] for r in rows),elapsed_seconds=time.time()-start,
        scope='Complete pure-rank same-coupling rectangular contractions, direct and exchanged zeroth/first Casimir moments. No half-filled parity identity is used.',
        environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')


if __name__=='__main__':main()
