"""Full same-J CAR coefficient Gram, independently checking coherent cross."""
from pathlib import Path
from fractions import Fraction as F
import hashlib, json, platform, time
import numpy as np
import sympy as s
from memory_quadratic_contraction import Contraction, anti
from memory_rank_blocks import casimir_matrix, projectors, exact_integer
from memory_parity_moments import product_table
from memory_coherent_cross import coherent_table, cross_BX_A
from cubic_linear_dynamics import weights

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT/'results/susskind-memory-central-cross'


def dot(a,b):
    bound = 2*max(float(abs(a).max()),float(abs(b).max()))**2*a.size
    assert bound < 2**61
    ar,ai = a.real.astype(np.int64),a.imag.astype(np.int64)
    br,bi = b.real.astype(np.int64),b.imag.astype(np.int64)
    return int(np.sum(ar*br+ai*bi,dtype=np.int64))


def main():
    OUT.mkdir(exist_ok=True)
    target = OUT/'car-exploration.json'
    if target.exists(): raise FileExistsError(target)
    start=time.time(); rows=[]
    for N,p,m in [(5,1,2),(5,2,3),(6,1,4),(6,2,4)]:
        o=Contraction(N,p,m)
        cas=casimir_matrix(N,m); Ps=projectors(cas,N,min(m,N-m))
        den=int(s.ilcm(*[d for _,d in Ps]))
        rank=sum(r*num*(den//d) for r,(num,d) in enumerate(Ps))
        rank2=sum(r*r*num*(den//d) for r,(num,d) in enumerate(Ps))
        def values(h,k):
            square=h@h
            rh=(rank@h.ravel()).reshape(h.shape)
            r2h=(rank2@h.ravel()).reshape(h.shape)
            sr=den*anti(h,rh)-den*(rank@square.ravel()).reshape(h.shape)
            xr=den*((N+1)*anti(h,rh)-anti(h,r2h))-2*rh@rh-den*den*(cas@square.ravel()).reshape(h.shape)/2
            return exact_integer(np.array([den*den*square,sr,xr,den*den*o.A(h,k)]))
        diagonal=[values(h,k) for h,k in zip(o.hs,o.ks)]
        G=[[F(0) for _ in range(4)] for _ in range(4)]
        count=len(o.rs)
        for i in range(count):
            for j in range(i,count):
                z=2*diagonal[i] if i==j else values(o.hs[i]+o.hs[j],o.ks[i]+o.ks[j])-diagonal[i]-diagonal[j]
                factor=F(int(o.rs[i]*o.rs[j])*(1 if i==j else 2),8*den**4*o.d)
                for a in range(4):
                    for b in range(a,4): G[a][b]+=factor*dot(z[a],z[b])
        for a in range(4):
            for b in range(a): G[a][b]=G[b][a]
        for a in range(4):
            for b in range(4): G[a][b]/=N**(int(a==2)+int(b==2))
        from memory_central_norm import exchanged,exact_gram
        from cubic_linear_dynamics import crossing
        omega=weights(N,m,p);prediction=F(0);xbpred=F(0)
        cf=F(m*(N-m)*(N+1),N)
        ct=F((m-1)*(N-m)*(N+1),N-1)
        g=F((m-p)**2,m);hh=2*p-m
        for a in range(p+1):
            for b in range(p+1):
                xa=a*(N+1-a);xb=b*(N+1-b)
                chi=crossing(N,m,a,b);xi=exchanged(N,m,a,b)
                un=-2*a*b+F(xa*xb,cf);uc=-2*a*b+F(xa*xb,ct)
                ec=-2*a*b*xi+m*(xi-chi)
                xbpred+=o.C**2*omega[a]*omega[b]*(un-2*a*b*chi)/N
                prediction+=o.C**2*omega[a]*omega[b]*(g*(uc+ec)+hh*(un-2*a*b*chi))/N
        oldgram=exact_gram(N,p,m)
        assert G[0][0]==o.C**2*oldgram[0][0]
        assert G[0][3]==o.C**2*oldgram[0][1]
        assert G[3][3]==o.C**2*oldgram[1][1]
        assert G[0][2]==xbpred
        print('candidate X--A difference',N,p,m,G[2][3]-prediction,flush=True)
        gm=s.Matrix(G)
        assert all(gm[:k,:k].det()>=0 for k in range(1,5))
        ids=list(gm.extract([0,3],[0,3]).rref()[1])
        old=[ [0,3][i] for i in ids]
        residual=s.factor(gm[1,1]-(gm.extract([1],old)*gm.extract(old,old).inv()*gm.extract(old,[1]))[0])
        assert residual>=0
        rows.append(dict(N=N,p=p,m=m,coordinate_pairs=count*(count+1)//2,
                         basis=['W2 H^2','S_N','X_N/N','W2 A'],
                         exact_gram=[[str(x) for x in row] for row in G],
                         candidate_XA=str(prediction),candidate_XA_difference=str(G[2][3]-prediction),
                         S_residual_after_B_A=str(residual)))
        print('CAR cross',N,p,'passed; residual',residual,flush=True)
    names=['code/explore_memory_central_cross.py','code/memory_central_norm.py','code/memory_coherent_cross.py',
           'code/memory_quadratic_contraction.py','code/memory_rank_blocks.py',
           'code/memory_parity_moments.py','code/cubic_linear_dynamics.py',
           'code/singlet_wick_moments.py','code/complex_syk_dictionary.py']
    data=dict(status='exploratory candidate check',rows=rows,elapsed_seconds=time.time()-start,
              method='All symmetric Gaussian coordinate pairs; exact integer CAR and Fraction 4x4 Gram. Independent comparison with the mixed Racah formula, including S and Casimir-deficit insertions.',
              scope='General-charge exact CAR Gram. Candidate mixed exchanged Casimir formula is exploratory; differences explicitly recorded.',
              environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__),
              input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(data,indent=2)+'\n')


if __name__=='__main__': main()
