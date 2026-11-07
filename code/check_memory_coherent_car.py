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
OUT = ROOT/'results/susskind-memory-coherent-cross'


def dot(a,b):
    bound = 2*max(float(abs(a).max()),float(abs(b).max()))**2*a.size
    assert bound < 2**61
    ar,ai = a.real.astype(np.int64),a.imag.astype(np.int64)
    br,bi = b.real.astype(np.int64),b.imag.astype(np.int64)
    return int(np.sum(ar*br+ai*bi,dtype=np.int64))


def main():
    target = OUT/'car-cross.json'
    if target.exists(): raise FileExistsError(target)
    start=time.time(); rows=[]
    for N,p,m in [(4,1,2),(6,1,3),(6,2,3)]:
        o=Contraction(N,p,m)
        cas=casimir_matrix(N,m); Ps=projectors(cas,N,m)
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
        _,_,pi=product_table(N); *_,mu=coherent_table(N,m)
        omega=weights(N,m,p); predicted=[F(0) for _ in range(3)]
        for a in range(p+1):
            for b in range(p+1):
                bx=cross_BX_A(N,p,a,b)
                predicted[0]+=o.C**2*omega[a]*omega[b]*bx[0]
                predicted[2]+=o.C**2*omega[a]*omega[b]*bx[1]
                gain=F((m-p)**2,m)
                direct=F(0)
                for r in range(m+1):
                    parity=0 if (a+b-r)%2 else 2
                    mur=mu[a,b][r] if r<m else F(0)
                    direct+=parity*(a+b-r)*(gain*mur+(2*p-m)*pi[a,b][r])
                predicted[1]+=o.C**2*omega[a]*omega[b]*direct
        assert [G[j][3] for j in range(3)]==predicted,(N,p,G,predicted)
        gm=s.Matrix(G)
        assert all(gm[:k,:k].det()>=0 for k in range(1,5))
        ids=list(gm.extract([0,3],[0,3]).rref()[1])
        old=[ [0,3][i] for i in ids]
        residual=s.factor(gm[1,1]-(gm.extract([1],old)*gm.extract(old,old).inv()*gm.extract(old,[1]))[0])
        assert residual>=0
        rows.append(dict(N=N,p=p,m=m,coordinate_pairs=count*(count+1)//2,
                         basis=['W2 H^2','S_N','X_N/N','W2 A'],
                         exact_gram=[[str(x) for x in row] for row in G],
                         predicted_cross=list(map(str,predicted)),
                         S_residual_after_B_A=str(residual)))
        print('CAR cross',N,p,'passed; residual',residual,flush=True)
    names=['code/check_memory_coherent_car.py','code/memory_coherent_cross.py',
           'code/memory_quadratic_contraction.py','code/memory_rank_blocks.py',
           'code/memory_parity_moments.py','code/cubic_linear_dynamics.py',
           'code/singlet_wick_moments.py','code/complex_syk_dictionary.py']
    data=dict(status='passed',rows=rows,elapsed_seconds=time.time()-start,
              method='All symmetric Gaussian coordinate pairs; exact integer CAR and Fraction 4x4 Gram. Independent comparison with the mixed Racah formula, including S and Casimir-deficit insertions.',
              scope='Finite half-filled systems; no DS conclusion from small systems alone.',
              environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__),
              input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(data,indent=2)+'\n')


if __name__=='__main__': main()
