"""Independent finite CAR Wick contraction plus exact covariance identities."""
from fractions import Fraction as F
from math import prod
from pathlib import Path
import json,hashlib,platform
import numpy as np
from cubic_linear_dynamics import weights,crossing,finite,b
from singlet_wick_moments import gue_basis
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-cubic-dynamics'

def main():
    OUT.mkdir(parents=True,exist_ok=True);path=OUT/'validation.json'
    if path.exists():raise FileExistsError(path)
    covariance=0;factorial=0;symmetry=0
    for N in range(3,15):
        for m in range(1,N):
            R=min(m,N-m);g=[]
            for r in range(R+1):
                w=weights(N,m,r);C=b(m,r)*b(N-m+r,r)
                g.append([(F(b(m-n,r)*b(N-m+r-n,r),C)-sum(w[a]*g[a][n] for a in range(r)))/w[r] for n in range(R+1)])
                for n in range(R+1):
                    assert g[r][n]==crossing(N,m,r,n);covariance+=1
                    assert abs(g[r][n])<=1;symmetry+=1
                for k in range(r+1):
                    lhs=sum(w[a]*prod((a-j)*(N+1-a-j) for j in range(k)) for a in range(r+1))
                    rhs=F(prod((r-j)*(N-m-j)*(N+1-j) for j in range(k)),prod(N-m+r-j for j in range(k)))
                    assert lhs==rhs;factorial+=1
            assert all(g[a][n]==g[n][a] for a in range(R+1) for n in range(R+1))
    records=[];maxentry=0;maxnorm=0
    for N,p in [(3,1),(4,2),(5,2)]:
        cs,vs,rs=gue_basis(N,p);dim=1<<N;charge=np.array([a.bit_count() for a in range(dim)]);Q=np.diag(charge)
        D=lambda Y:Q@Y-sum(c.conj().T@Y@c for c in cs)
        def twiceT(H):
            A=2*p*H@H-D(H@H)
            return 6*p*H@H@H-2*D(H@H@H)-3*(H@A+A@H)
        G=[]
        for V in vs:
            t0=twiceT(V)
            G.append(sum(int(r)*(twiceT(V+W)+twiceT(V-W)-2*t0) for W,r in zip(vs,rs)))
        G=np.array(G) # G_a=8 W1T3(V_a); covariance r_a/2.
        mean=np.zeros_like(Q,dtype=complex)
        norms=[0]*(N+1)
        for a,V in enumerate(vs):
            for c in range(a,len(vs)):
                K=V@G[c]-G[c]@V+vs[c]@G[a]-G[a]@vs[c]
                assert np.all(K.real==np.rint(K.real)) and np.all(K.imag==np.rint(K.imag))
                maxentry=max(maxentry,float(np.max(abs(K))));assert maxentry<2**20
                if a==c:mean+=int(rs[a])*K
                for m in range(p,N-p+1):
                    ix=np.flatnonzero(charge==m);block=K[np.ix_(ix,ix)]
                    real=block.real.astype(np.int64);imag=block.imag.astype(np.int64)
                    vv=int(np.sum(real*real+imag*imag));maxnorm=max(maxnorm,vv)
                    norms[m]+=(1 if a==c else 2)*int(rs[a]*rs[c])*vv
        assert np.all(mean==0)
        for m in range(p,N-p+1):
            direct=F(norms[m],512*b(N,m));formula=finite(N,p,m)
            assert direct==formula['M2'],(N,p,m,direct,formula)
            records.append(dict(N=N,p=p,m=m,direct_M2=str(direct),formula={k:str(v) for k,v in formula.items()}))
        print('CAR Gaussian contractions passed',N,p,flush=True)
    assert maxnorm<2**53
    # Shift-invariance and scalar rank: a shared scalar commutes with every channel.
    assert all(crossing(9,4,0,n)==1 for n in range(5))
    names=['code/check_cubic_dynamics.py','code/cubic_linear_dynamics.py','code/cubic_linear_channels.py','code/cubic_linear_limit.py','code/singlet_wick_moments.py','code/complex_syk_dictionary.py']
    path.write_text(json.dumps(dict(status='passed',exact_covariance_inversions=covariance,
      exact_factorial_moments=factorial,bounded_crossings=symmetry,car_records=records,
      max_integer_matrix_entry=maxentry,max_integer_norm_sum=maxnorm,
      arithmetic='Fraction identities; bounded Gaussian-integer CAR matrices with exact equality; Wick sum with Python integers',
      environment=dict(python=platform.python_version(),numpy=np.__version__),
      input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')
    print('passed',covariance,'crossings',factorial,'factorial moments',len(records),'CAR sectors',flush=True)
if __name__=='__main__':main()
