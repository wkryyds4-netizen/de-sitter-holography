"""Linear chaos of sym(H A), including exact centering and projection."""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib
import numpy as np
from cubic_linear_channels import choose,coeff
from cubic_linear_dynamics import weights
from singlet_wick_moments import gue_basis
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-cubic-cross'
def channel(N,p,m,a):
    C=choose(m,p)*choose(N-m+p,p);Cm=choose(m-1,p)*choose(N-m+p+1,p)
    kap=choose(m-a,p)*choose(N-m+p-a,p)
    z=F(2*C+kap)
    const=F(2*p-m)*C+m*Cm+F(p*p,m)*C
    J=F(p*(p*(N-m+1)-a*(N-a+1)),(m-a)*(N-m-a+1)) if m>a else F(0)
    # At m=a=p use the nonsingular identity instead of a removable ratio.
    zm=2*Cm+choose(m-1-a,p)*choose(N-m+1+p-a,p)
    ss=(F(3*p-m)*z+(m-p)*zm-coeff(N,p,m,a))/3
    if m>a:assert ss==const+kap*J
    return z,ss
def main():
    OUT.mkdir(parents=True,exist_ok=True);target=OUT/'linear-channels.json'
    if target.exists():raise FileExistsError(target)
    rows=[]
    for N,p in [(3,1),(4,2),(5,2),(6,2)]:
        cs,vs,rs=gue_basis(N,p);dim=1<<N;charges=np.array([a.bit_count() for a in range(dim)]);Q=np.diag(charges)
        D=lambda Y:Q@Y-sum(c.conj().T@Y@c for c in cs)
        def rawS2(H):
            A=2*p*H@H-D(H@H)
            return H@A+A@H
        for a in range(min(p,N//2)+1):
            O=np.eye(dim,dtype=complex)
            for i in range(a):O=O@cs[i].conj().T
            for i in range(a):O=O@cs[N-a+i]
            if a:O+=O.conj().T
            h=np.diag([choose(int(m)-a,p-a) for m in charges])@O
            s0=rawS2(h)
            lin8=sum(int(r)*(rawS2(h+v)+rawS2(h-v)-2*s0) for v,r in zip(vs,rs))
            for m in range(max(p,a),N-a+1):
                ix=np.flatnonzero(charges==m);block=h[np.ix_(ix,ix)];actual=lin8[np.ix_(ix,ix)]
                z,s=channel(N,p,m,a)
                assert np.all(actual*s.denominator==8*s.numerator*block)
                rows.append(dict(N=N,p=p,m=m,a=a,Z_linear=str(z),S_linear=str(s)))
        print('S channels',N,p,'passed',flush=True)
    names=['code/cubic_cross_channels.py','code/cubic_linear_channels.py','code/cubic_linear_dynamics.py','code/cubic_linear_limit.py','code/singlet_wick_moments.py','code/complex_syk_dictionary.py']
    target.write_text(json.dumps(dict(status='passed',checks=len(rows),records=rows,
      input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')
if __name__=='__main__':main()
