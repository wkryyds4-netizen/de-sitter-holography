"""Check charged covariance, channel weights and a nonconservation witness."""
from pathlib import Path
from math import comb
from fractions import Fraction
import json,hashlib
import numpy as np
import sympy as s
from singlet_wick_moments import gue_basis
from cubic_linear_channels import choose,coeff
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-cubic-linear'
def main():
    target=OUT/'validation.json'
    if target.exists():raise FileExistsError(target)
    charged=[];weight_checks=0
    for N,p in [(3,1),(4,2),(5,2),(6,2)]:
        cs,vs,rs=gue_basis(N,p);dim=1<<N;charges=np.array([a.bit_count() for a in range(dim)])
        for rank_s in range(1,min(p+1,(N+1)//2)+1):
            rank_r=rank_s-1;Y=np.eye(dim,dtype=complex)
            for i in range(rank_r):Y=Y@cs[i].conj().T
            for i in range(rank_s):Y=Y@cs[N-rank_s+i]
            mapped2=sum(int(r)*v@Y@v for v,r in zip(vs,rs))
            for m in range(rank_s,N-rank_r+1):
                rows=np.flatnonzero(charges==m-1);cols=np.flatnonzero(charges==m)
                expected=choose(m-rank_s,p)*choose(N-m+p-rank_r,p)
                assert np.all(mapped2[np.ix_(rows,cols)]==2*expected*Y[np.ix_(rows,cols)])
                charged.append(dict(N=N,p=p,m=m,r=rank_r,s=rank_s,eigenvalue=expected))
    for N in range(4,15):
        for p in range(1,min(3,N//2)+1):
            for m in range(p,N-p+1):
                C=choose(m,p)*choose(N-m+p,p)
                probabilities=[Fraction((choose(N,nu)**2-choose(N,nu-1)**2)*choose(m-nu,m-p)*choose(N-p-nu,m-p),choose(N,m)*C) for nu in range(p+1)]
                assert sum(probabilities)==1
                for nu in range(1,p+1):
                    ratio=Fraction(nu**2,(N-nu+1)**2)*(1-Fraction((nu-1)**2,(N-nu+2)**2))/(1-Fraction(nu**2,(N-nu+1)**2))*Fraction(m-nu+1,p-nu+1)*Fraction(N-p-nu+1,N-m-nu+1)
                    assert probabilities[nu-1]/probabilities[nu]==ratio
                weight_checks+=1
    # Exact rational proof of the charged low/high projection coefficients.
    N,m,nu,p=s.symbols('N m nu p')
    fu=s.cancel((N-m-nu+1)/(N-2*nu+1))
    vnorm=(N-m+1)*(m-nu)/(N-m-nu+1)
    fnorm=s.factor(m+vnorm-2*(m-nu))
    identities=[
      fnorm-nu*(N-2*nu+1)/(N-m-nu+1),
      fu*fnorm-nu,
      fu-(m-p)/(m-nu)*(fu-1)-(N-nu-p+1)/(N-2*nu+1)]
    assert all(s.factor(x)==0 for x in identities)
    # Replay deterministic nonconservation witness with exact integer-valued matrices.
    N,p=5,2;cs,vs,rs=gue_basis(N,p);dim=1<<N;charge=np.array([a.bit_count() for a in range(dim)]);Q=np.diag(charge)
    D=lambda Y:Q@Y-sum(c.conj().T@Y@c for c in cs)
    def t2(H):
        A=2*p*H@H-D(H@H)
        return 6*p*H@H@H-2*D(H@H@H)-3*(H@A+A@H)
    x=np.array([(7*i+3)%11-5 for i in range(len(vs))]);H=np.einsum('a,aij->ij',x,vs);t0=t2(H)
    lin8=sum(int(r)*(t2(H+v)+t2(H-v)-2*t0) for v,r in zip(vs,rs))
    comm8=H@lin8-lin8@H;ix=np.flatnonzero(charge==3);block=comm8[np.ix_(ix,ix)]
    assert np.all(block.real==np.rint(block.real)) and np.all(block.imag==np.rint(block.imag))
    rr=block.real.astype(np.int64);ii=block.imag.astype(np.int64)
    norm2=Fraction(int(np.sum(rr*rr+ii*ii)),64*len(ix));assert norm2>0
    target.write_text(json.dumps(dict(status='passed',charged_covariance_checks=len(charged),charged_records=charged,
      exact_weight_normalizations=weight_checks,symbolic_projection_identities=3,
      witness=dict(N=5,p=2,m=3,coefficients=x.tolist(),covariance='w=1',normalized_commutator_squared=str(norm2),
        meaning='[H,W1 T3] is not zero; subtracting c(Q)H leaves this commutator unchanged'),
      input_sha256={path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in ['code/check_cubic_linear.py','code/cubic_linear_channels.py','code/singlet_wick_moments.py','code/complex_syk_dictionary.py']}),indent=2)+'\n')
    print('passed',len(charged),'charged,',weight_checks,'weights; witness',norm2)
if __name__=='__main__':main()
