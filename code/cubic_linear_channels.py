"""Exact finite-sector channel coefficients for the first Gaussian chaos of T3."""
from pathlib import Path
from math import comb
from fractions import Fraction
import json,hashlib
import numpy as np
from singlet_wick_moments import gue_basis
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-cubic-linear'
def choose(n,k):return comb(n,k) if 0<=k<=n else 0
def coeff(N,p,m,nu):
    """Coupling covariance w=1; allowed nu<=min(p,m,N-m)."""
    C=lambda m:choose(m,p)*choose(N-m+p,p)
    kap=lambda m:choose(m-nu,p)*choose(N-m+p-nu,p)
    jc=choose(m-1,p)*choose(N-m+p,p)
    jl=choose(m-nu,p)*choose(N-m+p-nu+1,p)
    jh=choose(m-nu-1,p)*choose(N-m+p-nu,p)
    al=Fraction(nu*(N-nu-p+1),N-2*nu+1)
    base=-(m+2*p)*C(m-1)+(m-3*p)*C(m)+3*p*jc
    return Fraction(base+(m-p)*kap(m-1)-m*kap(m))+3*(al*jl+(p-al)*jh)
def main():
    target=OUT/'channels.json'
    if target.exists():raise FileExistsError(target)
    rows=[];maxentry=0
    for n,p in [(3,1),(4,2),(5,2),(6,2),(6,3)]:
        cs,vs,rs=gue_basis(n,p);dim=1<<n;charges=np.array([a.bit_count() for a in range(dim)])
        Q=np.diag(charges)
        def D(Y):return Q@Y-sum(c.conj().T@Y@c for c in cs)
        def twiceT(H):
            A=2*p*H@H-D(H@H)
            return 6*p*H@H@H-2*D(H@H@H)-3*(H@A+A@H)
        for nu in range(min(p,n//2)+1):
            O=np.eye(dim,dtype=complex)
            for i in range(nu):O=O@cs[i].conj().T
            for i in range(nu):O=O@cs[n-nu+i]
            if nu:O=O+O.conj().T
            f=np.diag([choose(int(m)-nu,p-nu) for m in charges]);H=f@O;t0=twiceT(H)
            linear8=sum(int(r)*(twiceT(H+v)+twiceT(H-v)-2*t0) for v,r in zip(vs,rs))
            assert np.all(linear8.real==np.rint(linear8.real)) and np.all(linear8.imag==np.rint(linear8.imag))
            maxentry=max(maxentry,float(np.max(abs(linear8))))
            assert maxentry<2**40
            for m in range(max(p,nu),n-nu+1):
                ix=np.flatnonzero(charges==m);h=H[np.ix_(ix,ix)];f8=linear8[np.ix_(ix,ix)]
                prediction=coeff(n,p,m,nu)
                assert np.all(f8*prediction.denominator==8*prediction.numerator*h),(n,p,m,nu,prediction,float(np.linalg.norm(f8-8*float(prediction)*h)))
                rows.append(dict(N=n,p=p,m=m,nu=nu,coefficient=str(prediction)))
        print('channels',n,p,'passed',flush=True)
    target.write_text(json.dumps(dict(status='passed',covariance='w=1; standard CAR; exact Gaussian-integer equality checks with bounded values',
      rows=rows,channel_checks=len(rows),max_integer_entry=maxentry,
      input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['code/cubic_linear_channels.py','code/singlet_wick_moments.py','code/complex_syk_dictionary.py']}),indent=2)+'\n')
if __name__=='__main__':main()
