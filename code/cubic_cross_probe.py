"""Exact first dynamic cross moments of F with cubic coupling chaos.

Use fully polarized cubic CAR tensors and only the six contractions
from the cubic block to F,H,H. No Gaussian sampling or rank cutoff.
"""
from pathlib import Path
from fractions import Fraction as Fr
import hashlib,json,platform
from math import comb
import numpy as np
from singlet_wick_moments import gue_basis
from cubic_linear_dynamics import weights
from cubic_linear_channels import coeff
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-cubic-cross'
def case(N,p,m):
    cs,vs,rs=gue_basis(N,p);dim=1<<N;charge=np.array([a.bit_count() for a in range(dim)]);Q=np.diag(charge)
    D=lambda Y:Q@Y-sum(c.conj().T@Y@c for c in cs)
    def t2(H):
        A=2*p*H@H-D(H@H)
        return 6*p*H@H@H-2*D(H@H@H)-3*(H@A+A@H)
    G=[sum(int(r)*(t2(V+W)+t2(V-W)-2*t2(V)) for W,r in zip(vs,rs)) for V in vs]
    ix=np.flatnonzero(charge==m);iy=np.flatnonzero(charge==m-1)
    h=vs[:,ix[:,None],ix];g=vs[:,iy[:,None],iy];cm=[c[np.ix_(iy,ix)] for c in cs]
    E=lambda Y:sum(c.conj().T@Y@c for c in cm)
    avg=sum(w*coeff(N,p,m,a) for a,w in enumerate(weights(N,m,p)))
    k=np.array([avg.denominator*x[np.ix_(ix,ix)]-8*avg.numerator*v for x,v in zip(G,h)])
    scale=8*avg.denominator;count=len(h)
    hp=h[:,None]@h[None,:]+h[None,:]@h[:,None]
    gp=g[:,None]@g[None,:]+g[None,:]@g[:,None]
    ep=E(gp)
    sums={name:[0,0,0] for name in ['T','S','Z']};maxentry=0.;maxtrace=0.
    for a in range(count):
        for b in range(count):
            Z6=h[a]@hp[b]+h[b]@hp[a]+h@hp[a,b]
            U6=E(g[a]@gp[b]+g[b]@gp[a]+g@gp[a,b])
            S12=2*(2*p-m)*Z6+h[a]@ep[b]+ep[b]@h[a]+h[b]@ep[a]+ep[a]@h[b]+h@ep[a,b]+ep[a,b]@h
            T12=2*(3*p-m)*Z6+2*U6-3*S12
            for name,C in [('T',T12),('S',S12),('Z',2*Z6)]:
                maxentry=max(maxentry,float(np.max(abs(C))));assert maxentry<2**20
                assert np.all(C.real==np.rint(C.real)) and np.all(C.imag==np.rint(C.imag))
                words=[k[a]@C@h[b]@h,k[a]@h[b]@C@h,k[a]@h[b]@h@C]
                for u,word in enumerate(words):
                    tr=np.trace(word,axis1=-2,axis2=-1)
                    maxtrace=max(maxtrace,float(np.max(abs(tr))));assert maxtrace<2**40
                    assert np.all(tr.real==np.rint(tr.real)) and np.all(tr.imag==np.rint(tr.imag))
                    # Imaginary parts cancel in the complete ensemble; preserve separately.
                    value=sum(int(r)*complex(int(z.real),int(z.imag)) for r,z in zip(rs,tr))
                    sums[name][u]+=int(rs[a]*rs[b])*value
    out={}
    for name,vals in sums.items():
        assert all(abs(x)<2**53 and x.imag==0 for x in vals)
        B=[Fr(int(x.real),16*scale*comb(N,m)) for x in vals]
        assert B[0]==B[2]
        out[name]=dict(two_arc_moments=list(map(str,B)),cross_M2=str(B[0]-2*B[1]+B[2]))
    return dict(N=N,p=p,m=m,covariance='w=1',moments=out,max_integer_cubic_entry=maxentry,max_integer_trace_abs=maxtrace)
def main():
    OUT.mkdir(parents=True,exist_ok=True);target=OUT/'finite-cross.json'
    if target.exists():raise FileExistsError(target)
    rows=[]
    for N,p,m in [(3,1,2),(4,2,2),(5,2,3)]:
        row=case(N,p,m);rows.append(row);print(N,p,m,row['moments'],flush=True)
    names=['code/cubic_cross_probe.py','code/cubic_linear_channels.py','code/cubic_linear_dynamics.py','code/cubic_linear_limit.py','code/singlet_wick_moments.py','code/complex_syk_dictionary.py']
    target.write_text(json.dumps(dict(status='computed',scope='finite-N ensemble cross moments F H^u W3(C) H^(2-u); exact bounded Gaussian-integer contractions',
      records=rows,environment=dict(python=platform.python_version(),numpy=np.__version__),
      input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')
if __name__=='__main__':main()
