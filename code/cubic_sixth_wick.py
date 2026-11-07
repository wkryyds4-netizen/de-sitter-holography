"""Finite CAR six-slot Wick diagrams, with same-coupling F at two slots."""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import json,hashlib
import numpy as np
import sympy as sy
from singlet_wick_moments import gue_basis
from cubic_linear_dynamics import finite,weights
from cubic_linear_channels import coeff
from singlet_coefficient import pairings
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-linear-kernel'
def crossings(pa):
    return sum(a<c<b<d or c<a<d<b for j,(a,b) in enumerate(pa) for c,d in pa[j+1:])
def main():
    OUT.mkdir(parents=True,exist_ok=True);path=OUT/'sixth-wick.json'
    if path.exists():raise FileExistsError(path)
    q=sy.Symbol('q');r=sy.Symbol('r');polys={};records=[];max_bound=0
    for u in range(5):
        v=4-u
        val=sum(q**crossings(pa)*r**sum(a<u<=b for a,b in pa) for pa in pairings(tuple(range(4))))
        polys[f'{u},{v}']=sy.expand(val)
    fourth=sy.factor(sum((-1)**k*comb(4,k)*polys[f'{4-k},{k}'] for k in range(5)))
    assert sy.factor(fourth.subs(r,q)-2*(1-q)**2*(3*q+5))==0
    for N,p,m in [(3,1,2),(4,2,2),(5,2,3)]:
        cs,vs,rs=gue_basis(N,p);dim=1<<N;charge=np.array([a.bit_count() for a in range(dim)]);Q=np.diag(charge)
        D=lambda Y:Q@Y-sum(c.conj().T@Y@c for c in cs)
        def t2(H):
            A=2*p*H@H-D(H@H)
            return 6*p*H@H@H-2*D(H@H@H)-3*(H@A+A@H)
        G=[sum(int(r)*(t2(V+W)+t2(V-W)-2*t2(V)) for W,r in zip(vs,rs)) for V in vs]
        avg=sum(w*coeff(N,p,m,a) for a,w in enumerate(weights(N,m,p)))
        ix=np.flatnonzero(charge==m);h=vs[:,ix[:,None],ix]
        k=np.array([avg.denominator*g[np.ix_(ix,ix)]-8*avg.numerator*v for g,v in zip(G,h)])
        scale=8*avg.denominator;cov={}
        for a,x in [('H',h),('F',k)]:
            for b,y in [('H',h),('F',k)]:
                tensor=np.einsum('a,aij,akl->ijkl',rs,x,y)
                assert np.all(tensor.imag==0) and np.all(tensor.real==np.rint(tensor.real))
                cov[a+b]=tensor.real.astype(np.int64)
        degree={};moments={}
        for order in [2,4]:
            for u in range(order+1):
                word='F'+'H'*u+'F'+'H'*(order-u);L=len(word);slots='abcdef'[:L];terms=[]
                for pa in pairings(tuple(range(L))):
                    tensors=[cov[word[a]+word[b]] for a,b in pa]
                    labels=[slots[a]+slots[(a+1)%L]+slots[b]+slots[(b+1)%L] for a,b in pa]
                    bound=comb(N,m)**L
                    for tensor in tensors:bound*=int(np.max(abs(tensor)))
                    assert bound<2**63;max_bound=max(max_bound,bound)
                    raw=int(np.einsum(','.join(labels)+'->',*tensors,optimize=True))
                    val=F(raw,2**len(pa)*scale**2*comb(N,m))
                    # Independent basis-coordinate contraction for the first nontrivial system.
                    # Degree two is independently known from the rank map, used below.
                    terms.append(dict(pairing=pa,value=str(val),mixed_pairs=sum(word[a]!=word[b] for a,b in pa)))
                total=sum(F(t['value']) for t in terms)
                moments[order,u]=total
                degree[f'{order},{u}']=dict(word=word,value=str(total),pairings=terms,
                  same_disorder_pairing_sum=str(sum(F(t['value']) for t in terms if t['mixed_pairs'])))
        ref=finite(N,p,m)
        assert moments[2,0]==ref['B20'] and moments[2,1]==ref['B11']
        M4=sum((-1)**k*comb(4,k)*moments[4,k] for k in range(5))
        assert M4>=0
        if ref['M2']==0:assert M4==0
        records.append(dict(N=N,p=p,m=m,variance=str(ref['norm']),M4=str(M4),
          normalized_M4=str(M4/ref['norm']),moments=degree))
        print('six-slot finite Wick',N,p,m,'M4',M4,flush=True)
    names=['code/cubic_sixth_wick.py','code/cubic_linear_dynamics.py','code/cubic_linear_channels.py','code/cubic_linear_limit.py','code/singlet_wick_moments.py','code/complex_syk_dictionary.py','code/singlet_coefficient.py']
    path.write_text(json.dumps(dict(status='passed',finite_records=records,
      predicted_probe_degree4={k:str(v) for k,v in polys.items()},predicted_M4_over_C2=str(fourth.subs(r,q).factor()),
      finite_scope='Exact CAR covariance tensors and all 15 Gaussian pairings per six-slot word; no claim finite size equals DS',
      int64_absolute_sum_bound=max_bound,
      input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')
if __name__=='__main__':main()
