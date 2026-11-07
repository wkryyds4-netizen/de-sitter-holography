"""All endpoint templates for the quartic residual two-energy kernel."""
from pathlib import Path
from itertools import combinations,permutations
from math import factorial
import json,hashlib,platform
import sympy as s
from singlet_spectral_graphs import graph

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-quartic'
q=s.Symbol('q');d3=q*q-q+1;d5=sum(q**j for j in range(5))
den=(1+q*q)*d3*d5
mpairs=list(combinations(range(4),2))
ww=[q**3*(1+q*q),-q*(1+q*q)*d5,den,
    q*q*(1+q+q*q)*(1+q*q),-q*(1+q*q)*d5,q**3*(1+q*q)]
constant=2*(q-1)**3*(q+1)*(1+q+q*q)*d3


def cross(a,b,c,d):
    a,b=sorted((a,b));c,d=sorted((c,d))
    return int(a<c<b<d or c<a<d<b)


def orbit(k,colors,split):
    A=list(range(4));B=list(range(4+split,8+split))
    ext=list(range(4,4+split))+list(range(8+split,8+len(colors)))
    ae=[e for e,c in zip(ext,colors) if c=='A'];be=[e for e,c in zip(ext,colors) if c=='B']
    seed=list(zip(A[:k],B[:k]))+list(zip(A[k:],ae))+list(zip(B[k:],be))
    words=set()
    for ap in permutations(A):
        for bp in permutations(B):
            pm=dict(zip(A,ap))|dict(zip(B,bp))
            words.add(tuple(sorted(tuple(sorted((pm.get(a,a),pm.get(b,b)))) for a,b in seed)))
    assert len(words)==576//factorial(k)
    return A,B,sorted(words)


def polynomial_dict(v):
    return {int(k[0]):int(c) for k,c in s.Poly(s.expand(v),q).terms()}


def accumulate(dest,poly,factor,shift):
    if factor:
        for power,value in poly.items():dest[power+shift]=dest.get(power+shift,0)+factor*value


def main():
    target=OUT/'orbit-certificate.json'
    if target.exists():raise FileExistsError(target)
    allw=ww+[constant]
    product=[[polynomial_dict(a*b) for b in allw] for a in allw]
    records=[];fails=[];nwords=checked=0
    for k in range(5):
        n=8-2*k
        for aset in combinations(range(n),4-k):
            colors=''.join('A' if j in aset else 'B' for j in range(n))
            for split in range(n+1):
                A,B,words=orbit(k,colors,split);nwords+=len(words)
                tensor={};delta_tensor={}
                for word_i,pa in enumerate(words):
                    pi={a:b for a,b in pa}|{b:a for a,b in pa}
                    cr=sum(cross(a,b,c,d) for j,(a,b) in enumerate(pa) for c,d in pa[j+1:])
                    aa=[2-4*cross(A[a],pi[A[a]],A[b],pi[A[b]]) for a,b in mpairs]+[1]
                    bb=[2-4*cross(B[a],pi[B[a]],B[b],pi[B[b]]) for a,b in mpairs]+[1]
                    t=tensor.setdefault(cr,[[0]*7 for _ in range(7)])
                    for i in range(7):
                        for j in range(7):t[i][j]+=aa[i]*bb[j]
                    dt=delta_tensor.setdefault(cr,[[0]*6 for _ in range(6)])
                    for i,(a,b) in enumerate(mpairs):
                        for j,(c,d) in enumerate(mpairs):
                            dd=4*int(pi[A[a]]==B[d] and pi[A[b]]==B[c]);dt[i][j]+=dd
                            if aset==tuple(range(4-k)) and split==0 and word_i in [0,len(words)-1]:
                                g=graph(8+n,pa,((A[a],A[b]),(B[c],B[d])))
                                assert (g['crossings'],g['alpha'],g['beta'],g['delta'])==(cr,aa[i],bb[j],dd)
                                checked+=1
                linear={};quadratic={}
                for cr,t in tensor.items():
                    for i in range(7):
                        for j in range(7):accumulate(quadratic,product[i][j],t[i][j],cr)
                for cr,t in delta_tensor.items():
                    for i in range(6):
                        for j in range(6):accumulate(linear,product[i][j],t[i][j],cr)
                linear={n:c for n,c in linear.items() if c};quadratic={n:c for n,c in quadratic.items() if c}
                passed=(not linear and not quadratic) if k<4 else True
                if not passed:fails.append(dict(bridges=k,colors=colors,split=split,linear=linear,quadratic=quadratic))
                records.append(dict(bridges=k,colors=colors,split=split,words=len(words),vanishes=not linear and not quadratic,
                                    linear=linear if k==4 or not passed else {},quadratic=quadratic if k==4 or not passed else {}))
        print('bridges',k,'templates',len(records),'counterexamples',len(fails),flush=True)
    names=['code/quartic_orbits.py','code/singlet_spectral_graphs.py','code/singlet_coefficient.py']
    result=dict(status='passed' if not fails else 'counterexamples',denominator=str(s.expand(den)),
                numerator_weights=list(map(str,ww)),constant_numerator=str(constant),templates=len(records),
                orbit_words=nwords,minor_checks=checked,records=records,counterexamples=fails,
                environment=dict(python=platform.python_version(),sympy=s.__version__),
                input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
