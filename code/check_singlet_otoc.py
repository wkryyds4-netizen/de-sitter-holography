"""Independent rational chord verification of BINT channel normalization."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from math import factorial
import hashlib,json,platform
import numpy as np
import sympy as sp
from singlet_coefficient import pairings
from singlet_multimark import crossing,probe_moment

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-singlet-otoc'


def compositions(n):
    for a in range(n+1):
        for b in range(n-a+1):
            for c in range(n-a-b+1):yield a,b,c,n-a-b-c


def direct(ks,q):
    k1,k2,k3,k4=ks
    word='H'*k4+'X'+'H'*k3+'X'+'H'*k2+'X'+'H'*k1+'X'
    hs=tuple(i for i,c in enumerate(word) if c=='H');xs=[i for i,c in enumerate(word) if c=='X']
    out=[]
    for xp in [((xs[0],xs[1]),(xs[2],xs[3])),((xs[0],xs[2]),(xs[1],xs[3])),((xs[0],xs[3]),(xs[1],xs[2]))]:
        value=F(0)
        for hp in pairings(hs):
            degree=sum(crossing(a,b) for a,b in combinations(hp,2))+4*crossing(*xp)
            degree+=2*sum(crossing(a,b) for a in hp for b in xp)
            value+=q**degree
        out.append(value)
    return out


def transfer(ks,q):
    n=sum(ks);L=n+2;r=q*q
    T=np.zeros((L,L),dtype=object)
    for j in range(L-1):T[j+1,j]=F(1);T[j,j+1]=(1-q**(j+1))/(1-q)
    powers=[np.eye(L,dtype=object)]
    for j in range(n):powers.append(powers[-1]@T)
    p1,p2,p3,p4=[powers[k] for k in ks]
    f=[F(1)];g=[F(1)]
    for j in range(1,L):f.append(f[-1]*(1-q**j));g.append(g[-1]*(1-r*r*q**(j-1)))
    crossed=F(0)
    for i in range(min(ks[0],ks[3])+1):
        for l in range(ks[0]-i+1):
            if not p1[l+i,0]:continue
            P=f[l+i]*g[i]/(f[l]*f[i])
            for m in range(ks[3]-i+1):
                if not p4[0,m+i]:continue
                for h in range(n+1):
                    crossed+=P*r**(i+h+l+m)*p4[0,m+i]*p3[m,h]*p2[h,l]*p1[l+i,0]
    # Eq (4.8): central power k2+k4; leaf powers k1,k3.
    def unc(k):
        a,b,c,d=k
        return sum(powers[c][0,l]*r**l*powers[b+d][l,m]*r**m*powers[a][m,0] for l in range(L) for m in range(L))
    return [unc(ks),q**4*crossed,unc((ks[1],ks[2],ks[3],ks[0]))]


def curvature():
    q=sp.symbols('q');r=q*q
    def moment(word):return sum(v*q**k[1] for k,v in probe_moment(word).items())
    result=sp.factor(2*moment('HHXXXX')+6*moment('HXXHXX')-8*moment('HXHXXX'))
    expected=2*(1-r)**2*(5+2*r+3*r*r)
    assert sp.expand(result-expected)==0
    return str(result)


def main():
    OUT.mkdir(exist_ok=True);target=OUT/'exact-channels.json'
    if target.exists():raise FileExistsError(target)
    rows=[]
    for q in [F(2,5),F(3,4)]:
        for n in [0,2,4,6]:
            for ks in compositions(n):
                a=direct(ks,q);b=transfer(ks,q)
                assert a==b,(q,ks,a,b)
                rows.append(dict(q=str(q),powers=ks,channels=list(map(str,a))))
        print('q',q,'all channels passed',flush=True)
    coeff=curvature()
    inputs=['code/check_singlet_otoc.py','code/singlet_coefficient.py','code/singlet_multimark.py']
    data=dict(status='passed',moment_cases=len(rows),rational_channel_comparisons=3*len(rows),
        commutator_t2_over_sigma2=coeff,records=rows,environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=sp.__version__),
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    target.write_text(json.dumps(data,indent=2)+'\n');print(coeff,flush=True)

if __name__=='__main__':main()
