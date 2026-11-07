"""Exact mixed-q commutator primary, and independent small-h covariance checks."""
from pathlib import Path
from collections import Counter
from functools import lru_cache
from itertools import combinations,product
from math import factorial
from fractions import Fraction
import sympy as s
import json,hashlib,platform
from singlet_coefficient import pairings
from singlet_multimark import crossing
ROOT=Path(__file__).resolve().parent.parent;OUT=ROOT/'results/susskind-commutator-primary'
q=s.Symbol('q')

@lru_cache(None)
def moment(word):
    slots={c:tuple(i for i,z in enumerate(word) if z==c) for c in 'HXRP'}
    if any(len(v)%2 for v in slots.values()):return {}
    total=Counter()
    for ps in product(*(pairings(slots[c]) for c in 'HXRP')):
        ch=[(p,d) for pairs,d in zip(ps,[1,2,3,5]) for p in pairs]
        e=sum(a*b*crossing(c,d) for (c,a),(d,b) in combinations(ch,2));total[e]+=1
    return dict(total)
def expr(poly):return s.factor(sum(v*q**k for k,v in poly.items()))
def addword(total,word,weight):
    for k,v in moment(word).items():total[k]+=weight*v
def finite_primary():
    rows=[]
    for n in range(0,9,2):
        for a in range(n+1):
            b=n-a;raw=Counter();cross=Counter()
            for left,wl in [('XR',1),('RX',-1)]:
                for right,wr in [('XR',1),('RX',-1)]:addword(raw,left+'H'*a+right+'H'*b,-wl*wr)
                for right in ['XR','RX']:addword(cross,left+'H'*a+right+'H'*b,wl)
            target=2*(1-q**6)*expr(moment('P'+'H'*a+'P'+'H'*b))
            assert s.expand(expr(raw)-target)==0 and expr(cross)==0
            rows.append(dict(H_left=a,H_right=b,raw_KK=str(expr(raw)),raw_KS=str(expr(cross))))
    high=[]
    for n in [2,4,6]:
        total=Counter()
        for bits in product(range(2),repeat=n):
            word=''.join('XR' if bit==0 else 'RX' for bit in bits)
            addword(total,word,(-1)**(sum(bits)+n//2))
        raw=expr(total);norm=s.factor(raw/(2*(1-q**6))**(n//2))
        high.append(dict(order=n,raw=str(raw),normalized=str(norm),q1_limit=str(s.limit(norm,q,1))))
    fourth=s.sympify(high[1]['normalized']);delta=s.factor(fourth-(2+q**25))
    h=s.Symbol('h',positive=True);lead=s.limit(delta.subs(q,s.exp(-h))/h**2,h,0)
    assert [r['q1_limit'] for r in high]==['1','3','15']
    return dict(two_arc_cases=rows,even_moments=high,fourth_probe_difference=str(delta),fourth_difference_over_h2_limit=str(lead),
        cubic_vertex='tau(K X R) = i sqrt((1-q^6)/2); tau(K R X) = -i sqrt((1-q^6)/2)')

def four(n,labels,times):
    arcs=[times[(j+1)%4]-times[j] for j in range(4)];out=Counter()
    for a in range(n+1):
      for b in range(n-a+1):
       for c in range(n-a-b+1):
        ks=[a,b,c,n-a-b-c];weight=factorial(n)
        for k,arc in zip(ks,arcs):weight=Fraction(weight,factorial(k))*arc**k
        if not weight:continue
        word=''.join(label+'H'*k for label,k in zip(labels,ks));addword(out,word,weight)
    return dict(out)
def hcoef(poly,n,power):
    k=n//2+power
    return Fraction((-1)**(n//2)*sum(v*(-d)**k for d,v in poly.items()),factorial(n)*factorial(k))
def covariance_checks():
    z=s.Symbol('z',real=True);rows=[];cancellations=0
    for ray in [(1,0,1,0),(1,1,0,0),(2,1,-1,0),(0,1,3,-1)]:
        a,b,c,d=[z*v for v in ray];dx=a-c;dy=b-d;A=a+c-b-d
        formula=12*s.sech(dx)**5*s.sech(dy)**7*(s.cosh(A)+s.sinh(dx)*s.sinh(dy))
        series=s.series(formula,z,0,9).removeO().expand()
        for n in range(0,9,2):
            total=Counter()
            for labels,ts,sign in [('XRXR',ray,-1),('XRRX',(ray[0],ray[1],ray[3],ray[2]),1),
                ('RXXR',(ray[1],ray[0],ray[2],ray[3]),1),('RXRX',(ray[1],ray[0],ray[3],ray[2]),-1)]:
                for k,v in four(n,labels,ts).items():total[k]+=sign*v
            zero=[hcoef(total,n,j) for j in range(-n//2,1)];assert all(x==0 for x in zero);cancellations+=len(zero)
            actual=hcoef(total,n,1);expected=Fraction(str(series.coeff(z,n)));assert actual==expected,(ray,n,actual,expected)
            rows.append(dict(ray=ray,power=n,coefficient=str(actual),zero_laurent_powers=list(range(-n//2,1))))
        print('covariance checked',ray,flush=True)
    return dict(records=rows,laurent_cancellations=cancellations)
def main():
    OUT.mkdir(parents=True,exist_ok=True);target=OUT/'exact.json'
    if target.exists():raise FileExistsError(target)
    primary=finite_primary();print('primary',primary['fourth_probe_difference'],flush=True)
    covariance=covariance_checks()
    inputs=['code/commutator_primary.py','code/singlet_coefficient.py','code/singlet_multimark.py']
    data=dict(status='passed',conventions='q=exp(-h), h=4lambda, u=2Js t; Hsrc variance=1/h; K=-i[X,Rhat]/sqrt(2(1-q^6))',
        primary=primary,covariance=covariance,environment=dict(python=platform.python_version(),sympy=s.__version__),
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    target.write_text(json.dumps(data,indent=2)+'\n');print('saved',len(primary['two_arc_cases']),len(covariance['records']),flush=True)
if __name__=='__main__':main()
