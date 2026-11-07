"""Exact chord polynomials verify all Laurent cancellations and h coefficients.

Independent of the real-time transfer-matrix implementation and of the
literature formula code. Hsrc has variance 1/h, q=exp(-h).
"""
from pathlib import Path
from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import factorial,comb
from fractions import Fraction
import json,hashlib,platform
import sympy as sp
from singlet_coefficient import pairings
from singlet_multimark import crossing
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-first-order'


def canonical(word):
    return min([word[i:]+word[:i] for i in range(len(word))]+[word[::-1][i:]+word[::-1][:i] for i in range(len(word))])


@lru_cache(None)
def moment(word):
    hs=tuple(i for i,x in enumerate(word) if x=='H');xs=tuple(i for i,x in enumerate(word) if x=='X')
    out=Counter()
    for hp in pairings(hs):
        hh=sum(crossing(a,b) for a,b in combinations(hp,2))
        for xp in pairings(xs):
            degree=hh+2*sum(crossing(a,b) for a in hp for b in xp)+4*sum(crossing(a,b) for a,b in combinations(xp,2))
            out[degree]+=1
    return dict(out)


def two(n):
    out=Counter()
    for k in range(n+1):
        weight=comb(n,k)*(-1)**k
        for degree,val in moment(canonical('H'*(n-k)+'X'+'H'*k+'X')).items():out[degree]+=weight*val
    return {k:v for k,v in out.items() if v}


def four(n,times):
    t1,t2,t3,t4=times;arcs=[t4-t3,t3-t2,t2-t1,t1-t4]
    out=Counter()
    for a in range(n+1):
      for b in range(n-a+1):
       for c in range(n-a-b+1):
        d=n-a-b-c;ks=[a,b,c,d]
        # Compute multinomial once using exact integer arithmetic.
        w=factorial(n)//(factorial(a)*factorial(b)*factorial(c)*factorial(d))
        for k,arc in zip(ks,arcs):w*=arc**k
        if not w:continue
        word='H'*d+'X'+'H'*c+'X'+'H'*b+'X'+'H'*a+'X'
        for degree,val in moment(canonical(word)).items():out[degree]+=w*val
    return dict(out)


def disconnected(n,times,twos):
    out=Counter()
    for xp in pairings(tuple(range(4))):
        (a,b),(c,d)=xp;dx=times[a]-times[b];dy=times[c]-times[d]
        for k in range(0,n+1,2):
            weight=comb(n,k)*dx**k*dy**(n-k)
            if not weight:continue
            for i,x in twos[k].items():
                for j,y in twos[n-k].items():out[i+j]+=weight*x*y
    return dict(out)


def hcoef(poly,n,power):
    k=n//2+power
    if k<0:return Fraction(0)
    return Fraction((-1)**(n//2)*sum(v*(-d)**k for d,v in poly.items()),factorial(n)*factorial(k))


def main():
    OUT.mkdir(exist_ok=True);path=OUT/'exact-moments.json'
    if path.exists():raise FileExistsError(path)
    u=sp.symbols('u',real=True)
    cases=[(1,0,1,0),(1,1,0,0),(2,1,-1,0),(0,1,3,-1)]
    twos={n:two(n) for n in range(0,11,2)}
    records=[];cancellations=0
    for times in cases:
        maxn=10 if times in cases[:2] else 8
        a,b,c,d=[v*u for v in times]
        g=lambda x:sp.sech(x)**4
        formula=-4*(g(a-b)*g(c-d)*sp.tanh(a-b)*sp.tanh(c-d)
            +g(b-c)*g(d-a)*sp.tanh(b-c)*sp.tanh(d-a)
            +g(a-c)*g(b-d)*sp.cosh(a+c-b-d)/(sp.cosh(a-c)*sp.cosh(b-d)))
        series=sp.series(formula,u,0,maxn+1).removeO().expand()
        for n in range(0,maxn+1,2):
            f=four(n,times);disc=disconnected(n,times,twos);p=Counter(f)
            for k,v in disc.items():p[k]-=v
            p={k:v for k,v in p.items() if v}
            zero={str(j):str(hcoef(p,n,j)) for j in range(-n//2,1)}
            assert all(Fraction(v)==0 for v in zero.values()),(times,n,zero)
            cancellations+=len(zero)
            exact=hcoef(p,n,1);expected=Fraction(str(series.coeff(u,n)))
            assert exact==expected,(times,n,exact,expected)
            records.append(dict(times=times,u_power=n,connected_q_polynomial=p,laurent_cancellations=zero,
                first_h_coefficient=str(exact),closed_form_coefficient=str(expected)))
        print('moments checked',times,maxn,flush=True)
    cform=8*(sp.cosh(2*u)+(4*sp.tanh(u)**2-1)*sp.sech(u)**8)
    curvature=sp.series(cform,u,0,12)
    payload=dict(status='passed',convention='q=e^-h, Hsrc variance 1/h, d=2, physical times t_j=times_j*u/(2 Js)',
        cases=len(cases),taylor_coefficients=len(records),laurent_cancellations=cancellations,
        commutator_per_h_series=str(curvature),records=records,
        environment=dict(python=platform.python_version(),sympy=sp.__version__),
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['code/check_singlet_first_order_moments.py','code/singlet_coefficient.py','code/singlet_multimark.py']})
    path.write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps(dict(status='passed',coefficients=len(records),laurent_cancellations=cancellations,series=str(curvature))))


if __name__=='__main__':main()
