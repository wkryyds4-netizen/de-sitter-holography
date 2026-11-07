"""Common-disorder primitive H/W/J moments with all local source masks.

C=1. W=(p-L-lambda)H; J=(2p-D-kappa)Z2 with
kappa=2lambda(1-q)/(1+q). The latter differs from the manuscript J2
only by its vanishing charge-dependent scalar counterterm.
Returns exact (lambda,q) numerator polynomials after multiplying by
(1+q)^(number of J insertions). No source order is silently discarded.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from math import comb

from memory_body_source_chords import exponent
from singlet_coefficient import pairings

def pmul(a,b):
    out=Counter()
    for (i,j),x in a.items():
        for (k,l),y in b.items():out[i+k,j+l]+=x*y
    return {k:v for k,v in out.items() if v}

def padd(*polys):
    out=Counter()
    for p in polys:out.update(p)
    return {k:v for k,v in out.items() if v}

def ppower(p,n):
    out={(0,0):F(1)}
    for _ in range(n):out=pmul(out,p)
    return out

def crossing(a,b):
    x,y=a;z,w=b
    return int(x<z<y<w or z<x<w<y)

def layout(word):
    cursor=0;blocks=[];wick=[];types=[]
    for letter in word:
        assert letter in 'HWJ'
        size=2 if letter=='J' else 1
        if letter!='H':
            blocks.append((cursor,cursor+size,'D' if letter=='J' else 'L'))
            types.append(letter)
            wick.append((cursor,cursor+size))
        cursor+=size
    return cursor,tuple(blocks),tuple(wick),tuple(types)

def shifted_weight(kappa,types):
    # Rescale every J source by (1+q), then subtract 2lambda(1-q).
    # Thus the final denominator is exactly (1+q)^number_J.
    marks=len(types);coeff={}
    for mask,value in kappa.items():
        if mask==0:continue
        nJ=sum(t=='J' and bool(mask>>i&1) for i,t in enumerate(types))
        coeff[mask]={(1,k):F(value)*comb(nJ,k) for k in range(nJ+1)}
    for i,t in enumerate(types):
        shift={(1,0):F(-1)} if t=='W' else {(1,0):F(-2),(1,1):F(2)}
        coeff[1<<i]=padd(coeff.get(1<<i,{}),shift)
    @lru_cache(None)
    def partition(mask):
        if not mask:return {(0,0):F(1)}
        first=mask&-mask;out={};sub=mask
        while sub:
            if sub&first and sub in coeff:
                out=padd(out,pmul(coeff[sub],partition(mask^sub)))
            sub=(sub-1)&mask
        return out
    return partition((1<<marks)-1)

@lru_cache(None)
def raw_moment(word):
    length,blocks,wicks,types=layout(word)
    if length%2:return {},dict(accepted_diagrams=0,negative_diagram_coefficients=0,maximum_source_order=0)
    out={};diagrams=negative=order=0
    for pa in pairings(tuple(range(length))):
        if any(a<=u<v<b for a,b in wicks for u,v in pa):continue
        kappa=exponent(length,pa,blocks)
        cr=int(-kappa.get(0,0)/4)
        weight=shifted_weight(kappa,types)
        out=padd(out,{(i,j+cr):v for (i,j),v in weight.items()})
        diagrams+=1;negative+=sum(v<0 for v in weight.values())
        order=max(order,max((mask.bit_count() for mask in kappa if mask),default=0))
    return out,dict(accepted_diagrams=diagrams,negative_diagram_coefficients=negative,
                    maximum_source_order=order)

def compressed_moment(word):
    positions={c:tuple(i for i,s in enumerate(word) if s==c) for c in 'HWJ'}
    if any(len(p)%2 for p in positions.values()):return {}
    pair_lists=[tuple(pairings(positions[c])) for c in 'HWJ']
    result=Counter()
    for selections in product(*pair_lists):
        chords=[(ch,2 if c=='J' else 1) for c,pa in zip('HWJ',selections) for ch in pa]
        power=sum(a*b*crossing(x,y) for (x,a),(y,b) in combinations(chords,2))
        result[0,power]+=1
    # J variance numerator is 4lambda(1+q)+16lambda^2 q.
    # Clearing one (1+q) for every J insertion leaves half of them.
    countJ=len(positions['J']);countW=len(positions['W'])
    factor=ppower({(1,0):4,(1,1):4,(2,1):16},countJ//2)
    factor=pmul(factor,ppower({(0,0):1,(0,1):1},countJ//2))
    factor=pmul(factor,{(countW//2,0):1})
    return pmul(dict(result),factor)

def evaluate(poly,q,lam,jcount=0):
    return sum((v*F(lam)**i*F(q)**j for (i,j),v in poly.items()),F(0))/(1+F(q))**jcount

def encode(poly):
    return [dict(lambda_degree=i,q_degree=j,coefficient=str(v)) for (i,j),v in sorted(poly.items())]
