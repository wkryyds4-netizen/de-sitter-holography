"""Multi-mark Poisson minors and independently paired mixed-probe words.

Sparse nilpotent variables are bit masks; q,lambda polynomials use integer
coefficient maps (lambda degree, q degree). No floating point in identities.
"""
from collections import Counter
from itertools import combinations
from math import comb
from singlet_coefficient import pairings


def crossing(a,b):
    x,y=a;z,w=b
    return int(x<z<y<w or z<x<w<y)


def nmul(a,b):
    out=Counter()
    for i,x in a.items():
        for j,y in b.items():
            if not i&j:out[i|j]+=x*y
    return {k:v for k,v in out.items() if v}


def expand(word):
    marks=[];length=0
    for c in word:
        if c=='X':marks.append((length,length+1));length+=2
        else:assert c=='H';length+=1
    return length,marks


def local(length,pairing,marks):
    pi={a:b for a,b in pairing}|{b:a for a,b in pairing}
    lookup={m:1<<j for j,m in enumerate(marks)}
    def entry(a,b):
        out={0:1 if a<=pi[b] else -1}
        mask=lookup.get((pi[b],a))
        if mask is not None:out[mask]=-2
        return out
    singles=[entry(i,i) for i in range(length)]
    assert all(len(s)==1 for s in singles)
    k=Counter({0:-comb(length,2)-length//2})
    for i,j in combinations(range(length),2):
        den=singles[i][0]*singles[j][0]
        for mask,v in nmul(singles[i],singles[j]).items():k[mask]+=den*v
        for mask,v in nmul(entry(i,j),entry(j,i)).items():k[mask]-=den*v
    cr=sum(crossing(a,b) for a,b in combinations(pairing,2))
    assert k[0]==-4*cr
    assert all(mask.bit_count()<=2 for mask,v in k.items() if v)
    alpha=[k[1<<i] for i in range(len(marks))]
    delta={(i,j):k[(1<<i)|(1<<j)] for i,j in combinations(range(len(marks)),2)}
    # Independent geometrical rule check against the determinant minors.
    for i,(a,b) in enumerate(marks):
        assert alpha[i]==2-4*crossing(tuple(sorted((a,pi[a]))),tuple(sorted((b,pi[b]))))
    for (i,j),v in delta.items():
        a,b=marks[i];c,d=marks[j]
        assert v==4*int(pi[a]==d and pi[b]==c)
    return cr,alpha,delta


def pmul(a,b):
    out=Counter()
    for (i,j),x in a.items():
        for (k,l),y in b.items():out[i+k,j+l]+=x*y
    return {k:v for k,v in out.items() if v}


def ppower(a,n):
    out={(0,0):1}
    for _ in range(n):out=pmul(out,a)
    return out


def mark_weight(alpha,delta):
    """(1+q)^m times the shifted derivative of exp(lambda*kappa)."""
    def rec(indices):
        if not indices:return {(0,0):1}
        i,*tail=indices
        # lambda*(alpha_i-2(1-q)/(1+q)) with denominator cleared.
        linear={(1,0):alpha[i]-2,(1,1):alpha[i]+2}
        out=Counter(pmul(linear,rec(tail)))
        for j in tail:
            v=delta.get((i,j),0)
            if v:
                pair={(1,0):v,(1,1):2*v,(1,2):v}
                out.update(pmul(pair,rec([x for x in tail if x!=j])))
        return {k:v for k,v in out.items() if v}
    return rec(list(range(len(alpha))))


def raw_moment(word):
    length,marks=expand(word);out=Counter();count=0
    if length%2:return {},0
    for p in pairings(tuple(range(length))):
        if any(m in p for m in marks):continue
        cr,alpha,delta=local(length,p,marks);count+=1
        for (i,j),v in mark_weight(alpha,delta).items():out[i,j+cr]+=v
    return {k:v for k,v in out.items() if v},count


def probe_moment(word):
    """Pair compressed H and X independently, weights q^(length product)."""
    hs=tuple(i for i,c in enumerate(word) if c=='H')
    xs=tuple(i for i,c in enumerate(word) if c=='X')
    if len(hs)%2 or len(xs)%2:return {}
    out=Counter()
    for a in pairings(hs):
        for b in pairings(xs):
            chords=[(p,1) for p in a]+[(p,2) for p in b]
            degree=sum(crossing(p,r)*i*j for (p,i),(r,j) in combinations(chords,2))
            out[0,degree]+=1
    return dict(out)


def predicted(word):
    m=word.count('X')
    if m%2:return {}
    # [(1+q)^2 W_d/sigma^4]^(m/2).
    pair={(1,0):4,(1,1):8,(1,2):4,(2,1):16,(2,2):16}
    return pmul(ppower(pair,m//2),probe_moment(word))
