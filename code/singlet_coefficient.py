"""Exact mode-polynomial reduction of the derivative-singlet second moment.

Integer determinants are local (at most 6 x 6). Coefficient arrays use Python
integers via dtype=object. No Hilbert-space matrices enter this implementation.
"""
from functools import lru_cache
from itertools import permutations
from math import comb
from fractions import Fraction
import numpy as np

ALL_MARKS=((((0,1),(4,5)),1),(((0,1),(3,4)),-1),
           (((1,2),(4,5)),-1),(((1,2),(3,4)),1))


def pairings(items):
    if not items:
        yield ()
    else:
        for j in range(1,len(items)):
            for tail in pairings(items[1:j]+items[j+1:]):
                yield ((items[0],items[j]),)+tail


def cross_pairings():
    for perm in permutations((3,4,5)):
        yield tuple(zip((0,1,2),perm))


def determinant(mat):
    n=len(mat)
    if not n:return 1
    if n==1:return mat[0][0]
    return sum((-1)**j*mat[0][j]*determinant([row[:j]+row[j+1:] for row in mat[1:]])
               for j in range(n) if mat[0][j])


@lru_cache(None)
def local_polynomial(length,pairing,marks):
    """Nonconstant terms of det(I+Z (2T) P), in 1,h,g,hg order."""
    pi={a:b for a,b in pairing}|{b:a for a,b in pairing}
    terms=[]
    for mask in range(1,1<<length):
        ss=[i for i in range(length) if mask>>i&1]
        values=[]
        for h,g in [(0,0),(1,0),(0,1),(1,1)]:
            mat=[]
            for a in ss:
                row=[]
                for b in ss:
                    value=1 if a<=pi[b] else -1
                    if (pi[b],a)==marks[0]:value-=2*h
                    if (pi[b],a)==marks[1]:value-=2*g
                    row.append(value)
                mat.append(row)
            values.append(determinant(mat))
        base,h,g,hg=values
        dual=(base,h-base,g-base,hg-h-g+base)
        if any(dual):
            terms.append((tuple((mask>>i)&1 for i in range(length)),dual))
    return tuple(terms)


def power_coefficients(p,length,terms):
    """[z_1^p...z_l^p hg] f^j, j=0,...,length*p, f has no constant."""
    cur=np.zeros((4,)+(p+1,)*length,dtype=object)
    cur[(0,)+(0,)*length]=1
    answer=[0]
    for degree in range(1,length*p+1):
        nxt=np.zeros_like(cur)
        for bits,(base,h,g,hg) in terms:
            src=tuple(slice(0,p+1-z) for z in bits)
            dst=tuple(slice(z,p+1) for z in bits)
            if base:nxt[(slice(None),)+dst]+=base*cur[(slice(None),)+src]
            if h:
                nxt[(1,)+dst]+=h*cur[(0,)+src]
                nxt[(3,)+dst]+=h*cur[(2,)+src]
            if g:
                nxt[(2,)+dst]+=g*cur[(0,)+src]
                nxt[(3,)+dst]+=g*cur[(1,)+src]
            if hg:nxt[(3,)+dst]+=hg*cur[(0,)+src]
        cur=nxt
        answer.append(int(cur[(3,)+(p,)*length]))
    return answer


def moment_coefficients(p,method='cross'):
    """Return c_j with mu2/w^3 = 2^(-6p) sum_j binom(N,j)c_j."""
    total=[0]*(6*p+1)
    if method=='cross':
        configurations=[(pair,marks,sgn) for pair in cross_pairings() for marks,sgn in ALL_MARKS]
    elif method=='all15':
        configurations=[(pair,marks,sgn) for pair in pairings(tuple(range(6)))
            for marks,sgn in [(((2,3),(4,5)),2),(((1,2),(4,5)),-2)]]
    else:raise ValueError(method)
    # Group identical local polynomials before raising them to powers.
    grouped={}
    for pair,marks,sgn in configurations:
        terms=local_polynomial(6,pair,marks)
        grouped[terms]=grouped.get(terms,0)+sgn
    for terms,sgn in grouped.items():
        if sgn:
            row=power_coefficients(p,6,terms)
            total=[a+sgn*b for a,b in zip(total,row)]
    return [(-1)**(3*p)*a for a in total]


def evaluate(n,p,coefficients):
    return Fraction(sum(comb(n,j)*v for j,v in enumerate(coefficients) if j<=n),2**(6*p))


def overlap_table():
    """Analytic Poisson-overlap exponents for all 24 marked cross diagrams."""
    rows=[]
    for pairing in cross_pairings():
        crossings=sum(a<c<b<d or c<a<d<b for i,(a,b) in enumerate(pairing) for c,d in pairing[i+1:])
        for marks,sgn in ALL_MARKS:
            terms=local_polynomial(6,pairing,marks)
            singles={bits.index(1):v for bits,v in terms if sum(bits)==1}
            assert all(v[1:]==(0,0,0) for v in singles.values())
            alpha=beta=delta=0
            # All fifteen pairs contribute -1, including identically zero
            # principal minors omitted from the sparse local polynomial.
            kappa=-15
            for bits,v in terms:
                if sum(bits)==2:
                    inds=[i for i,z in enumerate(bits) if z]
                    den=singles[inds[0]][0]*singles[inds[1]][0]
                    assert abs(den)==1
                    kappa+=v[0]*den
                    alpha+=v[1]*den;beta+=v[2]*den;delta+=v[3]*den
            assert kappa-3==-4*crossings
            rows.append(dict(pairing=pairing,marks=marks,sign=sgn,crossings=crossings,
                             alpha=alpha,beta=beta,delta=delta))
    return rows


def dual_product(a,b):
    return (a[0]*b[0],a[0]*b[1]+a[1]*b[0],a[0]*b[2]+a[2]*b[0],
            a[0]*b[3]+a[1]*b[2]+a[2]*b[1]+a[3]*b[0])


def fixed_word(n,p,creator_sets,pairing,marks):
    terms=dict(local_polynomial(6,pairing,marks))
    value=(1,0,0,0)
    for mode in range(n):
        bits=tuple(int(mode in subset) for subset in creator_sets)
        local=(1,0,0,0) if not any(bits) else terms.get(bits,(0,0,0,0))
        value=dual_product(value,local)
    return tuple(Fraction((-1)**(3*p)*v,2**(6*p)) for v in value)
