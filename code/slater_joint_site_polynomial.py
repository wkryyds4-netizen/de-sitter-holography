"""Exact local CAR polynomial for a Gaussian pairing in a fixed Slater state.

Every edge is oriented with right endpoint X=P^dag Q and left X^dag.
G: unrestricted; F: P in E,Q in O; B: P in O,Q in E.
No scientific main or file mutation on import.
"""
from itertools import combinations
from collections import defaultdict

def pairings(slots):
    if not slots:yield ();return
    a=slots[0]
    for j in range(1,len(slots)):
        for rest in pairings(slots[1:j]+slots[j+1:]):
            yield ((a,slots[j]),)+rest

def crossed(a,b):
    i,j=a;k,l=b
    return i<k<j<l or k<i<l<j

def local_coefficient(pi,kinds,occupied,mask):
    """Coefficient of a single-site membership mask, including JW sign."""
    k=len(pi)
    for e,kind in enumerate(kinds):
        P=bool(mask&(1<<(2*e)));Q=bool(mask&(1<<(2*e+1)))
        if kind=='F' and ((occupied and P) or (not occupied and Q)):return 0
        if kind=='B' and ((occupied and Q) or (not occupied and P)):return 0
    slots={}
    for e,(i,j) in enumerate(pi):slots[i]=(e,False);slots[j]=(e,True)
    state=occupied
    for slot in range(2*k-1,-1,-1):
        e,right=slots[slot]
        annih,create=(2*e+1,2*e) if right else (2*e,2*e+1)
        if mask&(1<<annih):
            if not state:return 0
            state=0
        if mask&(1<<create):
            if state:return 0
            state=1
    assert state==occupied
    parity=sum(((mask>>(2*i))&3).bit_count()*((mask>>(2*j))&3).bit_count()
        for i,j in combinations(range(k),2) if crossed(pi[i],pi[j]))
    return -1 if parity%2 else 1

def polynomial(pi,kinds,occupied):
    return {mask:v for mask in range(1<<(2*len(pi)))
            if (v:=local_coefficient(pi,kinds,occupied,mask))}

def coefficient(n,p,pi,kinds):
    """Exact coefficient of product z_i^p in F_E^(n/2) F_O^(n/2)."""
    dim=2*len(pi);base=p+1
    powers=[base**i for i in range(dim)]
    state={0:1}
    for occupied in (0,1):
        poly=polynomial(pi,kinds,occupied)
        terms=[(sum(powers[i] for i in range(dim) if mask>>i&1),
                tuple(i for i in range(dim) if mask>>i&1),v) for mask,v in poly.items()]
        for _ in range(n//2):
            result=defaultdict(int)
            for x,c in state.items():
                digits=[x//b%base for b in powers]
                for dx,indices,v in terms:
                    if all(digits[i]<p for i in indices):result[x+dx]+=c*v
            state={x:c for x,c in result.items() if c}
    return state.get(sum(p*powers[i] for i in range(dim)),0)

def penalty(pi,kinds):
    """Proposed/proved-in-companion-text all-forward pair penalty count."""
    total=0
    for i,j in combinations(range(len(pi)),2):
        a,b=pi[i];c,d=pi[j]
        if b<c or d<a:continue
        if crossed(pi[i],pi[j]):total+=1;continue
        inner=j if a<c<d<b else i
        total+=kinds[inner]=='F'
    return total
