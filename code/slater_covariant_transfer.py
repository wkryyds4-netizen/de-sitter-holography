"""Pure H/U/P transfer and marked-site utilities; no file mutations.

U=(p-half-Clifford-size-lambda)h/sqrt(lambda).
Parameters gamma=q**(1/4), d=sqrt(lambda); exact tests may vary them
independently, whereas the microscopic curve imposes q=exp(-4*d*d).
"""
from collections import defaultdict
from functools import lru_cache
from itertools import combinations
from slater_joint_site_polynomial import crossed,pairings,penalty

def marked_site(pi,occupied,membership):
    """Nilpotent coefficients for sources 2(p-L) at every endpoint.

Replacing n by I contributes 1 to the doubled-loss mark. All Jordan
signs are retained from the unmarked label regrouping.
"""
    slots={}
    for e,(i,j) in enumerate(pi):slots[i]=(e,False);slots[j]=(e,True)
    state={(occupied,0):1}
    for slot in range(2*len(pi)-1,-1,-1):
        e,right=slots[slot];P=bool(membership&(1<<(2*e)));Q=bool(membership&(1<<(2*e+1)))
        annih,create=(Q,P) if right else (P,Q)
        out=defaultdict(int)
        for (n,mask),v in state.items():
            if P and Q:out[n,mask|(1<<slot)]+=v
            nn=n
            if annih:
                if not nn:continue
                nn=0
            if create:
                if nn:continue
                nn=1
            out[nn,mask]+=v
        state=out
    sign=(-1)**sum(((membership>>(2*i))&3).bit_count()*((membership>>(2*j))&3).bit_count()
                   for i,j in combinations(range(len(pi)),2) if crossed(pi[i],pi[j]))
    return {mask:sign*v for (n,mask),v in state.items() if n==occupied and v}

def poly_mul(a,b):
    out=defaultdict(int)
    for i,x in a.items():
        for j,y in b.items():
            if not i&j:out[i|j]+=x*y
    return {i:x for i,x in out.items() if x}

def clean(v):return {w:c for w,c in v.items() if c}
def add(*vs):
    out=defaultdict(int)
    for v in vs:
        for w,c in v.items():out[w]+=c
    return clean(out)
def scale(c,v):return clean({w:c*x for w,x in v.items()})

def step(v,op,gamma,d):
    q=gamma**4;c={'H':gamma,'U':-d*gamma,'P':1};out=defaultdict(int)
    for w,x in v.items():
        n=len(w)
        if op in ('H','U'):out[w+op]+=x
        if op in ('A','V','P'):out[w+{'A':'H','V':'U','P':'P'}[op]]+=x
        if op in ('H','U','a','v','m') and n:
            color='H' if op in ('H','a') else 'U'
            for j,t in enumerate(w):
                weight=c[t]*q**(n-1) if op=='m' else (
                    q**(n-1-j) if t==color else (c[color]*q**(n-1) if t=='P' else 0))
                if weight:out[w[:j]+w[j+1:]]+=x*weight
    return clean(out)

def vector(word,gamma,d):
    v={'':1}
    for op in reversed(word):v=step(v,{'H':'H','U':'U','+':'P','-':'m'}[op],gamma,d)
    return v
def moment(word,gamma,d):return vector(word,gamma,d).get('',0)

@lru_cache(None)
def basis_inner(u,v,gamma,d):
    if len(u)!=len(v):return 0
    state={v:1}
    for letter in reversed(u):state=step(state,{'H':'a','U':'v','P':'m'}[letter],gamma,d)
    return state.get('',0)
def inner(u,v,gamma,d):
    return sum(x*y*basis_inner(a,b,gamma,d) for a,x in u.items() for b,y in v.items())

def paired_moment(word,gamma,d):
    total=0;q=gamma**4
    for pi in pairings(tuple(range(len(word)))):
        coefficient=1;kinds=[]
        for i,j in pi:
            a,b=word[i],word[j]
            if a in 'HU' and b in 'HU':
                if a!=b:coefficient=0;break
                kinds.append('G')
            else:
                if a=='+' or b=='-':coefficient=0;break
                kinds.append('F')
                if a in 'HU' or b in 'HU':
                    coefficient*=gamma*(-d if a=='U' or b=='U' else 1)
        if coefficient:total+=coefficient*q**penalty(pi,kinds)
    return total
