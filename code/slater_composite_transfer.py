"""Pure finite-memory transfer for Slater H/W/J/+/- moments.

W=sqrt(lambda)u, J=j. States are auxiliary weighted automaton states,
NOT an orthonormal Hilbert basis. No scientific driver runs on import.
"""
from collections import defaultdict, Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations

EMPTY = ('', ())

def geometric_exponent(pi, kinds, blocks):
    """Geometric local-source classification, independently of CAR matrices."""
    from slater_joint_site_polynomial import crossed, penalty
    edge={x:i for i,p in enumerate(pi) for x in p}
    mate={a:b for a,b in pi}|{b:a for a,b in pi}
    out={0:-4*penalty(pi,kinds)}
    for t,(a,b,kind) in enumerate(blocks):
        if kind=='L':v=int(kinds[edge[a]]=='G')
        else:
            assert b==a+2 and edge[a]!=edge[a+1]
            e,f=edge[a],edge[a+1];v=2-4*int(crossed(pi[e],pi[f]))
            if mate[a]>a+1 and mate[a+1]>a+1:
                younger=min((e,f),key=lambda x:pi[x][1])
                v*=int(kinds[younger]=='G')
            elif mate[a]<a and mate[a+1]<a:
                first=max((e,f),key=lambda x:pi[x][0])
                v*=int(kinds[first]=='G')
        out[1<<t]=v
    for t,s in combinations(range(len(blocks)),2):
        a,b,k=blocks[t];c,d,l=blocks[s];v=0
        if k==l=='L':v=int(edge[a]==edge[c] and kinds[edge[a]]=='G')
        elif k==l=='D':v=4*int(mate[a]==c+1 and mate[a+1]==c)
        out[(1<<t)|(1<<s)]=v
    return {i:v for i,v in out.items() if v}

def _append(state, color):
    w,links=state
    return w+color,links

def _link(state, mode):
    w,links=state;n=len(w)
    return w,tuple(sorted(links+((n-2,n-1,mode),)))

def _delete(state,j,lam,kap,allow_pair=False):
    w,links=state;factor=1;new=[]
    for a,b,mode in links:
        if j in (a,b):
            if mode=='R':
                if not allow_pair:return None,0
            else:
                factor*=2*lam*(1 if j==b else -1)*int(w[j]=='H')-kap
        else:new.append((a-int(a>j),b-int(b>j),mode))
    return (w[:j]+w[j+1:],tuple(new)),factor

def _opens(state,op):
    for c in {'H':'Hh','W':'Uw','+':'P','-':''}[op]:
        yield _append(state,c)

def _close(state,op,j,gamma,lam,kap):
    w,_=state;c=w[j]
    weight={'H':{'H':1,'P':gamma},'W':{'U':lam,'P':-lam*gamma},
            '-':{'h':gamma,'w':-lam*gamma,'P':1},'+':{}}[op].get(c,0)
    if not weight:return None,0
    new,factor=_delete(state,j,lam,kap)
    if not factor:return None,0
    exponent=len(w)-1-j if c in 'HU' else len(w)-1
    return new,weight*factor*gamma**(4*exponent)

@lru_cache(None)
def transitions(state,op,gamma,lam):
    q=gamma**4;kap=2*lam*(1-q)/(1+q);w,links=state;n=len(w)
    out=defaultdict(lambda:F(0))
    if op!='J':
        for st in _opens(state,op):out[st]+=1
        for j in range(n):
            st,c=_close(state,op,j,gamma,lam,kap)
            if c:out[st]+=c
    else:
        # Two singleton-assigned H endpoints, scanned right then left.
        for first in _opens(state,'H'):
            for second in _opens(first,'H'):out[_link(second,'S')]+=1
            for j in range(n): # never close the new endpoint internally
                st,c=_close(first,'H',j,gamma,lam,kap)
                if c:out[st]+=c*(-2*lam-kap)
        for j in range(n):
            first,c=_close(state,'H',j,gamma,lam,kap)
            if not c:continue
            for st in _opens(first,'H'):out[st]+=c*(2*lam-kap)
            for k in range(n-1):
                original=k+int(k>=j)
                st,d=_close(first,'H',k,gamma,lam,kap)
                if d:
                    alpha=(2 if j>original else -2)*int(w[max(j,original)]=='H')
                    out[st]+=c*d*(lam*alpha-kap)
        # A source pair is exactly two nested G chords joining two J blocks.
        out[_link(_append(_append(state,'H'),'H'),'R')]+=1
        for a,b,mode in links:
            if mode!='R':continue
            first,_=_delete(state,b,lam,kap,allow_pair=True)
            st,_=_delete(first,a,lam,kap)
            out[st]+=4*lam*q**((n-1-b)+(n-2-a))
    return tuple((s,c) for s,c in out.items() if c)

def step(vector,op,gamma,lam):
    gamma,lam=F(gamma),F(lam);out=defaultdict(lambda:F(0))
    for state,c in vector.items():
        for target,d in transitions(state,op,gamma,lam):out[target]+=c*d
    return {s:c for s,c in out.items() if c}

def vector(word,gamma,lam,max_open=None):
    v={EMPTY:F(1)}
    for op in reversed(word):
        v=step(v,op,gamma,lam)
        if max_open is not None:v={s:c for s,c in v.items() if len(s[0])<=max_open}
    return v

def moment(word,gamma,lam,stats=None):
    """Exact vacuum coefficient; discard only states impossible to close."""
    v={EMPTY:F(1)};remaining=sum(2 if c=='J' else 1 for c in word);peak=1
    left=Counter(word)
    for op in reversed(word):
        left[op]-=1
        remaining-=2 if op=='J' else 1
        v=step(v,op,gamma,lam)
        v={s:c for s,c in v.items() if len(s[0])<=remaining
           and s[0].count('h')+s[0].count('w')<=left['-']
           and s[0].count('U')<=left['W']
           and s[0].count('H')<=left['H']+2*left['J']
           and sum(mode=='R' for _,_,mode in s[1])<=left['J']}
        peak=max(peak,len(v))
    if stats is not None:stats.update(peak_states=peak,expanded_endpoints=sum(2 if c=='J' else 1 for c in word))
    return v.get(EMPTY,F(0))

def adjoint(word):return ''.join({'+':'-','-':'+'}.get(c,c) for c in reversed(word))

def gram(words,gamma,lam):
    return [[moment(adjoint(a)+b,gamma,lam) for b in words] for a in words]
