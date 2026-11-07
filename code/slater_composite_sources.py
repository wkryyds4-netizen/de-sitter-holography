"""Pure same-J Slater H/W/J/M source calculus; no file writes.

H=h, W=(p-L-lambda)h, J=(2p-D-kappa)Wick2(h*h).
Thus W=sqrt(lambda)u and J is the dimensionless centered curvature.
Output polynomial is in lambda,gamma with denominator (1+gamma^4)^#J.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from slater_joint_site_polynomial import pairings,crossed,penalty
from memory_body_source_chords import mm,loss,mul,inverse,I2,C2,CD2,ZERO
from bilocal_multisource_chords import shifted_weight

def layout(word):
    slots=[];blocks=[];types=[]
    for letter in word:
        a=len(slots)
        if letter=='J':slots.extend('HH');blocks.append((a,a+2,'D'));types.append('J')
        elif letter=='W':slots.append('H');blocks.append((a,a+1,'L'));types.append('W')
        else:
            assert letter in 'H+-';slots.append(letter)
    return tuple(slots),tuple(blocks),tuple(types)

@lru_cache(None)
def local(pairing,kinds,blocks,occupied,mask):
    for e,kind in enumerate(kinds):
        P=bool(mask>>(2*e)&1);Q=bool(mask>>(2*e+1)&1)
        if kind=='F' and ((occupied and P) or (not occupied and Q)):return {}
        if kind=='B' and ((occupied and Q) or (not occupied and P)):return {}
    endpoints={}
    for e,(i,j) in enumerate(pairing):endpoints[i]=(2*e+1,2*e);endpoints[j]=(2*e,2*e+1)
    starts={a:(k,b,kind) for k,(a,b,kind) in enumerate(blocks)}
    result={0:I2};pos=0
    while pos<2*len(pairing):
        mark,end,kind=starts.get(pos,(None,pos+1,None));value=I2
        for slot in range(pos,end):
            create,annih=endpoints[slot]
            if mask>>create&1:value=mm(value,CD2)
            if mask>>annih&1:value=mm(value,C2)
        factor={0:value}
        if mark is not None:factor[1<<mark]=tuple(-x for x in loss(value,kind))
        result={i|j:mm(a,b) for i,a in result.items() for j,b in factor.items()}
        pos=end
    sign=(-1)**sum(((mask>>(2*i))&3).bit_count()*((mask>>(2*j))&3).bit_count()
                   for i,j in combinations(range(len(pairing)),2) if crossed(pairing[i],pairing[j]))
    entry=3 if occupied else 0
    return {i:sign*a[entry] for i,a in result.items() if a[entry]}

@lru_cache(None)
def source_exponent(pairing,kinds,blocks):
    """Source coefficients after exact singleton cancellation; no internal pairs."""
    assert not any(a<=u<v<b for a,b,_ in blocks for u,v in pairing)
    singles=[local(pairing,kinds,blocks,i%2,1<<i) for i in range(2*len(pairing))]
    assert all(x.get(0)==1 for x in singles)
    answer=Counter()
    for i,j in combinations(range(2*len(pairing)),2):
        inv=inverse(mul(singles[i],singles[j]),len(blocks))
        for occupied in (0,1):
            term=local(pairing,kinds,blocks,occupied,(1<<i)|(1<<j))
            answer.update({mask:2*v for mask,v in mul(term,inv).items()})
    # The source-free exponent separately contains the falling-factorial
    # and H/M normalization. Return it for comparison with old all-G laws.
    k=len(pairing);nG=kinds.count('G')
    answer[0]-=2*k*(k-1)+2*nG
    assert answer[0]==-4*penalty(pairing,kinds)
    return {i:v for i,v in answer.items() if v}

@lru_cache(None)
def raw_moment(word):
    slots,blocks,types=layout(word);out=Counter();count=0
    if len(slots)%2:return {},0
    for pi in pairings(tuple(range(len(slots)))):
        if any(a<=i<j<b for a,b,_ in blocks for i,j in pi):continue
        kinds=[];mixed=0;valid=True
        for i,j in pi:
            left,right=slots[i],slots[j]
            if left==right=='H':kinds.append('G')
            elif left=='+' or right=='-':valid=False;break
            else:
                kinds.append('F');mixed+=int(left=='H' or right=='H')
        if not valid:continue
        kinds=tuple(kinds);kap=source_exponent(pi,kinds,blocks)
        weight=shifted_weight(kap,types)
        for (ld,qd),v in weight.items():out[ld,mixed+4*penalty(pi,kinds)+4*qd]+=v
        count+=1
    return {key:v for key,v in out.items() if v},count

def symbolic(word,gamma,lam):
    return sum(v*lam**i*gamma**j for (i,j),v in raw_moment(word)[0].items())/(1+gamma**4)**word.count('J')

def evaluate(word,gamma,lam):
    gamma,lam=F(gamma),F(lam)
    return sum((v*lam**i*gamma**j for (i,j),v in raw_moment(word)[0].items()),F(0))/(1+gamma**4)**word.count('J')
