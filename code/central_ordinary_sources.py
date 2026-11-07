"""Pure central-sector ordinary marked sources and same-copy regression.

No finite-N generator is executed on import. The central exponent is the
Fourier-conditioned full-trace exponent, not a restriction of a trace norm.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import comb,prod
from memory_body_source_chords import exponent,G_exponent
from memory_joint_regression import JointModel,Term

def internal_counts(pairing,blocks):
    return tuple(sum(a<=x<y<b for x,y in pairing) if kind=='D' else 0
                 for a,b,kind in blocks)

def central_exponent(length,pairing,blocks,gaussian_number_shift=False):
    base=G_exponent(length,pairing,blocks) if gaussian_number_shift else exponent(length,pairing,blocks)
    out=Counter(base);inside=internal_counts(pairing,blocks)
    for a,b in combinations(range(len(blocks)),2):
        out[(1<<a)|(1<<b)]-=4*inside[a]*inside[b]
    return {m:v for m,v in out.items() if v}

class CentralModel(JointModel):
    """Fixed-charge X=0 moments; inherited terms retain all chaos/mark regions."""
    def tilted_coefficient(self,linears,internals,connected,charge_power):
        if charge_power:return F(0)
        if not linears:return F(1)
        if len(linears)==1:return self.lam*linears[0]
        if len(linears)==2:
            return self.lam*(connected-4*internals[0]*internals[1])+self.lam**2*prod(linears)
        raise ValueError('Two-field Gram evaluator accepts at most two marks')

def central_singleton_coefficient(N,p,chords,internals):
    """Exact half-charge coefficient for disjoint singleton supports.

    Divided by the source-free coefficient. This isolates conditioning;
    it does not omit overlap corrections from a full model calculation.
    """
    remaining=N-2*chords*p;order=len(internals)
    assert remaining>=2*order and remaining%2==0
    center=remaining//2
    return F(p**order*prod(internals),comb(remaining,center))*sum(
        (-1)**(order-j)*comb(order,j)*comb(remaining,center-j)
        for j in range(order+1))

def limiting_singleton_coefficient(lam,internals):
    order=len(internals)
    if order%2:return F(0)
    return (-4*F(lam))**(order//2)*prod(range(1,order,2))*prod(internals)

def central_G2_zero_chaos(N,p):
    """E G2 / C_m at m=N/2, from the exact adjacent-sector covariance."""
    assert N%2==0 and p<=N//2
    return -F(p*(p+1),N//2+1)
