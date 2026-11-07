"""Pure fixed-degree marked-Wick calculus; never writes files on import.

Slots are normalized H factors. Disjoint D/L marks act before sector
restriction. A Wick region excludes internal Gaussian pairs; a chaos
region (a,b,d) retains exactly d edges exiting that region. Regions may
contain smaller marks. Center shifts are supplied separately in slot order.
"""
from collections import Counter
from fractions import Fraction as F
from slater_composite_sources import source_exponent,local
from slater_joint_site_polynomial import pairings,penalty
from memory_body_source_chords import mul,inverse

def validate(length,marks,wick_regions=(),chaos_regions=()):
    assert length>=0
    assert all(0<=a<b<=length and t in ('D','L') for a,b,t in marks)
    assert all(marks[i][1]<=marks[i+1][0] for i in range(len(marks)-1))
    assert all(0<=a<b<=length for a,b in wick_regions)
    assert all(0<=a<b<=length and 0<=d<=b-a and (b-a-d)%2==0
               for a,b,d in chaos_regions)

def allowed(pi,marks,wick_regions=(),chaos_regions=()):
    if any(a<=i<j<b for a,b in tuple(wick_regions)+tuple((a,b) for a,b,_ in marks)
           for i,j in pi):return False
    return all(sum(a<=i<j<b for i,j in pi)==(b-a-d)//2 for a,b,d in chaos_regions)

def exponential(kappa,marks,lam,shifts=None):
    """exp(lambda*K_nonconstant - sum shift_i z_i), all z_i^2=0."""
    shifts=(F(0),)*marks if shifts is None else tuple(map(F,shifts))
    assert len(shifts)==marks
    p=Counter({m:F(lam)*v for m,v in kappa.items() if m})
    for i,s in enumerate(shifts):p[1<<i]-=s
    p={m:v for m,v in p.items() if v}
    out=Counter({0:F(1)});power={0:F(1)}
    for degree in range(1,marks+1):
        power={m:v/F(degree) for m,v in mul(power,p).items()}
        out.update(power)
    return {m:v for m,v in out.items() if v}

def collision_part(pi,kinds,marks,edge_pair):
    """Source-dependent own-edge or two-edge exponent; constant omitted."""
    e,f=edge_pair
    pairs=[(2*e,2*e+1)] if e==f else [(i,j) for i in (2*e,2*e+1) for j in (2*f,2*f+1)]
    out=Counter()
    for i,j in pairs:
        si=local(pi,kinds,marks,i%2,1<<i)
        sj=local(pi,kinds,marks,j%2,1<<j)
        inv=inverse(mul(si,sj),len(marks))
        for occupied in (0,1):
            out.update({m:2*v for m,v in mul(local(pi,kinds,marks,occupied,(1<<i)|(1<<j)),inv).items()})
    return {m:v for m,v in out.items() if m and v}

def moment(length,marks,gamma,lam,wick_regions=(),chaos_regions=(),shifts=None,mode='trace'):
    """Exact rational moment: trace, selected=M B Mdag, or chi= (h-W)B(h-W).

    This is a diagram evaluator, not a finite-N convergence/tail theorem.
    """
    marks,wick_regions,chaos_regions=tuple(marks),tuple(wick_regions),tuple(chaos_regions)
    validate(length,marks,wick_regions,chaos_regions)
    assert mode in ('trace','selected','chi')
    gamma,lam=F(gamma),F(lam)
    shifts=(F(0),)*len(marks) if shifts is None else tuple(map(F,shifts))
    if length%2:return F(0)
    external=mode!='trace';offset=int(external);n=length+2*offset
    shifted=tuple((a+offset,b+offset,t) for a,b,t in marks)
    wr=tuple((a+offset,b+offset) for a,b in wick_regions)
    cr=tuple((a+offset,b+offset,d) for a,b,d in chaos_regions)
    if mode=='chi':
        blocks=((0,1,'L'),)+shifted+((n-1,n,'L'),)
        centers=(lam,)+shifts+(lam,)
    else:blocks=shifted;centers=shifts
    total=F(0)
    for pi in pairings(tuple(range(n))):
        if not allowed(pi,shifted,wr,cr):continue
        kinds=tuple('F' if mode=='selected' and (a==0 or b==n-1) else 'G' for a,b in pi)
        mixed=2 if mode=='selected' and (0,n-1) not in pi else 0
        coeff=exponential(source_exponent(pi,kinds,blocks),len(blocks),lam,centers)
        if mode=='chi':
            internal=((1<<len(marks))-1)<<1;right=1<<(len(blocks)-1)
            weight=sum(sign*coeff.get(mask,0) for mask,sign in
                       [(internal,1),(internal|1,-1),(internal|right,-1),(internal|1|right,1)])
        else:weight=coeff.get((1<<len(marks))-1,0)
        total+=gamma**(mixed+4*penalty(pi,kinds))*weight
    return total
