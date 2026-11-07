"""Joint charge/energy moments of whole-field G and neutral descendants.

C=1. q and lambda are independent positive formal parameters in the checks.
Fields are linear combinations of ordered terms. Wick exclusions, physical
mark intervals, number-filter intervals, and charge powers stay separate.
"""
from collections import Counter
from dataclasses import dataclass, replace
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import prod


@dataclass(frozen=True)
class Term:
    length: int
    wicks: tuple = ()
    marks: tuple = ()
    filters: tuple = ()  # (start,end,OU power,selected chaos), -1 means no selection
    xpower: int = 0

    def dagger(self):
        def reflect(items):
            return tuple(sorted((self.length-b,self.length-a,*tail) for a,b,*tail in items))
        return Term(self.length,reflect(self.wicks),reflect(self.marks),reflect(self.filters),self.xpower)


def pairings(slots):
    if not slots:
        yield ()
    else:
        for i in range(1,len(slots)):
            for tail in pairings(slots[1:i]+slots[i+1:]):
                yield ((slots[0],slots[i]),)+tail


def crosses(first,second):
    a,b=sorted(first);c,d=sorted(second)
    return a<c<b<d or c<a<d<b


def cr(chords):
    return sum(crosses(a,b) for a,b in combinations(chords,2))


def region(pa,a,b):
    internal=tuple(ch for ch in pa if all(a<=x<b for x in ch))
    opened=frozenset(ch for ch in pa if sum(a<=x<b for x in ch)==1)
    return internal,opened


class JointModel:
    def __init__(self,q,lam):
        self.q,self.lam=F(q),F(lam)
        assert 0<self.q<1 and self.lam>0
        self.stats=Counter()

    def add(self,*fields):
        out=Counter()
        for field in fields:out.update(field)
        return {term:F(value) for term,value in out.items() if value}

    def scale(self,field,value):
        return {term:F(value)*coefficient for term,coefficient in field.items() if value*coefficient}

    def H(self,n=0):
        return {Term(n):F(1)}

    def X(self,power=1):
        return {Term(0,xpower=power):F(1)}

    def charge(self,field,power=1):
        return {replace(term,xpower=term.xpower+power):value for term,value in field.items()}

    def hleft(self,field):
        def shift(items):return tuple((a+1,b+1,*tail) for a,b,*tail in items)
        return {Term(t.length+1,shift(t.wicks),shift(t.marks),shift(t.filters),t.xpower):v for t,v in field.items()}

    def hright(self,field):
        return {replace(t,length=t.length+1):v for t,v in field.items()}

    def qinteger(self,n):
        return sum((self.q**r for r in range(n)),F(0))

    def P(self,n):
        return prod((self.qinteger(k) for k in range(1,n+1)),start=F(1))

    def logarithmic_q_derivative(self,n):
        return sum((sum((r*self.q**r for r in range(k)),F(0))/self.qinteger(k)
                    for k in range(1,n+1)),F(0))

    @lru_cache(None)
    def hermite_coefficients(self,n):
        if n==0:return ((0,F(1)),)
        if n==1:return ((1,F(1)),)
        out=Counter({degree+1:value for degree,value in self.hermite_coefficients(n-1)})
        for degree,value in self.hermite_coefficients(n-2):out[degree]-=self.qinteger(n-1)*value
        return tuple(sorted((degree,value) for degree,value in out.items() if value))

    def Hermite(self,n):
        return {Term(degree):value for degree,value in self.hermite_coefficients(n)}

    @lru_cache(None)
    def a(self,n):
        if n==0:return ((0,F(1)),)
        out=Counter()
        for d,value in self.a(n-1):
            out[d+1]+=value
            if d:out[d-1]+=value*self.qinteger(d)
        return tuple(sorted(out.items()))

    def G(self,n):
        return {} if n<2 else {Term(n,marks=((0,n,'G'),)):F(1)}

    def g(self,n):
        answer={}
        for d,a in self.a(n):
            internal=(n-d)//2
            h=self.Hermite(d)
            constant=-4*self.lam*self.logarithmic_q_derivative(d)-2*self.lam*internal
            answer=self.add(answer,self.scale(h,a*constant),self.scale(self.charge(h),a*4*self.lam*internal))
        return answer

    def J(self,d):
        if d<2:return {}
        pure={Term(d,((0,d),),((0,d,'D'),)):F(1)}
        kappa=self.lam*(d*(d-1)-4*self.logarithmic_q_derivative(d))
        return self.add(pure,self.scale(self.Hermite(d),-kappa))

    def residual(self,n):
        return self.add(self.G(n),self.scale(self.g(n),-1))

    def resummed(self,n):
        return self.add(*(self.scale(self.J(d),a) for d,a in self.a(n)))

    def neutral(self,field):
        result={}
        for t,poly in field.terms.items():
            filters=((0,t.length,t.ou_power,t.selected_chaos),) if t.ou_power or t.selected_chaos>=0 else ()
            term=Term(t.length,t.wick_blocks,tuple((a,b,'L') for a,b in t.mark_blocks),filters)
            value=sum((coefficient*self.lam**i*self.q**j for (i,j),coefficient in poly.items()),F(0))
            result=self.add(result,{term:value})
        return result

    def gaussian_moment(self,n):
        if n<0 or n%2:return F(0)
        return prod(range(1,n,2),start=F(1))/(4*self.lam)**(n//2)

    def tilted_coefficient(self,linears,internals,connected,charge_power):
        # Derivatives in t of exp(t^2/(8 lambda)) times the mark polynomial.
        count=len(linears)
        if count==0:abc=(F(1),F(0),F(0))
        elif count==1:abc=(self.lam*linears[0],F(internals[0]),F(0))
        elif count==2:
            u,v=linears;i,j=internals
            abc=(self.lam*connected+self.lam**2*u*v,
                 self.lam*(u*j+v*i),F(i*j))
        else:raise ValueError('At most two marked insertions are supported')
        r=charge_power
        return (abc[0]*self.gaussian_moment(r)+abc[1]*r*self.gaussian_moment(r-1)
                +abc[2]*r*(r-1)*self.gaussian_moment(r-2))

    @lru_cache(None)
    def term_moment(self,left,right,a=0,b=0,no_bridge=False):
        dual=left.dagger();offset=dual.length+a;length=offset+right.length+b
        if length%2:return F(0)
        shift=lambda items:tuple((x+offset,y+offset,*tail) for x,y,*tail in items)
        wicks=dual.wicks+shift(right.wicks)
        marks=dual.marks+shift(right.marks)
        filters=dual.filters+shift(right.filters)
        assert len(marks)<=2
        total=F(0)
        for pa in pairings(tuple(range(length))):
            if any(x<=u<y and x<=v<y for x,y in wicks for u,v in pa):continue
            if no_bridge and any(u<dual.length and offset<=v<offset+right.length for u,v in pa):continue
            ou_weight=F(1);permitted=True
            for x,y,power,selected in filters:
                degree=len(region(pa,x,y)[1])
                if selected>=0 and degree!=selected:permitted=False;break
                ou_weight*=self.q**(power*degree)
            if not permitted:continue
            linears=[];internals=[];incidents=[];types=[]
            for x,y,kind in marks:
                inside,opened=region(pa,x,y)
                if kind=='L':
                    assert not inside
                    incident=frozenset(ch for ch in pa if any(x<=v<y for v in ch))
                    linear=(y-x)**2-4*cr(incident);internal=0
                else:
                    incident=opened;internal=len(inside);d=len(opened)
                    linear=-4*cr(opened)-2*internal
                    if kind=='D':linear+=d*(d-1)
                    else:assert kind=='G'
                linears.append(linear);internals.append(internal);incidents.append(incident);types.append(kind)
            connected=0
            if len(marks)==2:
                common=incidents[0]&incidents[1];d=len(common)
                if types[0]!='L' and types[1]!='L':connected=4*internals[0]*internals[1]+2*d*(d-1)-4*cr(common)
                elif types[0]=='L' and types[1]=='L':connected=d*d-4*cr(common)
                else:connected=d*(d-1)-4*cr(common)
            value=self.tilted_coefficient(linears,internals,connected,left.xpower+right.xpower)
            total+=self.q**cr(pa)*ou_weight*value
            self.stats['accepted_term_diagrams']+=1
            self.stats['nonzero_charge_connected_diagrams']+=bool(len(marks)==2 and internals[0]*internals[1])
            self.stats['negative_diagram_weights']+=value<0
        self.stats['term_moment_cache_misses']+=1
        return total

    def moment(self,left,right,a=0,b=0,no_bridge=False):
        return sum((x*y*self.term_moment(s,t,a,b,no_bridge)
                    for s,x in left.items() for t,y in right.items()),F(0))

    def mehler_moment(self,a,b,r):
        aa,bb=dict(self.a(a)),dict(self.a(b))
        return sum((aa[d]*bb[d]*F(r)**d*self.P(d) for d in aa.keys()&bb.keys()),F(0))

    def serialize(self,field):
        return [{'length':t.length,'wick_blocks':t.wicks,'mark_blocks':t.marks,
                 'filter_blocks':t.filters,'charge_power':t.xpower,'coefficient':str(v)}
                for t,v in sorted(field.items(),key=lambda item:repr(item[0]))]
