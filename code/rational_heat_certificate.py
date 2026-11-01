"""Independent elementary enclosure using only integer/Fraction arithmetic.

Taylor/geometric tails, alternating Machin series and integer square roots
enclose the analytic heat-cutoff bound. No floating/transcendental library
is used for the certificate. This is not a matrix-quadrature certificate.
"""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json,hashlib

BITS=120
SCALE=1<<BITS
def down(x):return F((x*SCALE).numerator//(x*SCALE).denominator,SCALE)
def up(x):return -down(-x)

class I:
    def __init__(self,a,b=None):
        if isinstance(a,I):self.lo,self.hi=a.lo,a.hi;return
        self.lo=down(F(a));self.hi=up(F(a if b is None else b));assert self.lo<=self.hi
    def __add__(a,b):
        b=I(b);return I(a.lo+b.lo,a.hi+b.hi)
    __radd__=__add__
    def __neg__(a):return I(-a.hi,-a.lo)
    def __sub__(a,b):return a+-I(b)
    def __rsub__(a,b):return I(b)+-a
    def __mul__(a,b):
        b=I(b);v=[a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi];return I(min(v),max(v))
    __rmul__=__mul__
    def reciprocal(a):
        assert a.lo*a.hi>0;return I(1/a.hi,1/a.lo)
    def __truediv__(a,b):return a*I(b).reciprocal()
    def __rtruediv__(a,b):return I(b)*a.reciprocal()
    def __pow__(a,n):
        assert isinstance(n,int) and n>=0
        if n==2:
            v=[a.lo*a.lo,a.hi*a.hi];return I(0 if a.lo<=0<=a.hi else min(v),max(v))
        result=I(1)
        for _ in range(n):result=result*a
        return result
    def packed(a):return {'lower':[a.lo.numerator,a.lo.denominator],'upper':[a.hi.numerator,a.hi.denominator]}

def root(a):
    a=I(a);assert a.lo>=0
    def floor_sqrt(x):return F(isqrt((x.numerator*SCALE*SCALE)//x.denominator),SCALE)
    return I(floor_sqrt(a.lo),floor_sqrt(a.hi)+F(1,SCALE))

def exp_endpoint(x):
    if x<0:return exp_endpoint(-x).reciprocal()
    k=0
    while x>F(1,8):x/=2;k+=1
    term=F(1);total=term;N=24
    for j in range(1,N+1):term*=x/j;total+=term
    first=term*x/(N+1)
    remainder=first/(1-x/F(N+2))
    value=I(total,total+remainder)
    for _ in range(k):value=value**2
    return value

def exponential(a):
    a=I(a);return I(exp_endpoint(a.lo).lo,exp_endpoint(a.hi).hi)

def log_small(x):
    assert 1<=x<=2
    z=(x-1)/(x+1);N=40;total=F(0)
    for j in range(N):total+=2*z**(2*j+1)/F(2*j+1)
    remainder=2*z**(2*N+1)/(F(2*N+1)*(1-z*z))
    return I(total,total+remainder)

def log_endpoint(x):
    assert x>0
    if x<1:return -log_endpoint(1/x)
    k=0
    while x>2:x/=2;k+=1
    return log_small(x)+k*log_small(F(2))

def logarithm(a):
    a=I(a);return I(log_endpoint(a.lo).lo,log_endpoint(a.hi).hi)

def pi_interval():
    def atan_inverse(n):
        total=sum(((-1)**k*F(1,F((2*k+1)*n**(2*k+1)))) for k in range(40))
        return I(total,total+F(1,81*n**81))
    return 16*atan_inverse(5)-4*atan_inverse(239)

def certify(epsilon):
    pi=pi_interval();eps=I(epsilon);rho=I('0.0001')
    poles=[(1,-3),(1,-1),(1,1),(1,3),(-3,1),(-1,1),(3,1)]
    density=lambda d,R:exponential(-(d-R)**2/pi)/(pi*pi)
    tail=I(F(1,4))*exponential(-64/pi)
    holes=sum(4*pi*rho*rho*density(pi*root(i*i+j*j),2*rho) for i,j in poles)
    difference=root(tail+holes)
    z=exponential(9-3*pi);remainder=z/((1-z)*(1-exponential(-2*pi)))
    constant=9*root(2)/pi+4+4*remainder
    inv_integrals=[]
    for i,j in poles+[(1,1)]:
        d=pi*root(i*i+j*j);R=I(F(1,2));c=I(F(3,5))
        inv_integrals.append(2*pi*density(d,R)*logarithm(R/rho)
            +1/((1-c)*d)**2+exponential(-(c*d)**2/pi)/(R*R))
    gradient=root(2/pi+3*constant**2+6*pi*pi*sum(inv_integrals))+root(2*tail)
    gradient+=root(sum(4*pi*density(pi*root(i*i+j*j),2*rho) for i,j in poles))
    vacuum=root(2*(pi*eps)**2/((1+pi*eps)*(1+2*pi*eps)))
    total=2*vacuum+difference+pi*root(eps)*gradient
    values=dict(pi=pi,tail=tail,holes=holes,cutoff_difference=difference,gradient_total=gradient,total=total)
    threshold={'0.00000001':F(1,100),'0.0000000001':F(982,1000000),'0.000000000001':F(117,1000000)}[epsilon]
    assert total.hi<threshold,(epsilon,total.packed())
    return dict(epsilon=epsilon,values={k:v.packed() for k,v in values.items()},
        strict_rational_upper=[threshold.numerator,threshold.denominator],proved_upper=True)

def main():
    rootdir=Path(__file__).resolve().parents[1];out=rootdir/'results/stage4/heat-matrix/rational-certificate.json'
    if out.exists():raise RuntimeError('Do not overwrite')
    rows=[certify(e) for e in ['0.00000001','0.0000000001','0.000000000001']]
    data=dict(status='passed',rows=rows,dyadic_bits=BITS,exp_degree=24,log_terms=40,machin_terms=40,
        method='Exact Fraction operations rounded outward to dyadic endpoints; explicit series tails and integer isqrt.',
        scope='Certifies evaluation of the stated analytic heat-cutoff upper bound only; no matrix element quadrature.',
        hashes={'code/rational_heat_certificate.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    out.write_text(json.dumps(data,indent=2));print('passed',[(r['epsilon'],r['strict_rational_upper']) for r in rows])
if __name__=='__main__':main()
