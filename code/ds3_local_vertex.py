"""Declared G=0 local scalar insertion in global dS3, radius one.

BD modes are KG normalized. Finite spherical-mode matrix elements only;
this module does not construct a physical observer direct-integral algebra.
"""
from itertools import product
import mpmath as mp
from sympy.physics.wigner import gaunt

def convolution(p,q):
    out=[mp.mpc(0)]*(len(p)+len(q)-1)
    for j,a in enumerate(p):
        for k,b in enumerate(q):out[j+k]+=a*b
    return out

def mode_polynomial(nu,l,sign):
    """D_l ... D_1 exp(i sign nu t) = exp(...) P(tanh t)."""
    p=[mp.mpc(1)]
    for j in range(1,l+1):
        q=[mp.mpc(0)]*(len(p)+1)
        for k,c in enumerate(p):
            q[k]+=1j*sign*nu*c
            q[k+1]-=(k+j)*c
            if k:q[k-1]+=k*c
        p=q
    return p

def branches(nu,l,conjugate=False):
    normal=(1j**l/mp.sqrt(2*nu*(1-mp.exp(-2*mp.pi*nu))
              *mp.fprod(nu**2+j**2 for j in range(1,l+1))))
    terms=[]
    for sign,weight in [(-1,1),(1,-mp.exp(-mp.pi*nu))]:
        p=[normal*weight*c for c in mode_polynomial(nu,l,sign)]
        freq=sign*nu
        if conjugate:p=[mp.conj(c) for c in p];freq=-freq
        terms.append((freq,p))
    return terms

def mode(nu,l,t):
    y=mp.tanh(t)
    return mp.fsum(mp.exp(1j*w*t)*mp.polyval(p[::-1],y)
                   for w,p in branches(nu,l))/mp.cosh(t)

def legendre_mode(nu,l,t):
    """Marolf-Morrison (16)-(17), d=2, sigma=-1+i nu, type-2 P."""
    norm=mp.sqrt(mp.gamma(l+1-1j*nu)*mp.gamma(l+1+1j*nu)/2)
    return norm/mp.sqrt(mp.cosh(t))*mp.legenp(-mp.mpf('.5')+1j*nu,
                      -l-mp.mpf('.5'),1j*mp.sinh(t),type=2)

def moment(n,w):
    """Integral sech(t)*tanh(t)^n*exp(i w t) dt over the real axis."""
    a=(1+1j*w)/2;b=(1-1j*w)/2
    return mp.fsum(mp.binomial(n,k)*2**k*(-1)**(n-k)*mp.beta(a+k,b)
                   for k in range(n+1))

def angular(ls):
    return mp.mpf(str(gaunt(*ls,0,0,0).evalf(mp.mp.dps+5)))

def vertex(nus,ls=(0,0,0),conjugates=(True,False,True)):
    gg=angular(ls);terms=[]
    for choice in product(*(branches(n,l,c) for n,l,c in zip(nus,ls,conjugates))):
        freq=sum(x[0] for x in choice);p=[mp.mpc(1)]
        for _,q in choice:p=convolution(p,q)
        terms.extend(gg*c*moment(k,freq) for k,c in enumerate(p))
    value=mp.fsum(terms)
    return value,mp.fsum(abs(x) for x in terms)/max(abs(value),mp.mpf('1e-100'))

def direct_vertex(nus,ls,conjugates,T=90,step=3):
    def f(t):
        factors=[mode(n,l,t) for n,l in zip(nus,ls)]
        return mp.cosh(t)**2*mp.fprod(mp.conj(x) if c else x
                                     for x,c in zip(factors,conjugates))
    points=[-T+i*step for i in range(int(2*T/step)+1)]
    if points[-1]!=T:points.append(T)
    return angular(ls)*mp.quad(f,points)

def tail_bound(nus,ls,conjugates,T):
    norms=[mp.fsum(abs(a) for _,p in branches(n,l,c) for a in p)
           for n,l,c in zip(nus,ls,conjugates)]
    return 4*mp.exp(-T)*abs(angular(ls))*mp.fprod(norms)

def cosh_denominator(nus):
    S=sum(nus)
    return mp.cosh(mp.pi*S/2)*mp.fprod(mp.cosh(mp.pi*(S-2*n)/2) for n in nus)

def spherical_vertex(nus,conjugates):
    """Closed l=m=0 result for every conjugation pattern."""
    base=mp.pi/(4*mp.sqrt(4*mp.pi))*mp.sqrt(
        mp.fprod(mp.sinh(mp.pi*n)/n for n in nus))/cosh_denominator(nus)
    count=sum(conjugates)
    if count in [0,3]:return -base
    selected=conjugates.index(True) if count==1 else conjugates.index(False)
    return base*mp.cosh(mp.pi*nus[selected])

def boundary_exponents(nus,conjugates):
    labels=[(-1 if c else 1)*n for n,c in zip(nus,conjugates)]
    S=sum(labels)
    return [-mp.mpf('.5')+1j*(S-2*n)/2 for n in labels]

def sphere_period(a):
    """Normalized S2 triple chordal integral; J(0,0,0)=1."""
    A=sum(a)
    return 4**A*mp.gamma(A+2)*mp.fprod(mp.gamma(x+1) for x in a)/mp.fprod(
        mp.gamma(A-x+2) for x in a)

def future_phase(nu,l,conjugate=False):
    c=1j**l*mp.fprod((-1j*nu-j)/mp.sqrt(nu**2+j**2) for j in range(1,l+1))
    return mp.conj(c) if conjugate else c

def boundary_ll0(nus,l,conjugates):
    """Boundary trilinear on Y_l0,Y_l0,Y_00 with future-data phases."""
    # P_l(1-z/2) = hyper([-l,l+1],[1],z/4).
    a=boundary_exponents(nus,conjugates);parts=[]
    for k in range(l+1):
        shifted=a.copy();shifted[2]+=k
        coeff=mp.rf(-l,k)*mp.rf(l+1,k)/(mp.factorial(k)**2*4**k)
        parts.append(coeff*sphere_period(shifted))
    phases=future_phase(nus[0],l,conjugates[0])*future_phase(nus[1],l,conjugates[1])
    return (4*mp.pi)**mp.mpf('1.5')*phases*mp.fsum(parts)

def reduced_coefficient(nus,conjugates):
    return spherical_vertex(nus,conjugates)/boundary_ll0(nus,0,conjugates)

def reduced_norm(nus,conjugates):
    count=sum(conjugates);factor=mp.mpf(1)
    if count in [1,2]:
        selected=conjugates.index(True) if count==1 else conjugates.index(False)
        factor=mp.cosh(mp.pi*nus[selected])**2
    return factor/(64*mp.pi**3*cosh_denominator(nus))
