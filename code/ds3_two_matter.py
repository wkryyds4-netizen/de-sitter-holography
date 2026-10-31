"""Two distinguishable matter species in the declared G=0 local observer model."""
import itertools
import mpmath as mp
import numpy as np
from ds3_group_pairing import group_density,triple_group,joint,state_density,emission000
from ds3_local_vertex import boundary_ll0

def pair_weight(a,q,b):
    return group_density(b)*triple_group([a,q,b])

def pair_weight_ds(a,q,b):
    period=boundary_ll0([a,q,b],0,(False,False,True))
    return b*b/(4*mp.pi**4)*abs(period)**2

def four_weight(a,b,c,p,q):
    return group_density(b)*triple_group([a,b,q])*triple_group([b,c,p])

def four_spherical_closed(nus):
    def primitive(k):
        x=mp.pi*k/2
        return 0 if not x else 1-x*mp.coth(x)
    return mp.fsum(mp.fprod(signs)*primitive(sum(s*n for s,n in zip(signs,nus)))
                   for signs in itertools.product([-1,1],repeat=4))/(16*mp.fprod(nus))

def chain(a,b,c,p,q):
    return joint(a,b,q)*joint(b,c,p)/state_density(b)

def balanced(a,b,p):
    A=mp.cosh(mp.pi*(a+b));B=mp.cosh(mp.pi*(a-b));C=mp.cosh(mp.pi*p)
    return (A-B)/(8*mp.pi*(A+C)*(B+C))

def clock_g(x):
    delta=1/(mp.sqrt(1+x*x)+x)
    return x*mp.exp(mp.pi*delta)/(1+mp.exp(-2*mp.pi*x))

def tilted(a,b,p):
    return clock_g(a)/clock_g(b)*balanced(a,b,p)

def k(w,p):
    return 1/(8*mp.pi*(mp.cosh(mp.pi*w)+mp.cosh(mp.pi*p)))

def lr(t,p):
    return (p if not t else mp.sin(p*t)/mp.sinh(t))/(4*mp.pi*mp.sinh(mp.pi*p))

def tilt_constant(s):
    return 2*mp.exp(mp.pi/(mp.sqrt(1+s*s)+s))

def tail_bound(R,s):
    B=tilt_constant(s);pi=mp.pi
    return mp.exp(-pi*R)/(4*pi*pi)*(
        (1+B)*(R/pi+1/pi**2)+B/s*(R*R/pi+2*R/pi**2+2/pi**3))

def limit_bound(M,s,p,q,T):
    assert M>=2*s
    pi=mp.pi
    L=2/M+2*pi/M**2+2*pi*mp.exp(-pi*M)
    C=mp.exp(M*L/2)
    eps=2*(2+mp.cosh(pi*p)+mp.cosh(pi*q))*mp.exp(-pi*M)
    return C*(L+2*T/M**2)/(2*pi**5)+eps*lr(0,p)*lr(0,q)+tail_bound(M/2,s)

def quadrature(M,s,p,q,tp,tq,order=64,D=14):
    """Conditional initial mass M, both creations shifted by i*pi/2.

    u=b-M and v=a-b are clipped to [-D,D], with exact observer cutoffs.
    Splits include both zero differences and the moving lower-bound corner.
    """
    xs,ws=np.polynomial.legendre.leggauss(order)
    def mesh(lo,hi,points):
        breaks=sorted(set([lo,hi]+[x for x in points if lo<x<hi]))
        return (np.concatenate([l+(r-l)*(xs+1)/2 for l,r in zip(breaks,breaks[1:])]),
                np.concatenate([ws*(r-l)/2 for l,r in zip(breaks,breaks[1:])]))
    def g(x):
        return x*np.exp(np.pi/(np.sqrt(1+x*x)+x))/(1+np.exp(-2*np.pi*x))
    def bal(a,b,p):
        A=np.cosh(np.pi*(a+b));B=np.cosh(np.pi*(a-b));C=np.cosh(np.pi*p)
        return (A-B)/((A+C)*(B+C)*8*np.pi)
    us,uw=mesh(max(s-M,-D),D,[0,s+D-M])
    result=0j
    for u,w in zip(us,uw):
        b=M+u
        vs,vw=mesh(max(s-b,-D),D,[0])
        a=b+vs
        ep=u*(b+M)/(np.sqrt(1+b*b)+np.sqrt(1+M*M))
        eq=vs*(a+b)/(np.sqrt(1+a*a)+np.sqrt(1+b*b))
        kernel=g(a)/g(M)*bal(a,b,q)*bal(b,M,p)
        result+=w*np.sum(vw*kernel*np.exp(-1j*(tp*ep+tq*eq)))
    return complex(result)
