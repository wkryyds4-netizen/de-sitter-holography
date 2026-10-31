"""Physical weights in the previously normalized scalar crossing channel.

G=0, ell=1, BD observers, KG matter, fixed positive observer threshold,
half-period regulated creation. M is a representation label, not energy.
No extra Plancherel density belongs in these Lebesgue-label integrals.
"""
import mpmath as mp
from ds3_two_matter import balanced, clock_g, k, four_weight
from ds3_group_pairing import emission000, group_density


def energy(x):
    return mp.sqrt(1+x*x)


def energy_difference(x,y):
    return (x-y)*(x+y)/(energy(x)+energy(y))


def labels(M,u,v,r):
    return M+u+v, M+u, M, M+v+r


def amplitude(M,u,v,p,q,r):
    a,b,c,d=labels(M,u,v,r)
    return clock_g(a)/clock_g(c)*mp.sqrt(
        balanced(a,b,q)*balanced(b,c,p)*balanced(a,d,p)*balanced(d,c,q))


def raw_amplitude(M,u,v,p,q,r):
    """Independent E/W/rho assembly, with the initial density divided out."""
    a,b,c,d=labels(M,u,v,r)
    es=emission000(a,b,q)*emission000(b,c,p)
    et=emission000(a,d,p)*emission000(d,c,q)
    return mp.conj(es)*et*mp.exp(mp.pi*energy_difference(a,c))/(
        group_density(c)*mp.sqrt(four_weight(a,b,c,p,q)*four_weight(a,d,c,q,p)))


def phase_energy(M,u,v,r,kind='crossed'):
    a,b,c,d=labels(M,u,v,r)
    if kind=='crossed':
        return energy_difference(a,b)-energy_difference(d,c)
    if kind=='exchange':
        return energy_difference(d,b)
    raise ValueError(kind)


def weight(M,u,v,p,q,r,t=0,kind='crossed'):
    return amplitude(M,u,v,p,q,r)*mp.exp(1j*t*phase_energy(M,u,v,r,kind))


def limiting_weight(u,v,p,q,r,t=0,kind='crossed'):
    phase=-r if kind=='crossed' else v+r-u
    return mp.sqrt(k(v,q)*k(u,p)*k(u-r,p)*k(v+r,q))*mp.exp(1j*t*phase)


def log_k_derivatives(x,p):
    C=mp.cosh(mp.pi*p);B=mp.cosh(mp.pi*x)
    return -mp.pi*mp.sinh(mp.pi*x)/(B+C), -mp.pi**2*(1+C*B)/(B+C)**2


def limit_r_derivatives(u,v,p,q,r,t=0,kind='crossed'):
    val=limiting_weight(u,v,p,q,r,t,kind)
    lp,lpp=log_k_derivatives(u-r,p);lq,lqq=log_k_derivatives(v+r,q)
    first=(-lp+lq)/2+1j*t*(-1 if kind=='crossed' else 1)
    second=(lpp+lqq)/2
    return val,val*first,val*(first**2+second)
