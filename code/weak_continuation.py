"""Literal local weak integral and a Cauchy-pole subtraction.

This is a compact test pairing of one chiral factor, not a physical state
or a bounded pentagon. The Cauchy kernel identity is standard mathematics.
"""
from functools import lru_cache
import mpmath as m

HC=m.mpc(1,1)

def point(u):
    a=u.real;b=u.imag
    x=m.pi*(a*a+a+b*b)/b;y=-m.pi*(a+1)/b
    jac=m.pi**2*((a+1)**2+b*b)/b**3
    return x,y,y+1j*m.pi,jac

@lru_cache(None)
def constants(H,dps):
    q=m.exp(m.pi*1j*H);t=m.exp(-m.pi*1j/H)
    return q,t,-m.qp(t*t,t*t)/m.qp(q*q,q*q)

def normalized_direct(H,u):
    x,y,w,_=point(u);q,t,R=constants(H,m.mp.dps);z=x+H*y
    phi=m.qp(-t*m.exp(z/H),t*t)/m.qp(-q*m.exp(z),q*q)
    return w*phi/R

def regular(H,u):
    _,_,w,_=point(u);q,t,_=constants(H,m.mp.dps);d=w*(H-u)
    if not d:
        derivative=m.nsum(lambda k:q**(2*k)/(1-q**(2*k)),[1,m.inf])-m.nsum(lambda k:t**(2*k)/(1-t**(2*k)),[1,m.inf])/H
        return w*(derivative-m.mpf('.5'))
    e=m.qp(t*t*m.exp(d/H),t*t)/m.qp(t*t,t*t)*m.qp(q*q,q*q)/m.qp(q*q*m.exp(d),q*q)
    return w*(e/m.expm1(d)-1/d)

@lru_cache(None)
def nodes(n,dps):
    x,w=m.gauss_quadrature(n,'legendre')
    return [((x[i]+1)/2,w[i]/2) for i in range(n)]

def cauchy_exact(H,radius,power=8):
    z=H-HC;r2=radius**2
    if not z:return m.mpc(0)
    s=min(m.mpf(1),abs(z)**2/r2)
    return m.pi*r2/(power+1)*(1-(1-s)**(power+1))/z

def fixed_integral(H,nr,nt,radius,power=8,subtract=False):
    terms=[]
    for s,wr in nodes(nr,m.mp.dps):
        for k in range(nt):
            u=HC+radius*m.sqrt(s)*m.exp(2j*m.pi*k/nt)
            kernel=regular(H,u) if subtract else normalized_direct(H,u)
            terms.append(wr*(1-s)**power*kernel)
    return m.fsum(terms)*m.pi*radius**2/nt

def weak_subtracted(H,nr,nt,radius,power=8):
    return cauchy_exact(H,radius,power)+fixed_integral(H,nr,nt,radius,power,True)

def pole_centered_weak(H,nr,nt,radius,power=8):
    """Direct full product. Area Jacobian cancels its integrable 1/r pole."""
    offset=H-HC
    if abs(offset)>=radius:raise ValueError('Interior evaluation only')
    terms=[]
    for k in range(nt):
        e=m.exp(2j*m.pi*k/nt);projection=(offset*m.conj(e)).real
        edge=-projection+m.sqrt(projection**2+radius**2-abs(offset)**2)
        for s,wr in nodes(nr,m.mp.dps):
            r=edge*s;u=H+r*e
            g=max(m.mpf(0),1-abs(u-HC)**2/radius**2)**power
            terms.append(wr*r*edge*g*normalized_direct(H,u))
    return m.fsum(terms)*2*m.pi/nt
