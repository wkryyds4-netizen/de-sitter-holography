"""Spacetime L2 spectral normalization in dS3, radius one; not group averaging."""
import mpmath as mp
from ds3_local_vertex import mode_polynomial, mode, reduced_coefficient

def d_l(nu,l):
    return mp.sqrt(mp.fprod(nu**2+j**2 for j in range(1,l+1)))

def continuum(nu,l,t,sign=-1):
    p=mode_polynomial(nu,l,sign)
    return 1j**l*mp.exp(1j*sign*nu*t)*mp.polyval(p[::-1],mp.tanh(t))/(mp.sqrt(2*mp.pi)*d_l(nu,l))

def mu(nu):
    return nu*mp.tanh(mp.pi*nu)/mp.pi

def bd_vector(nu):
    q=mp.exp(-mp.pi*nu)
    return mp.matrix([1,-q])/mp.sqrt(1+q*q)

def orth_vector(nu):
    q=mp.exp(-mp.pi*nu)
    return mp.matrix([q,1])/mp.sqrt(1+q*q)

def gram(nu):
    a=1/mu(nu);b=-mp.pi/(nu*mp.sinh(mp.pi*nu))
    return mp.matrix([[a,b],[b,a]])

def projector(nu):
    v=bd_vector(nu)
    return v*v.T

def bound_polynomial(l,k):
    # sech(t)^k * polynomial(tanh(t)); normalized Darboux seed and chain.
    assert 1<=k<=l
    p=[1/mp.sqrt(mp.beta(mp.mpf('.5'),k))]
    for j in range(k+1,l+1):
        out=[mp.mpf(0)]*(len(p)+1)
        for n,c in enumerate(p):
            out[n+1]-=(n+k+j)*c
            if n:out[n-1]+=n*c
        p=[c/mp.sqrt(j*j-k*k) for c in out]
    return p

def bound(l,k,t):
    return mp.sech(t)**k*mp.polyval(bound_polynomial(l,k)[::-1],mp.tanh(t))

def sech_moment(p,n,nu):
    return 2**(p-1)*mp.fsum(mp.binomial(n,j)*2**j*(-1)**(n-j)*
        mp.beta((p+1j*nu)/2+j,(p-1j*nu)/2) for j in range(n+1))

def sech_transform(p,l,nu,sign=-1):
    # Direct finite-polynomial integral, independent of the operator completeness proof.
    poly=mode_polynomial(nu,l,sign)
    return (-1j)**l*mp.fsum(mp.conj(c)*sech_moment(p,n,-sign*nu)
        for n,c in enumerate(poly))/(mp.sqrt(2*mp.pi)*d_l(nu,l))

def sech_bound_overlap(p,l,k):
    return mp.fsum(c*mp.beta(mp.mpf(n+1)/2,(p+k)/2)
        for n,c in enumerate(bound_polynomial(l,k)) if n%2==0)

def observer_coefficient(nus,conjugates):
    # Slots (out,in,probe); only observer legs change from KG to spacetime normalization.
    return mp.sqrt(mu(nus[0])*mu(nus[1]))*reduced_coefficient(nus,conjugates)
