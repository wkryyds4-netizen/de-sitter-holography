"""Actual hard-threshold limiting curve; empirical quadrature, no finite-M data.

The exponential packet is integrated out exactly. All transverse max-folds
are split explicitly. No asymptotic coefficient is used to generate a curve.
"""
from pathlib import Path
from functools import lru_cache
import numpy as np
from scipy.special import roots_legendre, loggamma, erf
from preparation_asymptotics import f, grid, cj, NU, SIGMA, BETA

ROOT = Path(__file__).resolve().parents[1]

@lru_cache(None)
def roots(n):
    return roots_legendre(n)

def knots(radius):
    return np.array(sorted(set([-radius, radius] +
        [v for v in [-8,-6,-4,-3,-2,-1,0,1,2,3,4,6,8] if abs(v)<radius])))

def rule(lo, hi, n, breaks):
    if hi <= lo:
        return np.empty(0), np.empty(0)
    edges = sorted(set([float(lo), float(hi)] + [float(v) for v in breaks if lo<v<hi]))
    return grid(edges,n)

def partial(b, chi, n, radius, lower, tilt=0.):
    """For an array of b, integrate f(v)*exp(tilt*v) below/above b.

    Fixed smooth-kernel panels are clipped at each moving endpoint. The
    endpoint itself is never sampled as an indicator inside a panel.
    """
    b,restore=np.unique(np.atleast_1d(b),return_inverse=True)
    edges=knots(radius)
    lo=np.broadcast_to(edges[:-1],(len(b),len(edges)-1)).copy()
    hi=np.broadcast_to(edges[1:],lo.shape).copy()
    if lower:
        hi=np.minimum(hi,b[:,None]); hi=np.maximum(hi,lo)
    else:
        lo=np.maximum(lo,b[:,None]); lo=np.minimum(lo,hi)
    x,w=roots(n)
    v=(lo+hi)[:,:,None]/2+(hi-lo)[:,:,None]*x/2
    value=f(v,chi)*np.exp(tilt*v)
    return np.sum(value*(hi-lo)[:,:,None]*w/2,axis=(1,2))[restore]

def jq_slice(chi, n=24, radius=10.):
    """Exact slice geometry of exp[-beta max(0,-u-c/2,-v-c/2,-u-v)]."""
    u,wu=rule(-radius,radius,n,list(knots(radius))+[-chi/2,chi/2])
    a=np.maximum(0.,-u-chi/2)
    b=np.maximum(-chi/2,-u)
    turn=b-a
    inner=(np.exp(-BETA*b)*partial(turn,chi,n,radius,True,BETA)
           +np.exp(-BETA*a)*partial(turn,chi,n,radius,False))
    return np.dot(wu*f(u,chi),inner)

def jq_ordered(chi,n=24,radius=10.):
    """Independent geometry: use symmetry to integrate only v<=u.

    In this half-plane the last two competitors switch at u=chi/2;
    the zero-weight boundary is v=-chi/2 or v=-u.
    """
    u,wu=rule(-radius,radius,n,list(knots(radius))+[0.,chi/2,-chi/2])
    d=np.where(u>=chi/2,chi/2,u)
    turn=np.minimum(u,-d)
    lower=partial(turn,chi,n,radius,True,BETA)
    # Integrate the upper segment directly, avoiding two large primitives.
    edges=knots(radius)
    lo=np.maximum(edges[:-1][None,:],(-d)[:,None])
    hi=np.minimum(edges[1:][None,:],u[:,None]);hi=np.maximum(hi,lo)
    x,w=roots(n)
    v=(lo+hi)[:,:,None]/2+(hi-lo)[:,:,None]*x/2
    upper=np.sum(f(v,chi)*(hi-lo)[:,:,None]*w/2,axis=(1,2))
    return 2*np.dot(wu*f(u,chi),np.exp(BETA*d)*lower+upper)

def plateau(n=32,radius=10.):
    """Independent positive path integrals for N, X_infty and D_edge."""
    x,wx=rule(0.,radius,n,knots(radius))
    kpos=f(x,0).real/(16*np.pi**3)
    kneg=f(-x,0).real/(16*np.pi**3)
    gain=np.dot(wx,kpos);loss=np.dot(wx,kneg*np.exp(-BETA*x))
    mix=[];difference=[]
    for xx in x:
        u,wu=rule(0.,radius,n,list(knots(radius))+[xx])
        weight=np.exp(-BETA*np.maximum(0.,u-xx))
        diff=weight*(-np.expm1(-BETA*np.minimum(xx,u)))
        ku=f(-u,0).real/(16*np.pi**3)
        mix.append(np.dot(wu,ku*weight));difference.append(np.dot(wu,ku*diff))
    mixed=np.dot(wx*kpos,mix)
    norm=gain*gain+gain*loss+mixed+loss*loss
    cross=(gain+loss)**2
    edge=2*np.dot(wx*kpos,difference)
    return dict(N=float(norm),Xinf=float(cross),Dedge=float(edge),
                difference_identity_error=float(abs(2*(norm-cross)-edge)))

def prefactor(chi):
    chi=np.asarray(chi,dtype=complex)
    return np.exp(1j*chi*np.log(2)+loggamma(-1j*chi/2)
        -loggamma((1+1j*chi)/2)-SIGMA**2*chi**2/2)/(1024*np.pi**6.5)

def fourier(chi, weights, jq, xinf, times, subtraction_width=1.):
    """Return complex X plus directly formed early difference Xinf-X."""
    a=subtraction_width
    amp=prefactor(chi)*jq
    reg=amp-1j*xinf*np.exp(-a*chi*chi)/(2*np.pi*chi)
    phase=np.exp(1j*np.asarray(times)[:,None]*chi[None,:])
    integ=phase@(weights*reg)
    erfv=erf(np.asarray(times)/(2*np.sqrt(a)))
    cross=xinf/2*(1-erfv)+integ
    diff=xinf/2*(1+erfv)-integ
    absolute_scale=np.sum(np.abs(weights*reg))
    return cross,diff,float(absolute_scale)
