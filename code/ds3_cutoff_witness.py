"""A finite-K group-averaged witness for the actual zero-time order difference."""
import itertools
import numpy as np
from scipy.special import roots_legendre
from ds3_crossed_packet import probability


def odd_kernel(a,b,p):
    # Positive difference evaluated without subtracting close reciprocals.
    A=np.cosh(np.pi*(a+b));B=np.cosh(np.pi*(a-b));C=np.cosh(np.pi*p)
    return 2*np.sinh(np.pi*a)*np.sinh(np.pi*b)/((A+C)*(B+C))


def coefficient(p):return np.sqrt(np.sinh(np.pi*p)/p)/(4*np.sqrt(np.pi))


def four_spherical(a,q,c,p):
    a,c=np.broadcast_arrays(a,c);total=np.zeros(a.shape)
    for signs in itertools.product([-1,1],repeat=4):
        x=np.pi*(signs[0]*a+signs[1]*q+signs[2]*c+signs[3]*p)/2
        with np.errstate(invalid='ignore',divide='ignore'):
            val=np.where(abs(x)<1e-4,-x*x/3+x**4/45-2*x**6/945,1-x/np.tanh(x))
        total+=np.prod(signs)*val
    return total/(16*a*q*c*p)


def witness(order=24,s=.2,p=.6,q=.9,a0=1.1,a1=1.4,c0=.4,c1=.8):
    x,w=roots_legendre(order)
    a=a0+(a1-a0)*(x+1)/2;aw=w*(a1-a0)/2
    c=c0+(c1-c0)*(x+1)/2;cw=w*(c1-c0)/2
    b=s*(x+1)/2;bw=w*s/2
    f2=probability(c,c0,c1)/(2*c*c/np.pi)
    aa=a[:,None,None];cc=c[None,:,None];bb=b[None,None,:]
    low=np.sum((odd_kernel(aa,bb,q)*odd_kernel(bb,cc,p)-
                odd_kernel(aa,bb,p)*odd_kernel(bb,cc,q))*bw,axis=-1)
    pref=coefficient(p)*coefficient(q)*np.sqrt(np.cosh(np.pi*c)[None,:]/np.cosh(np.pi*a)[:,None])
    pref*=np.exp(np.pi*(np.sqrt(1+a*a)[:,None]-np.sqrt(1+c*c)[None,:])/2)
    density=aw[:,None]*cw[None,:]*f2[None,:]
    overlap=-np.sum(density*pref*low)
    testnorm=np.sum(density*four_spherical(a[:,None],q,c[None,:],p))
    return dict(order=order,s=s,p=p,q=q,a0=a0,a1=a1,c0=c0,c1=c1,
                overlap=float(overlap),test_norm_squared=float(testnorm),
                order_difference_lower_bound=float(overlap*overlap/testnorm),
                sign_condition=float(np.cosh(np.pi*p)*np.cosh(np.pi*q)-np.sinh(np.pi*s)**2))
