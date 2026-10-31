"""dS2 quantized-observer packet kernels from arXiv:2411.08090v1.

Equations 4.31, 4.35, 5.8, 5.6, 6.4 and 6.17; L=1, G->0.
The right two-point is a positive spectral quadratic form.  The LR
propagation kernel is real symmetric but is NOT assumed positive.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.special import loggamma, roots_legendre
from scipy.integrate import quad_vec

PI = np.pi


def mass(s):
    return np.sqrt(np.asarray(s, dtype=float)**2 + 0.25)


def s_of_p(p, m0):
    if m0 < 0.5:
        raise ValueError('principal-series lower mass must be >= 1/2')
    return np.sqrt((m0 + np.asarray(p, dtype=float)/(2*PI))**2 - 0.25)


def p_of_s(s, m0):
    return 2*PI*(mass(s)-m0)


def dp_ds(s):
    return 2*PI*np.asarray(s, dtype=float)/mass(s)


def rho_pl(s):
    s = np.asarray(s, dtype=float)
    return s*np.tanh(PI*s)/(4*PI**2)


def logcosh(x):
    x = np.abs(np.asarray(x, dtype=float))
    return x + np.log1p(np.exp(-2*x)) - np.log(2.0)


def logsinh_positive(x):
    x = np.asarray(x, dtype=float)
    with np.errstate(divide='ignore', invalid='raise'):
        return x + np.log(-np.expm1(-2*x)) - np.log(2.0)


def log_gamma8(nu, a, b):
    """All eight signs, evaluated as four complex-conjugate pairs."""
    if not np.isreal(nu) or float(nu) <= 0:
        raise ValueError('this implementation uses real principal-series nu>0')
    a, b = np.broadcast_arrays(np.asarray(a, float), np.asarray(b, float))
    out = np.zeros_like(a)
    for sa in (-1, 1):
        for sb in (-1, 1):
            out += 2*np.real(loggamma(0.25 + 0.5j*(nu+sa*a+sb*b)))
    return out


def log_hh_kernel(nu, a, b):
    a, b = np.broadcast_arrays(np.asarray(a, float), np.asarray(b, float))
    with np.errstate(divide='ignore'):
        return (np.log(a)+np.log(np.tanh(PI*a))+np.log(b)
                +np.log(np.tanh(PI*b))-np.log(64)-7*np.log(PI)
                +log_gamma8(nu,a,b))


def hh_kernel(nu, a, b):
    """Unnormalized HH spectral density D, NOT divided by Z_HH."""
    return np.exp(log_hh_kernel(nu,a,b))


def log_right_kernel(nu, a, b):
    """V(a|b); a spectral weight, not a normalized probability density."""
    a, b = np.broadcast_arrays(np.asarray(a, float), np.asarray(b, float))
    with np.errstate(divide='ignore'):
        return (np.log(a)+np.log(np.tanh(PI*a))+2*logcosh(PI*b)
                -np.log(16)-5*np.log(PI)+log_gamma8(nu,a,b))


def right_kernel(nu, a, b):
    return np.exp(log_right_kernel(nu,a,b))


def log_lr_kernel(nu, a, b):
    a, b = np.broadcast_arrays(np.asarray(a, float), np.asarray(b, float))
    with np.errstate(divide='ignore'):
        return (np.log(a)+logsinh_positive(PI*a)+np.log(b)
                +logsinh_positive(PI*b)-np.log(64)-7*np.log(PI)
                +log_gamma8(nu,a,b))


def lr_kernel(nu, a, b):
    """Raw ds' ds kernel from (4.35), not the kernel in flat dp measure."""
    return np.exp(log_lr_kernel(nu,a,b))


def lr_flat_kernel(nu, p_final, p_initial, m0):
    a, b = np.broadcast_arrays(s_of_p(p_final,m0),s_of_p(p_initial,m0))
    with np.errstate(divide='ignore'):
        log_jac = 0.5*(np.log(rho_pl(a))+np.log(dp_ds(a))
                       +np.log(rho_pl(b))+np.log(dp_ds(b)))
    return np.exp(log_lr_kernel(nu,a,b)-log_jac)


def gauss_interval(lo, hi, order):
    x, w = roots_legendre(order)
    return 0.5*(lo+hi)+0.5*(hi-lo)*x, 0.5*(hi-lo)*w


def composite_gauss(lo, hi, order, max_width=1.0):
    edges = np.linspace(lo, hi, int(np.ceil((hi-lo)/max_width))+1)
    pieces = [gauss_interval(a,b,order) for a,b in zip(edges[:-1],edges[1:])]
    return np.concatenate([x for x,w in pieces]),np.concatenate([w for x,w in pieces])


@dataclass(frozen=True)
class DoublePacket:
    c1: float = 3.0
    c2: float = 9.0
    width: float = 1.0

    def __post_init__(self):
        if not (self.width>0 and self.c1>self.width and self.c2-self.c1>2*self.width):
            raise ValueError('positive, disjoint lobe supports are required')

    @property
    def intervals(self):
        return [(self.c1-self.width,self.c1+self.width),
                (self.c2-self.width,self.c2+self.width)]

    def lobe(self,p,index,derivative=False):
        p = np.asarray(p,float)
        c = (self.c1,self.c2)[index]
        x = (p-c)/self.width
        if derivative:
            val = -np.sqrt(4/(3*self.width))*PI/(2*self.width)*np.sin(PI*x)
        else:
            val = np.sqrt(4/(3*self.width))*np.cos(PI*x/2)**2
        return np.where(np.abs(x)<=1,val,0.0)

    def amplitude(self,p,sign=1,derivative=False):
        if sign not in (-1,1):
            raise ValueError('sign must be +/-1')
        return (self.lobe(p,0,derivative)+sign*self.lobe(p,1,derivative))/np.sqrt(2)

    def probability(self,p,sign=1):
        return self.amplitude(p,sign)**2

    def nodes(self,order):
        pairs = [gauss_interval(lo,hi,order) for lo,hi in self.intervals]
        return np.concatenate([x for x,w in pairs]),np.concatenate([w for x,w in pairs])

    def physical_wavefunction(self,s,m0,sign=1):
        s=np.asarray(s,float)
        p=p_of_s(s,m0)
        with np.errstate(divide='ignore',invalid='ignore'):
            value=np.sqrt(dp_ds(s)/rho_pl(s))*self.amplitude(p,sign)
        return np.where(p>0,value,0.0)


def right_two_point(packet,m0,nu,sigma,times,p_order=32,a_order=24,tail_extra=12.0):
    """Mass-cutoff Q_m0 at EVERY insertion; no cutoff re-normalization."""
    if sigma <= 0:
        raise ValueError('production task requires proper-time Gaussian sigma>0')
    p,wp=packet.nodes(p_order)
    b=s_of_p(p,m0)
    s0=float(s_of_p(0,m0))
    amax=float(s_of_p(packet.intervals[-1][1],m0)+nu+tail_extra)
    a,wa=composite_gauss(s0,amax,a_order)
    delta=mass(a)[None,:]-mass(b)[:,None]
    weight=(right_kernel(nu,a[None,:],b[:,None])*np.exp(-(sigma*delta)**2)
            *(wp*packet.probability(p))[:,None]*wa[None,:])
    values=np.array([np.sum(weight*np.exp(1j*t*delta)) for t in times])
    return {'values':values,'a_cutoff':amax,'p_nodes':len(p),'a_nodes':len(a),
            'weight_mass':float(weight.sum()),'p_mass':float(np.dot(wp,packet.probability(p))),
            'minimum_weight':float(weight.min())}


def right_two_point_adaptive(packet,m0,nu,sigma,times,epsabs=2e-10):
    """Independent nested adaptive quadrature, with the a tail to infinity."""
    times=np.asarray(times,float)
    s0=float(s_of_p(0,m0))
    max_inner_error=[0.0]
    def outer(p):
        b=float(s_of_p(p,m0))
        mb=float(mass(b))
        def inner(a):
            delta=float(mass(a))-mb
            weight=float(np.exp(log_right_kernel(nu,a,b)-(sigma*delta)**2))
            return weight*np.exp(1j*times*delta)
        breaks=sorted(set([s0]+[x for x in (b-nu,b,b+nu,b+8) if x>s0]))
        total=np.zeros(len(times),complex)
        error=0.0
        for left,right in zip(breaks,breaks[1:]+[np.inf]):
            val,err=quad_vec(inner,left,right,epsabs=epsabs/20,epsrel=2e-10)
            total+=val
            error+=float(err)
        max_inner_error[0]=max(max_inner_error[0],error)
        return float(packet.probability(p))*total
    total=np.zeros(len(times),complex)
    outer_error=0.0
    for lo,hi in packet.intervals:
        val,err=quad_vec(outer,lo,hi,epsabs=epsabs,epsrel=2e-9)
        total+=val
        outer_error+=float(err)
    return {'values':total,'outer_estimated_error':outer_error,
            'maximum_inner_estimated_error':max_inner_error[0]}


def left_right_packet_matrix(packet,m0,nu,order=40):
    pairs=[gauss_interval(lo,hi,order) for lo,hi in packet.intervals]
    mat=np.zeros((2,2))
    for i,(p,wp) in enumerate(pairs):
        for j,(q,wq) in enumerate(pairs):
            kernel=lr_flat_kernel(nu,p[:,None],q[None,:],m0)
            mat[i,j]=np.einsum('i,ij,j->',wp*packet.lobe(p,i),kernel,wq*packet.lobe(q,j))
    signs={str(s):float(np.array([1,s])@mat@np.array([1,s])/2) for s in (-1,1)}
    return {'lobe_matrix':mat,'expectations':signs,
            'contrast':signs['1']-signs['-1']}