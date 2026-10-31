"""dS2 large observer-mass kernel and common-clock noncrossing limit."""
import numpy as np
from scipy.special import loggamma
from quantized_observer import composite_gauss

def stationary_transition(delta,nu):
    delta=np.asarray(delta,float)
    lg=2*np.real(loggamma(.25+.5j*(nu+delta))+loggamma(.25+.5j*(nu-delta)))
    return np.exp(lg-np.pi*delta-np.log(16)-3*np.log(np.pi))

def packet_survival(packet,q):
    q=np.asarray(q,float)
    cdf=0.
    for center in (packet.c1,packet.c2):
        x=np.clip((q-center)/packet.width,-1,1)
        cdf+=.5*((x+1)/2+2/(3*np.pi)*np.sin(np.pi*x)+np.sin(2*np.pi*x)/(12*np.pi))
    return 1-cdf

class ClockLimit:
    def __init__(self,packet,nu,sigma,order=28,p_order=48,extra=14.):
        p,wp=packet.nodes(p_order)
        self.outer=wp*packet.probability(p)
        self.r=p/(2*np.pi)
        self.z,self.w=composite_gauss(0,packet.intervals[-1][1]/(2*np.pi)+nu+extra,order)
        self.de=self.z[None,:]-self.r[:,None]
        self.dg=self.z[None,:]-self.z[:,None]
        self.ext=stationary_transition(self.de,nu)*np.exp(-(sigma*self.de)**2)*self.w[None,:]
        self.grid=stationary_transition(self.dg,nu)*np.exp(-(sigma*self.dg)**2)*self.w[None,:]
    def h_ext(self,t):return np.sum(self.ext*np.exp(1j*t*self.de),axis=1)
    def h_grid(self,t):return np.sum(self.grid*np.exp(1j*t*self.dg),axis=1)
    def values(self,times):
        a,b,c,d=times
        return dict(C2=self.outer@self.h_ext(0),
                    AABB=self.outer@(self.h_ext(a-b)*self.h_ext(c-d)),
                    ABBA=self.outer@((self.ext*np.exp(1j*(a-d)*self.de))@self.h_grid(b-c)))

def energy_tail_moments(packet,nu,sigma,order=28,radius=14.):
    """Independent proper-frequency quadrature using survival K(2*pi*q).
    Split the AABB max diagonal exactly, not a square indicator grid.
    """
    def nodes(lo,hi,breaks=()):
        edges=sorted(set([lo,hi]+[float(x) for x in breaks if lo<x<hi]))
        pieces=[composite_gauss(a,b,order) for a,b in zip(edges[:-1],edges[1:])]
        return np.concatenate([x for x,w in pieces]),np.concatenate([w for x,w in pieces])
    endpoints=np.array(packet.intervals).ravel()/(2*np.pi)
    energy,w=nodes(-radius,radius,[0,*endpoints])
    r=stationary_transition(-energy,nu)*np.exp(-(sigma*energy)**2)
    K=packet_survival(packet,2*np.pi*np.maximum(energy,0))
    c2=np.dot(w,r*K)
    cdf=[];positive=[]
    for e in energy:
        x,wx=nodes(-radius,float(e))
        cdf.append(np.dot(wx,stationary_transition(-x,nu)*np.exp(-(sigma*x)**2)))
        pos,wp=nodes(0,radius,[-e,*(endpoints-e)])
        rp=stationary_transition(-pos,nu)*np.exp(-(sigma*pos)**2)
        positive.append(np.dot(wp,rp*packet_survival(packet,2*np.pi*np.maximum(e+pos,0))))
    adjacent=2*np.dot(w,r*K*np.asarray(cdf))
    neg,wn=nodes(-radius,0)
    f0=np.dot(wn,stationary_transition(-neg,nu)*np.exp(-(sigma*neg)**2))
    nested=np.dot(w,r*(K*f0+np.asarray(positive)))
    return dict(C2=float(c2),AABB=float(adjacent),ABBA=float(nested),radius=radius,order=order)
