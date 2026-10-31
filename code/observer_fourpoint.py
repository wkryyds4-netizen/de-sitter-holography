"""Exact two-flavor uncrossed tasks of 2411.08090v1 (6.23)-(6.25).

dS2, L=1, G->0, same real principal-series nu for independent A/B.
Times are ordered LEFT TO RIGHT. These are complete AABB / ABBA
observables, NOT the sum or the ABAB OTOC of a single species.
"""
import numpy as np
from quantized_observer import (mass, s_of_p, right_kernel, composite_gauss)

class FourPointQuadrature:
    def __init__(self, packet, m0, nu, sigma, projected=True,
                 p_order=32, mass_order=20, tail_extra=12.):
        if sigma<=0: raise ValueError("Positive proper-time Gaussian width required")
        self.packet,self.m0,self.nu,self.sigma=packet,m0,nu,sigma
        p,wp=packet.nodes(p_order)
        self.b=s_of_p(p,m0)
        self.outer_weights=wp*packet.probability(p)
        self.lower=float(s_of_p(0,m0)) if projected else 0.
        self.upper=float(s_of_p(packet.intervals[-1][1],m0)+nu+tail_extra)
        self.x,wx=composite_gauss(self.lower,self.upper,mass_order)
        mb,mx=mass(self.b),mass(self.x)
        self.delta_external=mx[None,:]-mb[:,None]
        self.delta_grid=mx[None,:]-mx[:,None]
        self.external=(right_kernel(nu,self.x[None,:],self.b[:,None])
                       *np.exp(-(sigma*self.delta_external)**2)*wx[None,:])
        self.transition=(right_kernel(nu,self.x[None,:],self.x[:,None])
                         *np.exp(-(sigma*self.delta_grid)**2)*wx[None,:])
        self._external_cache={}
        self._grid_cache={}
        self.parameters=dict(m0=m0,nu=nu,sigma_proper=sigma,
                             projected=bool(projected),mass_lower=self.lower,
                             mass_upper=self.upper,p_order=p_order,mass_order=mass_order,
                             tail_extra=tail_extra,p_nodes=len(p),mass_nodes=len(self.x))

    def h_external(self,u):
        if u not in self._external_cache:
            self._external_cache[u]=np.sum(self.external*np.exp(1j*u*self.delta_external),axis=1)
        return self._external_cache[u]

    def h_grid(self,u):
        if u not in self._grid_cache:
            self._grid_cache[u]=np.sum(self.transition*np.exp(1j*u*self.delta_grid),axis=1)
        return self._grid_cache[u]

    def two_point(self,u):
        return self.outer_weights@self.h_external(u)

    def aabb(self,times):
        t1,t2,t3,t4=times
        return self.outer_weights@(self.h_external(t1-t2)*self.h_external(t3-t4))

    def abba(self,times):
        t1,t2,t3,t4=times
        outer=self.external*np.exp(1j*(t1-t4)*self.delta_external)
        inner=self.h_grid(t2-t3)
        return self.outer_weights@(outer@inner)

    def values(self,times):
        return dict(AABB=self.aabb(times),ABBA=self.abba(times))

    def zero_diagnostics(self):
        h=self.h_external(0).real
        c2=float(self.outer_weights@h)
        d=float(self.aabb([0,0,0,0]).real)
        n=float(self.abba([0,0,0,0]).real)
        variance=float(self.outer_weights@((h-c2)**2))
        return dict(C2=c2,AABB=d,ABBA=n,C2_squared=c2*c2,
                    adjacent_excess=d-c2*c2,variance_of_conditional_C2=variance,
                    nested_minus_adjacent=n-d,
                    normalized_adjacent=d/(c2*c2),normalized_nested=n/(c2*c2))
