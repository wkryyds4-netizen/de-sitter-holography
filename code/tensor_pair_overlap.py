"""Finite Hermite pair cubature via weighted Cauchy-Schwarz and Mehler.

Uses DLMF 18.18.28 and Johansson 1802.07942 section 2. These are
standard bounds, not new identities. All saved coefficients are dyadic.
"""
import numpy as np
from flint import arb,acb,acb_mat
from tensor_overlap import GaussianOverlap
from chirp_residual_certificate import real_exact,complex_exact,norm_squared

def fourier(c):
    # Multiplication by these four phases only swaps/signs dyadic components.
    return c*np.array([1,-1j,-1,1j],dtype=complex)[np.indices(c.shape).sum(axis=0)%4]

def energies(c):
    n=c.shape[0];diag=[arb(0) for _ in range(2*n-1)]
    x=[arb(0) for _ in range(n)];y=[arb(0) for _ in range(n)]
    for i in range(n):
        for j in range(n):
            e=real_exact(c[i,j].real)**2+real_exact(c[i,j].imag)**2
            diag[i+j]+=e;x[i]+=e;y[j]+=e
    return diag,x,y

def weighted(coef,r):
    result=arb(0)
    for c in reversed(coef):result=result/r+c
    return result

class PairOverlap(GaussianOverlap):
    def __init__(self,v,w,terms=8,order=32,unit=False):
        assert v.shape==w.shape
        super().__init__(w,terms,order,unit)
        self.V=acb_mat(self.n,self.n,[complex_exact(z) for z in v.flat])
        self.nv=norm_squared(v);self.nw=norm_squared(w)
        self.ev=energies(v);self.ew=energies(w)
        self.weighted=[]
        for k in range(1,9):
            r=1-arb(1)/(2**k)
            pref=(weighted(self.ev[0],r)*weighted(self.ew[0],r)).sqrt()/(self.pi*(1-r*r))
            self.weighted.append((r,pref))

    def amplitude(self,x,y):
        candidates=[]
        for r,pref in self.weighted:
            a=(1-r)/(1+r)
            bound=pref*(-a*(x.real*x.real+y.real*y.real)
                        +(x.imag*x.imag+y.imag*y.imag)/a).exp()
            candidates.append(bound.upper())
        return min(candidates,key=float)

    def mass_tail(self,energy,L):
        candidates=[]
        for r,_ in self.weighted:
            a=(1-r)/(1+r)
            e=weighted(energy[1],r)+weighted(energy[2],r)
            candidates.append((e*(a.sqrt()*L).erfc()/(1-r)).upper())
        return min(candidates,key=float)

    def tails(self,L):
        tv=self.mass_tail(self.ev,arb(L));tw=self.mass_tail(self.ew,arb(L))
        return tv,tw,(tv*tw).sqrt()

    def phase_bound(self,z,t):
        """z complex, t real. Bound ratios using their deviation from unitarity.

        |1+exp(u+iv)|/exp(u) = sqrt((exp(-u)-1)^2+4exp(-u)cos(v/2)^2).
        A zero-containing lower bound rejects holomorphy for the whole box.
        """
        assert t.imag.is_zero()
        x=self.scale*z.real;b=self.scale*z.imag;y=self.scale*t.real
        result=arb(1)
        def positive_lower(value):
            lo=value.lower()
            return lo if lo>0 else arb(0)
        def denominator(u,v):
            q=(-u).exp()
            # arb ** 2 may return NaN on a zero-containing interval.
            # Only nonnegative exact lower endpoints are squared here.
            a=positive_lower(abs(q-1));c=positive_lower(abs((v/2).cos()))
            return (a*a+4*positive_lower(q)*c*c).sqrt()
        for j in range(self.terms):
            a=(2*j+1)*self.pi
            d1=denominator(x-a,b+y);d2=denominator(y-a-b,x)
            if not (d1>0 and d2>0):return arb('inf')
            result*=1+2*abs(b.sin())/d1
            result*=1+abs((2*b).expm1())/d2
        return result.upper()

    def bound(self,rect,depth=0):
        a,b,c,d=map(arb,rect);hx=(b-a)/2;hy=(d-c)/2;mx=(a+b)/2;my=(c+d)/2
        xr=acb(arb(mx,hx));yr=acb(arb(my,hy))
        direct=(4*hx*hy*self.amplitude(xr,yr)).upper()
        xe=acb(arb(mx,5*hx/4),arb(0,3*hx/4))
        ye=acb(arb(my,5*hy/4),arb(0,3*hy/4))
        px=arb(1) if self.unit else self.phase_bound(xe,yr)
        py=arb(1) if self.unit else self.phase_bound(ye,xr)
        kind='direct';error=direct
        if px.is_finite() and py.is_finite():
            quad=(100*hx*hy/arb(2)**(2*self.order)*
                  (self.amplitude(xe,yr)*px+self.amplitude(xr,ye)*py)).upper()
            if quad<direct:kind='gauss';error=quad
        return dict(rect=rect,depth=depth,kind=kind,error=error)

    def quadrature(self,rect):
        a,b,c,d=rect
        hx,_,wx,px,mx,rx=self.basis(a,b)
        hy,_,wy,py,my,ry=self.basis(c,d)
        w=hx.transpose()*self.C*hy
        v=hx.transpose()*self.V*hy
        total=acb(0)
        for i in range(self.order):
            for j in range(self.order):
                phase=acb(1)
                if not self.unit:
                    for k in range(self.terms):
                        phase*=(1+rx[k][i]*my[j])/(1+rx[k][i]*py[j])
                        phase*=(1+ry[k][j]*mx[i])/(1+ry[k][j]*px[i])
                total+=wx[i]*wy[j]*w[i,j]*v[i,j]*phase
        assert total.is_finite()
        return total
