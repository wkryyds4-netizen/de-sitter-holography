"""Validated tensor Gauss-Legendre cubature for Gaussian/Hermite overlaps.

Uses Johansson (1802.07942, section 2): C_rho < 50 for rho > 1.1.
Here rho=2. Real cells can always be enclosed using |F_N|=1.
"""
from pathlib import Path
import sys,heapq
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'tmp/arb-runtime'))
from flint import arb,acb,acb_mat
from chirp_residual_certificate import complex_exact,norm_squared
from residual_overlap_rectangle import phase_product

class GaussianOverlap:
    def __init__(self,w,terms=3,order=24,unit=False):
        self.n=w.shape[0];self.norm=norm_squared(w).sqrt()
        self.C=acb_mat(self.n,self.n,[complex_exact(z).conjugate() for z in w.flat])
        self.terms=terms;self.order=order;self.unit=unit
        self.pi=arb.pi();self.scale=self.pi.sqrt()
        self.nodes=[arb.legendre_p_root(order,k,weight=True) for k in range(order)]
        self.cache={}

    def amplitude(self,x,y):
        # Bound |conj(W) * vacuum| on complex x,y, using conjugate coefficients.
        candidates=[]
        for k in range(1,9):
            r=1-arb(1)/(2**k)
            bound=self.norm*r**(-(self.n-1))/(self.pi*(1-r*r).sqrt())
            bound*=(-(x.real*x.real+y.real*y.real)/(1+r)
                    +(x.imag*x.imag+y.imag*y.imag)/(1-r)).exp()
            candidates.append(bound.upper())
        return min(candidates,key=float)

    def bound(self,rect,depth=0):
        a,b,c,d=map(arb,rect);hx=(b-a)/2;hy=(d-c)/2;mx=(a+b)/2;my=(c+d)/2
        xr=acb(arb(mx,hx));yr=acb(arb(my,hy));area=4*hx*hy
        direct=(area*self.amplitude(xr,yr)).upper()
        xe=acb(arb(mx,5*hx/4),arb(0,3*hx/4))
        ye=acb(arb(my,5*hy/4),arb(0,3*hy/4))
        bounds=[]
        for x,y in [(xe,yr),(xr,ye)]:
            phase=acb(1) if self.unit else phase_product(self.scale*x,self.scale*y,self.terms)
            if not phase.is_finite():break
            bounds.append(self.amplitude(x,y)*abs(phase))
        kind='direct';error=direct
        if len(bounds)==2:
            # (Ix-Qx)Iy + Qx(Iy-Qy), with positive GL weights summing to 2.
            quad=(2*hx*hy*arb(50)/arb(2)**(2*self.order)*sum(bounds,arb(0))).upper()
            if quad < direct:kind='gauss';error=quad
        return dict(rect=rect,depth=depth,kind=kind,error=error)

    def basis(self,a,b):
        key=(a,b)
        if key in self.cache:return self.cache[key]
        a,b=arb(a),arb(b);half=(b-a)/2;mid=(b+a)/2
        x=[mid+half*t for t,w in self.nodes]
        weights=[half*w for t,w in self.nodes]
        h=[[self.pi**(-arb(1)/4)*(-t*t/2).exp() for t in x]]
        if self.n>1:h.append([arb(2).sqrt()*t*v for t,v in zip(x,h[0])])
        for k in range(1,self.n-1):
            h.append([(arb(2)/(k+1)).sqrt()*xj*v-(arb(k)/(k+1)).sqrt()*old
                      for xj,v,old in zip(x,h[-1],h[-2])])
        matrix=acb_mat(self.n,self.order,[z for row in h for z in row])
        exp_p=[(acb(0,1)*self.scale*t).exp() for t in x]
        exp_m=[(acb(0,-1)*self.scale*t).exp() for t in x]
        radial=[[(self.scale*t-(2*j+1)*self.pi).exp() for t in x] for j in range(self.terms)]
        result=(matrix,h[0],weights,exp_p,exp_m,radial)
        self.cache[key]=result
        return result

    def quadrature(self,rect):
        a,b,c,d=rect
        hx,vx,wx,px,mx,rx=self.basis(a,b)
        hy,vy,wy,py,my,ry=self.basis(c,d)
        w=hx.transpose()*self.C*hy
        total=acb(0)
        for i in range(self.order):
            for j in range(self.order):
                phase=acb(1)
                if not self.unit:
                    for k in range(self.terms):
                        phase*=(1+rx[k][i]*my[j])/(1+rx[k][i]*py[j])
                        phase*=(1+ry[k][j]*mx[i])/(1+ry[k][j]*px[i])
                total+=wx[i]*wy[j]*w[i,j]*vx[i]*vy[j]*phase
        assert total.is_finite()
        return total

    def partition(self,L,tolerance,progress=None,max_leaves=30000):
        first=self.bound((-float(L),float(L),-float(L),float(L)))
        queue=[(-float(first['error']),0,first)];counter=0;total=first['error']
        while not total < tolerance:
            if len(queue)>=max_leaves:raise RuntimeError('Partition budget exhausted; no certificate')
            _,_,old=heapq.heappop(queue)
            a,b,c,d=old['rect'];x=(a+b)/2;y=(c+d)/2
            if old['depth']>=40:raise RuntimeError('Dyadic depth limit; no certificate')
            total-=old['error']
            for rect in [(a,x,c,y),(x,b,c,y),(a,x,y,d),(x,b,y,d)]:
                new=self.bound(rect,old['depth']+1);counter+=1
                heapq.heappush(queue,(-float(new['error']),counter,new));total+=new['error']
            if progress and counter%200==0:progress(len(queue),total)
        cells=[r[2] for r in queue]
        # Recompute from final leaves; interval subtraction only guided refinement.
        actual=sum((c['error'] for c in cells),arb(0))
        assert actual < tolerance
        return cells,actual
