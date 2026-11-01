"""Literal KS scalar multiplier and two finite-dimensional realizations.

Numerical pentagon residuals are diagnostics, NOT a counterexample to a
continuum theorem. The finite representations do not satisfy exact Weyl CCR.
"""
import numpy as np
from scipy.fft import fft2, ifft2, fftfreq
from scipy.linalg import eigh_tridiagonal

def phase(x, y, h, cutoff=38.0):
    x,y=np.broadcast_arrays(x,y)
    if h<=0: raise ValueError('h must be positive')
    arg=np.zeros(x.shape);tail=0.;terms=[]
    for a,b,r,sign in [(x,h*y,np.pi*h,-1),(y,-x/h,np.pi/h,1)]:
        n=max(1,int(np.ceil((float(np.max(a))+cutoff)/(2*r))))
        for j in range(n):
            v=a-(2*j+1)*r;e=np.exp(-abs(v))
            z=np.where(v>0,np.exp(1j*b)+e,1+e*np.exp(1j*b))
            arg+=sign*np.angle(z)
        first=np.exp(float(np.max(a))-(2*n+1)*r)
        tail+=first/((1-np.exp(-2*r))*(1-first));terms.append(n)
    return np.exp(2j*arg),dict(terms=terms,phase_tail_bound=2*tail,
        caveat='Exact-arithmetic product tail only; no rounding or spectral truncation bound.')

def multiplier(x,y,h,flat=False):
    if flat:return np.exp(-1j*y*np.logaddexp(0,x)/np.pi),dict(terms=[],phase_tail_bound=0.)
    return phase(x,y,h)

def measure(a,B,chirp,f):
    ab=a*B(f);ba=B(a*f);rhs=B(chirp.conj()*B(chirp*a*f))
    z=np.vdot(ab,ba);err=ab-rhs;pr=np.vdot(f,err)
    gram=np.array([[np.vdot(ab,ab),z],[z.conjugate(),np.vdot(ba,ba)]])
    return dict(norms=[float(np.linalg.norm(v)) for v in [f,ab,ba,rhs]],
        overlap=[float(z.real),float(z.imag)],D=float(np.linalg.norm(ab-ba)**2),
        gram_eigenvalues=np.linalg.eigvalsh(gram).tolist(),
        pentagon_residual=float(np.linalg.norm(err)),
        projected_residual=[float(pr.real),float(pr.imag)],
        residual_status='unresolved continuum reproduction; no pass threshold inferred')

def fft_run(n,L,h,flat=False):
    t=np.linspace(-L,L,n,endpoint=False);p=2*np.pi**2*fftfreq(n,d=2*L/n)
    x=t[:,None];y=t[None,:]
    a,ta=multiplier(x,y,h,flat);b,tb=multiplier(p[None,:],p[:,None],h,flat)
    chirp=np.exp(1j*x*y/np.pi)
    def B(f):return ifft2(b*fft2(f,norm='ortho'),norm='ortho')
    f=np.exp(-(x*x+y*y)/(2*np.pi));f/=np.linalg.norm(f)
    out=measure(a,B,chirp,f)
    out.update(method='periodic_fft',n=n,L=L,h=h,flat=flat,product_tails=[ta,tb],
        state='normalized sampled exp[-(x^2+y^2)/(2*pi)]',
        state_boundary_mass=float(np.sum(abs(f[(abs(x)+abs(y))>L/2])**2)))
    return out

def hermite_run(n,h,flat=False):
    z,Q=eigh_tridiagonal(np.zeros(n),np.sqrt(np.arange(1,n)/2))
    T=(Q.T*((-1j)**np.arange(n)))@Q
    t=np.sqrt(np.pi)*z;x=t[:,None];y=t[None,:]
    a,ta=multiplier(x,y,h,flat);b,tb=multiplier(y,x,h,flat)
    chirp=np.exp(1j*x*y/np.pi)
    def B(f):return T.conj().T@(b*(T@f@T.T))@T.conj()
    f=Q[0,:,None]*Q[0,None,:]
    out=measure(a,B,chirp,f)
    out.update(method='hermite_dvr',n=n,h=h,flat=flat,product_tails=[ta,tb],
        state='oscillator ground state with x=sqrt(pi)*z',
        fourier_unitarity_max=float(np.max(abs(T.conj().T@T-np.eye(n)))))
    return out
