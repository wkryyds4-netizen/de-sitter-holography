"""Exploratory matrix elements of the fixed bounded heat-regulated word.

Neither FFT periodization nor Hermite quadrature is certified here.
Heat contractions must not be treated as unitary projection leakage.
"""
import numpy as np
from scipy.fft import fft2,ifft2,fftfreq
from scipy.linalg import eigh_tridiagonal
from unitary_flip import multiplier

def pack(z):return [float(z.real),float(z.imag)]

def fft_family(n,L,h,epsilons,flat=False):
    x=np.linspace(-L,L,n,endpoint=False)[:,None];y=x.T
    p=2*np.pi**2*fftfreq(n,d=2*L/n)
    a,ta=multiplier(x,y,h,flat);b,tb=multiplier(p[None,:],p[:,None],h,flat)
    momentum_square=p[:,None]**2+p[None,:]**2
    chirp=np.exp(1j*x*y/np.pi)
    f=np.exp(-(x*x+y*y)/(2*np.pi));f/=np.linalg.norm(f)
    ff=fft2(f,norm='ortho');af=a*f;maf=chirp*af;fm=fft2(maf,norm='ortho')
    rows=[]
    for epsilon in epsilons:
        symbol=b*np.exp(-epsilon*momentum_square)
        left=a*ifft2(symbol*ff,norm='ortho')
        middle=chirp.conj()*ifft2(symbol*fm,norm='ortho')
        right=ifft2(symbol*fft2(middle,norm='ortho'),norm='ortho')
        lhs=np.vdot(f,left);rhs=np.vdot(f,right);defect=lhs-rhs
        rows.append(dict(method='periodic_fft',n=n,L=L,h=h,flat=flat,epsilon=epsilon,
            lhs=pack(lhs),rhs=pack(rhs),defect=pack(defect),defect_modulus=float(abs(defect)),
            word_vector_residual=float(np.linalg.norm(left-right)),
            left_norm=float(np.linalg.norm(left)),right_norm=float(np.linalg.norm(right)),
            product_tails=[ta,tb]))
    return rows

def hermite_family(n,q,h,epsilons,flat=False):
    nodes,Q=eigh_tridiagonal(np.zeros(q),np.sqrt(np.arange(1,q)/2));W=Q[:n]
    x=np.sqrt(np.pi)*nodes[:,None];y=x.T
    a,ta=multiplier(x,y,h,flat);b,tb=multiplier(y,x,h,flat);chirp=np.exp(1j*x*y/np.pi)
    d=(-1j)**(np.arange(n)[:,None]+np.arange(n)[None,:])
    def mul(v,symbol):return W@(symbol*(W.T@v@W))@W.T
    f=np.zeros((n,n),complex);f[0,0]=1
    af=mul(f,a);maf=mul(af,chirp);rows=[]
    for epsilon in epsilons:
        symbol=b*np.exp(-epsilon*(x*x+y*y))
        def B(v):return d.conj()*mul(d*v,symbol)
        left=mul(B(f),a);right=B(mul(B(maf),chirp.conj()))
        lhs=left[0,0];rhs=right[0,0];defect=lhs-rhs
        rows.append(dict(method='oversampled_hermite',n=n,q=q,h=h,flat=flat,epsilon=epsilon,
            lhs=pack(lhs),rhs=pack(rhs),defect=pack(defect),defect_modulus=float(abs(defect)),
            word_vector_residual=float(np.linalg.norm(left-right)),
            left_norm=float(np.linalg.norm(left)),right_norm=float(np.linalg.norm(right)),
            product_tails=[ta,tb]))
    return rows
