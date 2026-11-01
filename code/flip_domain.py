"""Oversampled Hermite projection, not unitary nodewise DVR.

The resulting leakage estimates depend on quadrature. They are NOT
certified continuum bounds. Scalar multiplier is the frozen prior input.
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal
from unitary_flip import multiplier

def project_run(n,q,h,flat=False):
    if q<=n:raise ValueError('Use independent oversampling, q>n')
    nodes,Q=eigh_tridiagonal(np.zeros(q),np.sqrt(np.arange(1,q)/2))
    W=Q[:n];x=np.sqrt(np.pi)*nodes[:,None];y=x.T
    a,ta=multiplier(x,y,h,flat);b,tb=multiplier(y,x,h,flat)
    m=np.exp(1j*x*y/np.pi)
    d=(-1j)**(np.arange(n)[:,None]+np.arange(n)[None,:])
    def mul(v,f):return W@(f*(W.T@v@W))@W.T
    def A(v):return mul(v,a)
    def B(v):return d.conj()*mul(d*v,b)
    def M(v):return mul(v,m)
    def Mi(v):return mul(v,m.conj())
    f=np.zeros((n,n),complex);f[0,0]=1
    def chain(ops):
        vectors=[f.copy()];leaks=[];norms=[1.]
        for op in ops:
            v=op(vectors[-1]);no=float(np.linalg.norm(v));loss=norms[-1]**2-no**2
            leaks.append(float(np.sqrt(max(0,loss))));norms.append(no);vectors.append(v)
        return vectors,leaks,norms
    lv,ll,ln=chain([B,A]);rv,rl,rn=chain([A,M,B,Mi,B])
    delta=lv[-1]-rv[-1];r=float(np.linalg.norm(delta))
    result=dict(n=n,q=q,h=h,flat=flat,residual=r,left_leaks=ll,right_leaks=rl,
        left_norms=ln,right_norms=rn,quadrature_dependent_margin=r-sum(ll)-sum(rl),
        vacuum_projection=[float(delta[0,0].real),float(delta[0,0].imag)],
        row_orthogonality_max=float(np.max(abs(W@W.T-np.eye(n)))),product_tails=[ta,tb],
        state='normalized oscillator vacuum; physical x=sqrt(pi)*node',
        scope='Oversampled approximations to true Galerkin matrix elements. Quadrature errors not certified; positive margin is not a continuum lower bound.')
    return result,dict(left=np.stack(lv),right=np.stack(rv))
