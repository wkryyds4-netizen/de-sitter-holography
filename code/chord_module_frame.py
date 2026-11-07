"""Pure orthogonal one-matter frame and moving coherent-core actions.

Parameters: epsilon=-log(q), r=q**d. n is the common left oscillator,
k labels its vacuum multiplicity. No output or scientific run on import.
"""
import math
import numpy as np
from scipy.special import gammaln


def qi(eps, x):
    return -np.expm1(-eps*np.asarray(x))/(-np.expm1(-eps))


def coherent_descendant(m, j, cutoff):
    """Exact coefficients of D(-sqrt(m))|j>, j=0,1,2, then truncated."""
    n=np.arange(cutoff+1)
    c=(-1.)**n*np.exp((-m+n*np.log(m)-gammaln(n+1))/2)
    if j==0:return c
    if j==1:return (m-n)/math.sqrt(m)*c
    if j==2:return ((n-m)**2-n)/(math.sqrt(2)*m)*c
    raise ValueError(j)


def act(eps, d, name, v, limit=False):
    """Operators in a rectangular orthonormal (n,k) grid.

    The caller pads support before applying words. T denotes Btr;
    Q=E-P. Coefficients use expm1 and a cancellation-free E diagonal.
    """
    if name=='Q':return act(eps,d,'E',v,limit)-act(eps,d,'P',v,limit)
    n=np.arange(v.shape[0])[:,None]
    k=np.arange(v.shape[1])[None,:]
    out=np.zeros_like(v,dtype=complex)
    if limit:
        up=np.sqrt((k+1)*(2*d+k))
        if name=='E':return (d+k)*v
        if name in ('B','T'):
            sign=1 if name=='B' else -1
            out[:,1:]+=sign*up[:,:-1]/2*v[:,:-1]
            out[:,:-1]+=sign*up[:,:-1]/2*v[:,1:]
        elif name=='P':
            out[:,1:]+=1j*up[:,:-1]/2*v[:,:-1]
            out[:,:-1]+=-1j*up[:,:-1]/2*v[:,1:]
        else:raise ValueError(name)
        return out
    nn=qi(eps,n);R=np.exp(-eps*(d+k));omR=-np.expm1(-eps*(d+k))
    c=np.exp(-eps*(n+k+d)/2)
    root=np.sqrt(-np.expm1(-eps*(n+k+d)))
    su=np.sqrt(qi(eps,k+1)*(-np.expm1(-eps*(2*d+k))))
    sd2=qi(eps,k)*(-np.expm1(-eps*(2*d+k-1)))
    if name=='T':
        a=omR*np.sqrt(qi(eps,n+1))/(2*np.sqrt(eps))
        out[1:,:]+=a[:-1,:]*v[:-1,:]
        out[:-1,:]+=a[:-1,:]*v[1:,:]
        b=su/(2*np.sqrt(eps))
        out[:,1:]-=b[:,:-1]*v[:,:-1]
        out[:,:-1]-=b[:,:-1]*v[:,1:]
        return out
    if name=='B':
        out=(((-np.expm1(-2*eps*(d+k)))*nn-sd2)/(2*c*root))*v
        forward=-R*su*np.sqrt(nn)/(2*c*root)
    elif name=='E':
        out=((qi(eps,d+k)-omR**2*nn/2-sd2/2)/c)*v
        forward=-R*su*np.sqrt(nn)/(2*c)
    elif name=='P':
        forward=-1j*su*np.sqrt(nn)/(2*root)
    else:raise ValueError(name)
    # (n,k)->(n-1,k+1), and the conjugate reverse edge.
    edge=forward[1:,:-1]
    out[:-1,1:]+=edge*v[1:,:-1]
    out[1:,:-1]+=np.conjugate(edge)*v[:-1,1:]
    return out


def word(eps,d,letters,v,limit=False):
    for name in reversed(letters):v=act(eps,d,name,v,limit)
    return v


def residuals(eps,d,v):
    out={}
    for w in ('B','P','E','T','Q','BE','PT','EBP','EEEE'):
        a=word(eps,d,w,v);b=word(eps,d,w,v,True)
        err=np.linalg.norm(a-b);scale=np.linalg.norm(b)
        out[w]=dict(absolute=float(err),target_norm=float(scale),relative=float(err/scale))
    for a,b,c in [('B','P','E'),('E','P','B'),('B','E','P')]:
        first=word(eps,d,a+b,v);second=word(eps,d,b+a,v);last=1j*act(eps,d,c,v)
        norm=np.linalg.norm(first-second-last)
        out[f'comm_{a}{b}']=dict(absolute=float(norm),
          scale=float(np.linalg.norm(last)),
          cancellation=float((np.linalg.norm(first)+np.linalg.norm(second)+np.linalg.norm(last))/norm))
    cas=word(eps,d,'EE',v)-word(eps,d,'BB',v)-word(eps,d,'PP',v)-d*(d-1)*v
    out['casimir']=dict(absolute=float(np.linalg.norm(cas)))
    out['boost_match']=dict(absolute=float(np.linalg.norm(act(eps,d,'B',v)+act(eps,d,'T',v))))
    return out


def raw_basis(q,r,n,k):
    """Unnormalized (aLdag)^n product_j(aRdag-r q^j aLdag)|seed>."""
    vec={(0,0):q*0+1}
    for j in range(k):
        out={}
        for (l,z),a in vec.items():
            out[l,z+1]=out.get((l,z+1),0)+a
            out[l+1,z]=out.get((l+1,z),0)-r*q**j*a
        vec=out
    return {(l+n,z):a for (l,z),a in vec.items()}


def raw_ladder(q,r,vec,side,creation=False):
    out={}
    for (l,k),a in vec.items():
        if creation:
            u=(l+(side=='L'),k+(side=='R'))
            out[u]=out.get(u,0)+a
        else:
            if l:
                u=(l-1,k);w=(1-q**l)/(1-q)*(1 if side=='L' else r*q**k)
                out[u]=out.get(u,0)+a*w
            if k:
                u=(l,k-1);w=(1-q**k)/(1-q)*(1 if side=='R' else r*q**l)
                out[u]=out.get(u,0)+a*w
    return {u:a for u,a in out.items() if a}
