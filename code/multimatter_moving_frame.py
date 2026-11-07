"""Pure ordered two-matter triangular frame. No runs or writes on import."""
import numpy as np
from chord_module_frame import qi

def norm(v):
    return float(np.sqrt(np.sum(np.abs(v)**2)))

def shift(v,axis,coeff,up=False):
    out=np.zeros_like(v,dtype=complex)
    src=[slice(None)]*v.ndim; dst=src.copy()
    if up: src[axis]=slice(None,-1);dst[axis]=slice(1,None)
    else: src[axis]=slice(1,None);dst[axis]=slice(None,-1)
    out[tuple(dst)]=(v*coeff)[tuple(src)]
    return out

def action(eps,ds,name,v,limit=False):
    n=np.arange(v.shape[0])[:,None,None]
    k=np.arange(v.shape[1])[None,:,None]
    l=np.arange(v.shape[2])[None,None,:]
    w=ds[0]+ds[1]+k+l
    if limit:
        if name=='E':return w*v
        if name=='Q':return action(eps,ds,'E',v,True)-action(eps,ds,'P',v,True)
        out=np.zeros_like(v,dtype=complex)
        for axis,z,d in ((1,k,ds[0]),(2,l,ds[1])):
            up=shift(v,axis,np.sqrt((z+1)*(2*d+z)),True)
            down=shift(v,axis,np.sqrt(z*(2*d+z-1)))
            if name in ('B','T'):out+=(up+down)/2*(1 if name=='B' else -1)
            elif name=='P':out+=1j*(up-down)/2
            else:raise ValueError(name)
        return out
    if name=='Q':return action(eps,ds,'E',v)-action(eps,ds,'P',v)
    R=np.exp(-eps*w);omR=-np.expm1(-eps*w)
    R2=np.exp(-eps*(ds[1]+l))
    c=np.exp(-eps*(n+w)/2);root=np.sqrt(-np.expm1(-eps*(n+w)))
    def a(x,up=False):return shift(x,0,np.sqrt(qi(eps,n+(1 if up else 0))),up)
    def C(x,up=False):
        out=np.zeros_like(x,dtype=complex)
        for axis,z,d,fac in ((1,k,ds[0],R2),(2,l,ds[1],1)):
            zt=z+(1 if up else 0)
            coeff=np.sqrt(qi(eps,zt)*(-np.expm1(-eps*(2*d+zt-1))))
            out+=shift(x*fac,axis,coeff,up)
        return out
    if name=='T':return (omR*(a(v)+a(v,True))-C(v)-C(v,True))/(2*np.sqrt(eps))
    if name=='P':return 1j*(a(C(v),True)-C(a(v),True))/(2*root)
    cross=a(R*C(v),True)+C(R*a(v),True)+C(C(v),True)
    if name=='E':return ((qi(eps,w)-omR**2*qi(eps,n)/2)*v-cross/2)/c
    if name=='B':return ((-np.expm1(-2*eps*w))*qi(eps,n)*v-cross)/(2*c*root)
    raise ValueError(name)

def word(eps,ds,letters,v,limit=False):
    for name in reversed(letters):v=action(eps,ds,name,v,limit)
    return v

def exact_ladder(q,ds,vec,segment,up=False):
    """Monic orthogonal basis; A_j=a0 prod R + sum bi prod R."""
    out={}
    def qint(z):return (1-q**z)/(1-q)
    for state,c in vec.items():
        for i in range(segment+1):
            z=state[i]
            if not up and z==0:continue
            coeff=q**sum(ds[t-1]+state[t] for t in range(i+1,segment+1))
            if not up:
                coeff*=qint(z)
                if i:coeff*=1-q**(2*ds[i-1]+z-1)
            target=list(state);target[i]+=1 if up else -1;target=tuple(target)
            out[target]=out.get(target,0)+c*coeff
    return {u:c for u,c in out.items() if c}

def raw_image(q,ds,state):
    v={(0,0,0):1}
    for segment,count in enumerate(state):
        for _ in range(count):v=exact_ladder(q,ds,v,segment,True)
    return v

def exact_inner(q,ds,u,v):
    out=0
    for z,a in u.items():
        if z not in v:continue
        weight=1
        for i,number in enumerate(z):
            for j in range(1,number+1):
                weight*=(1-q**j)/(1-q)
                if i:weight*=1-q**(2*ds[i-1]+j-1)
        out+=a*v[z]*weight
    return out
