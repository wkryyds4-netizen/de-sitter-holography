"""Pure finite CAR/Hermite and block-convolution checks for sector returns."""
from math import comb, factorial, sqrt
import numpy as np
from scipy.linalg import expm
from complex_syk_dictionary import fermions, creators


def system(n, p, m):
    cs = fermions(n)
    creation = creators(cs, p)
    scale = comb(m,p)*comb(n-m+p,p)
    full = []
    for a, ca in enumerate(creation):
        full.append(ca@ca.conj().T/sqrt(scale))
        for cb in creation[:a]:
            x = ca@cb.conj().T
            full.extend([(x+x.conj().T)/sqrt(2*scale),
                         1j*(x-x.conj().T)/sqrt(2*scale)])
    hi = [s for s in range(1 << n) if s.bit_count()==m]
    lo = [s for s in range(1 << n) if s.bit_count()==m-1]
    upper = [b[np.ix_(hi,hi)] for b in full]
    lower = [b[np.ix_(lo,lo)] for b in full]
    c = [a[np.ix_(lo,hi)] for a in cs]
    return upper, lower, c


def add(a,b,factor=1):
    out={k:x.copy() for k,x in a.items()}
    for k,x in b.items():
        out[k]=out.get(k,0)+factor*x
    return out


def step(poly, linear, right=False):
    """Multiply by a linear Gaussian field in probabilists' Hermite basis."""
    out={}
    for key,x in poly.items():
        for j,b in enumerate(linear):
            y=x@b if right else b@x
            up=list(key);up[j]+=1;up=tuple(up)
            out[up]=out.get(up,0)+y
            if key[j]:
                down=list(key);down[j]-=1;down=tuple(down)
                out[down]=out.get(down,0)+key[j]*y
    return out


def norm(poly, dimension):
    return sqrt(sum(np.prod([factorial(k) for k in key],dtype=float)
                    *float(np.vdot(x,x).real)/dimension for key,x in poly.items()))


def inverse_rows(n,p,m,maximum):
    upper,lower,c=system(n,p,m)
    d=len(upper[0]);dl=len(lower[0]);zero=(0,)*len(upper)
    powers_u=[{zero:np.eye(d,dtype=complex)}]
    powers_l=[{zero:np.eye(dl,dtype=complex)}]
    for _ in range(maximum+2):
        powers_u.append(step(powers_u[-1],upper))
        powers_l.append(step(powers_l[-1],lower))
    K=[]
    for u,l in zip(powers_u,powers_l):
        out={key:(p*sum(key)-m)*x for key,x in u.items()}
        for key,x in l.items():
            out[key]=out.get(key,0)+sum(a.conj().T@x@a for a in c)
        K.append(out)
    assert norm(K[0],d)<1e-12 and norm(K[1],d)<1e-12
    rows=[]
    for k in range(maximum+1):
        gamma={}
        for a in c:
            force=[b@a-a@h for b,h in zip(lower,upper)]
            term=step(step(powers_l[k],[f.conj().T for f in force]),force,True)
            gamma=add(gamma,term)
        phi={key:sum(h@x@h for h in upper) for key,x in powers_u[k].items()}
        direct=add(gamma,phi,-2*p)
        inverse=add(add(add(K[k+2],step(K[k+1],upper),-1),
                        step(K[k+1],upper,True),-1),
                    step(step(K[k],upper),upper,True))
        error=norm(add(direct,inverse,-1),d)
        scale=max(1.,norm(direct,d),norm(inverse,d))
        low={key:x for key,x in direct.items() if sum(key)<k+2}
        assert error/scale<2e-12,(n,p,m,k,error,scale)
        rows.append(dict(degree=k,hermite_coefficients=len(direct),
                         direct_norm=norm(direct,d),lower_chaos_norm=norm(low,d),
                         absolute_error=error,relative_error=error/scale))
    return dict(N=n,p=p,m=m,C_m=1,coordinates=len(upper),maximum=maximum,rows=rows)


def schatten(x,r,d):
    return float((np.sum(np.linalg.svd(x,compute_uv=False)**r)/d)**(1/r))


def convolution(generators, bridges, t):
    """Exact block exponential realizes an ordered convolution (float64)."""
    sizes=[len(a) for a in generators]
    cuts=np.cumsum([0]+sizes)
    matrix=np.zeros((cuts[-1],cuts[-1]),complex)
    for j,a in enumerate(generators):
        matrix[cuts[j]:cuts[j+1],cuts[j]:cuts[j+1]]=1j*a
    for j,b in enumerate(bridges):
        matrix[cuts[j]:cuts[j+1],cuts[j+1]:cuts[j+2]]=b
    return expm(t*matrix)[:sizes[0],cuts[-2]:]


def dilation_rows(n,p,m,seed,times):
    upper,lower,c=system(n,p,m)
    rng=np.random.default_rng(seed)
    g=rng.normal(size=len(upper))
    h=sum(x*b for x,b in zip(g,upper)); low=sum(x*b for x,b in zip(g,lower))
    d=len(h);V=np.vstack(c)/sqrt(m);K=np.kron(np.eye(n),low)
    Q=np.eye(len(K))-V@V.conj().T
    a=(1-p/m)*h;R=Q@K@V;D=Q@K@Q
    DF=K@V-V@h;G=R.conj().T@R
    A=sum((low@ci-ci@h).conj().T@(low@ci-ci@h) for ci in c)
    structural={
        'isometry':np.linalg.norm(V.conj().T@V-np.eye(d)),
        'charge_compression':np.linalg.norm(V.conj().T@K@V-a),
        'force_shift':np.linalg.norm(R-DF-(p/m)*V@h),
        'force_gram':np.linalg.norm(G-A/m+(p/m)**2*h@h),
    }
    assert max(structural.values())<1e-11
    rows=[]
    for t in times:
        ek=expm(1j*t*K)
        lam=R.conj().T@ek@R;mem=R.conj().T@expm(1j*t*D)@R
        j1=R.conj().T@convolution([K,a,K],[R,R.conj().T],t)@R
        e1=R.conj().T@convolution([D,a,K,a,K],[R,R.conj().T,R,R.conj().T],t)@R
        residual=mem-lam-j1-e1
        gamma=sum((low@ci-ci@h).conj().T@expm(1j*t*low)@(low@ci-ci@h) for ci in c)
        ferr=np.linalg.norm(DF.conj().T@ek@DF-gamma/m)
        rel=np.linalg.norm(residual)/max(1,np.linalg.norm(mem),np.linalg.norm(lam))
        assert rel<2e-12 and ferr<2e-11
        bounds=[]
        for r in (1,2,4):
            actual_j=schatten(j1,r,d);bound_j=t*t/2*schatten(G,2*r,d)**2
            actual_e=schatten(e1,r,d);bound_e=t**4/24*schatten(G,3*r,d)**3
            assert actual_j<=bound_j+1e-11 and actual_e<=bound_e+1e-11
            bounds.append(dict(r=r,J1=actual_j,J1_bound=bound_j,E1=actual_e,E1_bound=bound_e))
        rows.append(dict(t=t,volterra_relative_error=rel,force_sandwich_error=float(ferr),bounds=bounds))
    return dict(N=n,p=p,m=m,C_m=1,seed=seed,source_dimension=d,
                dilation_dimension=len(K),structural_errors={k:float(v) for k,v in structural.items()},rows=rows)
