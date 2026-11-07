"""Pure finite and rational checks for charged OU and bilocal sector interfaces."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import sqrt
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import block_diag, null_space
import sector_memory_return as mr
from slater_joint_site_polynomial import pairings
from slater_composite_sources import source_exponent


def padd(a,b,factor=F(1)):
    out=dict(a)
    for k,v in b.items():out[k]=out.get(k,F(0))+factor*v
    return {k:v for k,v in out.items() if v}


def fock(word,q):
    state={():F(1)}
    for color in reversed(word):
        out=defaultdict(F)
        for tensor,v in state.items():
            out[(color,)+tensor]+=v
            for j,c in enumerate(tensor):
                if c==color:out[tensor[:j]+tensor[j+1:]]+=v*q**j
        state=dict(out)
    return state


@lru_cache(None)
def counts(word):
    def pairs(points):
        if not points:yield ();return
        a=points[0]
        for j,b in enumerate(points[1:],1):
            if word[a]==word[b]:
                for rest in pairs(points[1:j]+points[j+1:]):yield ((a,b),)+rest
    out=defaultdict(int)
    for pi in pairs(tuple(range(len(word)))):
        cross=sum(a<c<b<d or c<a<d<b for i,(a,b) in enumerate(pi) for c,d in pi[i+1:])
        out[cross]+=1
    return dict(out)


def moment(word,q):
    return sum((F(n)*q**k for k,n in counts(tuple(word)).items()),F(0))


def colored_checks():
    source_rows=[]
    for length in (2,4,6,8):
        for pi in pairings(tuple(range(length))):
            cross=sum(a<c<b<d or c<a<d<b for i,(a,b) in enumerate(pi) for c,d in pi[i+1:])
            src=source_exponent(pi,('G',)*len(pi),())
            assert src=={0:F(-4*cross)} or (not src and cross==0),(pi,src)
            source_rows.append(dict(length=length,crossings=cross,constant=int(src.get(0,0))))
    rows=[]
    blocks=[abc for abc in product(range(3),repeat=3) if sum(abc)<=2]
    def block(abc,color):
        a,b,c=abc
        return (color,)+(0,)*a+(color,)+(0,)*b+(color,)+(0,)*c+(color,)
    for q in (F(1,4),F(2,3)):
        colored=0
        for length in (0,2,4,6):
            for word in product(range(3),repeat=length):
                assert fock(word,q).get((),F(0))==moment(word,q)
                colored+=1
        hermites=[{0:F(1)},{1:F(1)}]
        for j in range(1,6):
            hermites.append(padd({k+1:v for k,v in hermites[-1].items()},hermites[-2],
                                -sum((q**i for i in range(j)),F(0))))
        projected={}
        for abc in blocks:
            poly={}
            for tensor,v in fock(block(abc,1),q).items():
                if all(c==0 for c in tensor):poly=padd(poly,hermites[len(tensor)],v)
            projected[abc]=poly
        replica=0
        for abc,xyz,gap in product(blocks,blocks,range(2)):
            word=tuple(reversed(block(abc,1)))+(0,)*gap+block(xyz,2)
            expected=sum((a*b*moment((0,)*(i+j+gap),q)
                          for i,a in projected[abc].items() for j,b in projected[xyz].items()),F(0))
            assert moment(word,q)==expected
            replica+=1
        rows.append(dict(q=str(q),colored_words=colored,whole_block_replica=replica))
    return dict(source_rows=source_rows,parameters=rows)


def evaluate(poly,g):
    maximum=max(sum(k) for k in poly)
    hs=[]
    for x in g:
        values=[1.,x]
        for j in range(1,maximum):values.append(x*values[-1]-j*values[-2])
        hs.append(values)
    return sum(np.prod([hs[i][v] for i,v in enumerate(key)])*x for key,x in poly.items())


def ou_sandwich(h,basis,alpha,n):
    beta=sqrt(1-alpha**2);z=(0,)*len(basis)
    polynomial={tuple(int(i==j) for i in range(len(basis))):b.copy() for j,b in enumerate(basis)}
    for _ in range(n):
        fixed={k:alpha*x@h for k,x in polynomial.items()}
        polynomial=mr.add(fixed,mr.step(polynomial,[beta*b for b in basis],True))
    return sum(polynomial.get(tuple(int(i==j) for i in range(len(basis))),0)@b
               for j,b in enumerate(basis))


def charged_rows(n,p,m,maximum,seed):
    upper,lower,c=mr.system(n,p,m);d=len(upper[0]);zero=(0,)*len(upper)
    g=np.random.default_rng(seed).normal(size=len(upper))
    h=sum(x*b for x,b in zip(g,upper));low=sum(x*b for x,b in zip(g,lower))
    alpha=1-p/m;eps=1-alpha**2
    powers=[{zero:np.eye(d,dtype=complex)}]
    for _ in range(maximum+2):powers.append(mr.step(powers[-1],upper))
    op=[evaluate({key:alpha**sum(key)*x for key,x in a.items()},g) for a in powers]
    transfer=[sum(ci.conj().T@np.linalg.matrix_power(low,k)@ci for ci in c)/m
              for k in range(maximum+3)]
    e=[m*(a-b) for a,b in zip(transfer,op)]
    V=np.vstack(c)/sqrt(m);K=np.kron(np.eye(n),low);R=(np.eye(len(K))-V@V.conj().T)@K@V
    rows=[]
    for k in range(maximum+1):
        ou=eps*ou_sandwich(h,upper,alpha,k)
        qou=op[k+2]-alpha*(h@op[k+1]+op[k+1]@h)+alpha**2*h@op[k]@h
        lhs=m*(R.conj().T@np.linalg.matrix_power(K,k)@R-ou)
        rhs=e[k+2]-alpha*(h@e[k+1]+e[k+1]@h)+alpha**2*h@e[k]@h
        errors=[np.linalg.norm(ou-qou)/max(1,np.linalg.norm(ou)),
                np.linalg.norm(lhs-rhs)/max(1,np.linalg.norm(lhs),np.linalg.norm(rhs))]
        assert max(errors)<3e-12,(n,p,m,k,errors)
        rows.append(dict(degree=k,OU_relative_error=float(errors[0]),charged_relative_error=float(errors[1]),
                         charged_norm=float(np.linalg.norm(lhs)/sqrt(d))))
    return dict(N=n,p=p,m=m,C_m=1,seed=seed,rows=rows)


def number_rows():
    rows=[]
    for p,x in product((2,4,8,16,32),(-1,0,1)):
        n=4*p*p;m=n//2+x*p;lam=F(p*p,n);alpha=1-F(p,m)
        a=4*lam*F(abs(m-n//2),n)+8*lam**2/(3*p)
        for k in (0,1,2,3,8,32,80):
            rem=m*(1-alpha**k)-p*k+lam*k*(k-1)
            bound=a*(k+1)**3
            assert abs(rem)<=bound
            rows.append(dict(N=n,p=p,m=m,k=k,remainder=str(rem),bound=str(bound)))
    return rows


class Propagator:
    def __init__(self,H,V):
        self.e,u=np.linalg.eigh(H);self.w=u.conj().T@V
    def __call__(self,t):
        return self.w.conj().T@(np.exp(1j*t*self.e)[:,None]*self.w)


def integrate(left,middle,right,t,order):
    nodes,weights=leggauss(order);nodes=(nodes+1)/2;weights=weights/2
    out=0
    for x,wx in zip(nodes,weights):
        l=left(t*(1-x))
        for y,wy in zip(nodes,weights):
            out=out+t*t*x*wx*wy*l@middle(t*x*(1-y))@right(t*x*y)
    return out


def bilocal_rows(kind):
    if kind=='CAR':
        n,p,m,seed=4,2,3,819
        up,lo,c=mr.system(n,p,m);g=np.random.default_rng(seed).normal(size=len(up))
        h=sum(x*b for x,b in zip(g,up));low=sum(x*b for x,b in zip(g,lo))
        V=np.vstack(c)/sqrt(m);K=np.kron(np.eye(n),low)
        bath=null_space(V.conj().T);R=bath.conj().T@K@V;D=bath.conj().T@K@bath
        alpha=1-p/m
    else:
        m=3;seed=820;alpha=.64
        h=np.diag([-.4,-.4,.7]).astype(complex)
        rng=np.random.default_rng(seed);raw=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
        D=(raw+raw.conj().T)/5
        R=.17*(rng.normal(size=(4,3))+1j*rng.normal(size=(4,3)))
        K=np.block([[alpha*h,R.conj().T],[R,D]])
        V=np.vstack([np.eye(3),np.zeros((4,3))])
    d=len(h);evals,unitary=np.linalg.eigh(h)
    groups=[]
    for i,e in enumerate(evals):
        if groups and abs(e-evals[groups[-1][0]])<1e-10:groups[-1].append(i)
        else:groups.append([i])
    projectors=[unitary[:,ids]@unitary[:,ids].conj().T for ids in groups]
    pin=lambda x:sum(P@x@P for P in projectors)
    D0=block_diag(*([D]*len(projectors)));R0=np.vstack([R@P for P in projectors])
    K0=np.block([[alpha*h,R0.conj().T],[R0,D0]])
    V0=np.vstack([np.eye(d),np.zeros((len(D0),d))])
    prop=Propagator(K,V);static=Propagator(K0,V0);memory=Propagator(D,R)
    phase=Propagator(h,np.eye(d));Q=lambda t:m*(memory(t)-pin(memory(t)))
    gram=R.conj().T@R
    structure=[]
    for t in (-1.,0.,.35,.9,1.):
        f=prop(t);f0=static(t)
        comm=np.linalg.norm(h@f0-f0@h)
        assert np.linalg.norm(f,2)<=1+1e-12 and np.linalg.norm(f0,2)<=1+1e-12 and comm<1e-11
        structure.append(dict(t=t,F_op=float(np.linalg.norm(f,2)),F0_op=float(np.linalg.norm(f0,2)),commutator=float(comm)))
    rows=[]
    for t in (-.55,.35,.9):
        center=phase(-t/2);B=m*center@prop(t)@center;B=B-pin(B)
        values=[]
        for order in (10,18):
            exact=integrate(static,Q,prop,t,order)
            filt=-center@integrate(phase,Q,phase,t,order)@center
            identity=m*(prop(t)-static(t))+exact
            values.append((filt,float(np.linalg.norm(identity))))
        delta=np.linalg.norm(values[1][0]-values[0][0]);rel=values[1][1]/max(1,np.linalg.norm(B))
        assert delta<2e-10 and rel<3e-11
        bounds=[]
        for r in (1,2,4):
            gn=mr.schatten(gram,2*r,d);hn=mr.schatten(h,2*r,d)
            bound=m*t**4/3*gn**2+4*m*abs(1-alpha)*abs(t)**3/3*gn*hn
            observed=mr.schatten(B-values[1][0],r,d)
            assert observed<=bound+1e-11
            bounds.append(dict(r=r,observed=observed,coarse_finite_bound=bound))
        rows.append(dict(t=t,quadrature_difference=float(delta),volterra_residual=rel,
                         bilocal_norm=mr.schatten(B,2,d),filter_pinching_norm=mr.schatten(pin(values[1][0]),2,d),bounds=bounds))
    degenerate_non_scalar=0.
    for ids in groups:
        if len(ids)>1:
            block=unitary[:,ids].conj().T@pin(memory(.35))@unitary[:,ids]
            degenerate_non_scalar=max(degenerate_non_scalar,float(np.linalg.norm(block-np.trace(block)/len(ids)*np.eye(len(ids)))))
    if kind!='CAR':assert degenerate_non_scalar>1e-4
    assert max(row['bilocal_norm'] for row in rows)>1e-8
    return dict(kind=kind,seed=seed,m=m,alpha=alpha,source_dimension=d,
                energy_block_sizes=[len(x) for x in groups],degenerate_block_non_scalar=degenerate_non_scalar,
                structure=structure,rows=rows)
