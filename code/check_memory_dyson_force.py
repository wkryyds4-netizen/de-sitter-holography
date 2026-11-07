"""Noncommuting 2x2 Gaussian fixture for exact Wick words and finite Dyson error."""
from pathlib import Path
from itertools import product
import hashlib,json,platform
import numpy as np
from numpy.polynomial.hermite_e import hermegauss
from numpy.polynomial.legendre import leggauss

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-uniform-bound'
I=np.eye(2,dtype=complex)
P=np.array([I,[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
B=P/2
H=np.array([0.,.4,0.,.7])


def exponent(components,s):
    components=np.asarray(components)
    r=np.linalg.norm(components[...,1:],axis=-1)
    vector=np.einsum('...a,aij->...ij',components[...,1:],P[1:])
    return np.exp(1j*s*components[...,0])[...,None,None]*(np.cos(s*r)[...,None,None]*I+
           1j*s*np.sinc(s*r/np.pi)[...,None,None]*vector)


def gaussian(order):
    x,w=hermegauss(order);w=w/np.sqrt(2*np.pi)
    indices=np.array(list(product(range(order),repeat=4)))
    g=x[indices];weights=np.prod(w[indices],axis=1)
    return g/2,weights,np.einsum('na,aij->nij',g/2,P)


def phi(x):return sum(b@x@b for b in B)


def wick4(A,D,E):
    crossing=sum(b@A@c@D@b@E@c for b in B for c in B)
    return phi(A)@D@phi(E)+phi(A@phi(D)@E)+crossing


def integral(alpha,s,order):
    x,w=leggauss(order);x=(x+1)/2;w=w/2
    result=np.zeros((2,2),complex)
    for a,wa in zip(x,w):
        for b,wb in zip(x,w):
            result+=s*s*a*wa*wb*wick4(exponent(alpha*H,s*(1-a)),
                                     exponent(alpha*H,s*a*(1-b)),exponent(alpha*H,s*a*b))
    return result


def hs(x):return float(np.linalg.norm(x)/np.sqrt(2))


if __name__=='__main__':
    dest=OUT/'dyson-force-checks.json'
    if dest.exists():raise FileExistsError(dest)
    assert np.linalg.norm(sum(b@b for b in B)-I)<1e-14
    _,w,u=gaussian(3);words=[]
    for ts in [(.2,-.5,.7),(-.8,.4,-.1),(1.3,.9,-.6)]:
        A,D,E=[exponent(H,t) for t in ts]
        direct=np.einsum('n,nij->ij',w,u@A@u@D@u@E@u)
        err=hs(direct-wick4(A,D,E));assert err<2e-12
        words.append(dict(times=ts,complete_wick_word_error=err))
    quadratures={order:gaussian(order) for order in (10,14)}
    rows=[]
    for eps in (.04,.16,.36):
        alpha=np.sqrt(1-eps)
        for s in (-.7,.7,1.4):
            vals=[]
            for order in (10,14):
                comps,w,u=quadratures[order]
                ex=exponent(alpha*H+np.sqrt(eps)*comps,s)
                vals.append(eps*np.einsum('n,nij->ij',w,u@ex@u))
            second=integral(alpha,s,18);second_lo=integral(alpha,s,12)
            proxy=eps*phi(exponent(alpha*H,s))-eps**2*second
            residual=hs(vals[1]-proxy)
            analytic_bound=eps**3*abs(s)**4*np.sqrt(10395)/24
            qdiff=hs(vals[1]-vals[0]);tdiff=hs(second-second_lo)
            assert qdiff<2e-11 and tdiff<2e-12
            assert residual<=analytic_bound+2e-11
            assert np.linalg.norm(second,2)<=1.5*s*s+1e-12
            rows.append(dict(epsilon=eps,alpha=alpha,time=s,remainder_norm=residual,
                             remainder_over_epsilon_cubed=residual/eps**3,
                             analytic_remainder_bound=analytic_bound,
                             Gaussian_order_difference=qdiff,time_order_difference=tdiff))
    result=dict(fixture='Full four-coordinate Hermitian 2x2 Gaussian series, B=Pauli/2, C=1; deterministic H=.4 sigma_x+.7 sigma_z. A noncommuting identity/expansion fixture, not a DS simulation.',
                exact_degree_four_quadrature_words=words,finite_time_rows=rows,
                scope='Three complete fourth-moment Wick comparisons and nine finite-real-time remainder checks. Quadrature order comparisons are reproducibility checks, not interval-certified quadrature errors. Uniform remainder is proved analytically.',
                environment=dict(python=platform.python_version(),numpy=np.__version__),
                input_sha256={'code/check_memory_dyson_force.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    with dest.open('x',encoding='utf-8') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='passed',Wick_words=3,finite_time_rows=len(rows),max_Gaussian_order_difference=max(x['Gaussian_order_difference'] for x in rows))))
