"""Exact rank-product probability and DS bounds for output-rank deficit.

The covariance maps of independent isotropic rank channels compose.
Racah orthogonality inverts their product without selecting multiplicities.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import hashlib,json,platform,time
import sympy as s
from cubic_linear_dynamics import crossing

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-rank-blocks'


def product_distribution(N,m,a,b):
    R=min(m,N-m);D=comb(N,m)
    dims=[comb(N,j)**2-(comb(N,j-1)**2 if j else 0) for j in range(R+1)]
    chi=[[crossing(N,m,j,k) for k in range(R+1)] for j in range(R+1)]
    for j in range(R+1):
        for k in range(R+1):
            expect=F(D*D,dims[j]) if j==k else F(0)
            assert sum(dims[l]*chi[j][l]*chi[k][l] for l in range(R+1))==expect
    prob=[F(dims[r],D*D)*sum(dims[l]*chi[a][l]*chi[b][l]*chi[r][l] for l in range(R+1)) for r in range(R+1)]
    assert sum(prob)==1 and min(prob)>=0
    assert all(prob[r]==0 for r in range(a+b+1,R+1))
    return prob


def exact_moments(order=3):
    N,m,x,y,z=s.symbols('N m x y z')
    P=[s.Integer(1)]
    for j in range(1,order+1):P.append(s.expand(P[-1]*(z-(j-1)*(N+2-j))))
    es=[s.Integer(1)];rows=[]
    for n in range(1,order+1):
        cs=[s.cancel(s.rf(-n,j)*s.rf(n-N-1,j)/(s.rf(-m,j)*s.rf(m-N,j)*s.rf(-N-1,j)*s.factorial(j))) for j in range(n+1)]
        cx=sum(cs[j]*P[j].subs(z,x) for j in range(n+1))
        cy=sum(cs[j]*P[j].subs(z,y) for j in range(n+1))
        value=s.cancel((cx*cy-sum(cs[j]*es[j] for j in range(n)))/cs[n]);es.append(value)
        print('Exact factorial Casimir moment',n,'ready',flush=True)
    return (N,m,x,y,z),P,es


def main():
    OUT.mkdir(exist_ok=True);target=OUT/'product-rank-moments.json'
    if target.exists():raise FileExistsError(target)
    start=time.time();symbols,P,es=exact_moments(3);N,m,x,y,z=symbols
    checks=[]
    for nn,mm,aa,bb in [(5,2,1,2),(6,3,2,2),(8,4,2,3),(10,5,2,2)]:
        prob=product_distribution(nn,mm,aa,bb);subs={N:nn,m:mm,x:aa*(nn+1-aa),y:bb*(nn+1-bb)}
        for j in range(min(3,mm,nn-mm)+1):
            empirical=sum(s.Rational(pr.numerator,pr.denominator)*P[j].subs({N:nn,z:r*(nn+1-r)}) for r,pr in enumerate(prob))
            assert s.simplify(empirical-es[j].subs(subs))==0
        checks.append(dict(N=nn,m=mm,a=aa,b=bb,checked_moment_order=min(3,mm,nn-mm),probabilities=list(map(str,prob))))
    # Delta x = x_(a+b)-x_r. Derive exact moments by the triangular P_j basis.
    top=s.symbols('top');delta=[]
    for k in range(1,4):
        expr=s.Poly((top-z)**k,z);weights=[s.Integer(0)]*(k+1)
        for j in range(k,-1,-1):
            weights[j]=expr.coeff_monomial(z**j);expr=s.Poly(s.expand(expr.as_expr()-weights[j]*P[j]),z)
        assert expr.is_zero
        delta.append(s.cancel(sum(weights[j]*es[j] for j in range(k+1))))
    # Fixed rank deficits u,v, half filling; general lambda is kept symbolic.
    p,lam,u,v=s.symbols('p lam u v');nn=p*p/lam;aa=p-u;bb=p-v
    replace={N:nn,m:nn/2,x:aa*(nn+1-aa),y:bb*(nn+1-bb),top:(aa+bb)*(nn+1-aa-bb)}
    limits=[]
    for k,expr in enumerate(delta,1):
        scaled=s.cancel(expr.subs(replace,simultaneous=True)/nn**k)
        numer,denom=scaled.as_numer_denom();num=s.Poly(numer,p);den=s.Poly(denom,p)
        assert num.degree()<=den.degree()
        limit=s.factor(num.LC()/den.LC()) if num.degree()==den.degree() else s.Integer(0)
        limits.append(str(limit));print('DS delta-Casimir moment',k,limit,flush=True)
    expected=[2*lam,2*lam+4*lam**2,2*lam+12*lam**2+8*lam**3]
    assert all(s.simplify(s.sympify(value)-ex)==0 for value,ex in zip(limits,expected))
    names=['code/rank_product_moments.py','code/cubic_linear_dynamics.py']
    result=dict(status='passed',finite_checks=checks,factorial_polynomials=list(map(str,P)),
                exact_factorial_expectations=list(map(str,es)),exact_delta_casimir_moments=list(map(str,delta)),
                ds_scaled_delta_casimir_moments=limits,
                scope='Independent isotropic rank-a/rank-b matrix product; fixed input deficits a=p-u,b=p-v; N=p^2/lambda,m=N/2. Moments through third do not determine a full output-deficit law.',
                elapsed_seconds=time.time()-start,environment=dict(python=platform.python_version(),sympy=s.__version__),
                input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
