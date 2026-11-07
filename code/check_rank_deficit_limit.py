"""Uniform input-rank scaling and independent high-precision deficit bounds."""
from pathlib import Path
import hashlib,json,platform,time
import sympy as s
import mpmath as mp
from memory_rank_blocks import ROOT,OUT


def main():
    target=OUT/'deficit-limit-checks.json'
    if target.exists():raise FileExistsError(target)
    source=OUT/'product-rank-moments.json';old=json.loads(source.read_text())
    N,m,x,y,top=s.symbols('N m x y top');p,lam,alpha,beta,rho=s.symbols('p lam alpha beta rho')
    expr=[s.sympify(v,locals={str(v):v for v in [N,m,x,y,top]}) for v in old['exact_delta_casimir_moments']]
    nn=p*p/lam;aa=alpha*p;bb=beta*p
    replace={N:nn,m:rho*nn,x:aa*(nn+1-aa),y:bb*(nn+1-bb),top:(aa+bb)*(nn+1-aa-bb)}
    limits=[];degrees=[]
    theta=s.factor(alpha*beta*lam*(1/(rho*(1-rho))-2))
    for k,v in enumerate(expr,1):
        rational=s.cancel(v.subs(replace,simultaneous=True)/nn**k)
        a,b=rational.as_numer_denom();a=s.Poly(a,p);b=s.Poly(b,p)
        assert a.degree()<=b.degree()
        limit=s.factor(a.LC()/b.LC()) if a.degree()==b.degree() else s.Integer(0)
        limits.append(limit);degrees.append([a.degree(),b.degree()])
        print('Uniform filling/rank moment',k,limit,flush=True)
    expected=[theta,theta**2+theta,theta**3+3*theta**2+theta]
    assert all(s.factor(a-b)==0 for a,b in zip(limits,expected))
    # Direct positive rational probability is checked independently in the parent.
    # Here evaluate exact moment expressions with two precisions, including
    # finite-p deviations and a rigorous-in-exact-arithmetic Markov bound.
    functions=[s.lambdify((N,m,x,y,top),v,'mpmath') for v in expr]
    rows=[]
    for pp in [16,32,64,128]:
        for deficits in [(0,0),(1,2)]:
            results=[]
            for dps in [80,120]:
                mp.mp.dps=dps;ll=mp.mpf('0.25');n=4*pp*pp;mm=n//2;a=pp-deficits[0];b=pp-deficits[1];t=a+b
                args=list(map(mp.mpf,[n,mm,a*(n+1-a),b*(n+1-b),t*(n+1-t)]))
                moments=[f(*args) for f in functions];A=n+1-2*t
                assert A>0 and min(moments)>=0
                # D <= Delta/A; each E D^j is bounded by these moments.
                upper=[moments[j]/mp.mpf(A)**(j+1) for j in range(3)]
                results.append([mp.nstr(z,65) for z in upper])
            assert max(abs(mp.mpf(a)/mp.mpf(b)-1) for a,b in zip(*results))<mp.mpf('1e-60')
            rows.append(dict(p=pp,N=n,m=mm,input_deficits=list(deficits),lambda_value='0.25',
                             D_moment_upper_bounds=results[1],limit_first_two_moments=['0.5','0.75'],
                             second_moment_relative_bound_excess=mp.nstr(mp.mpf(results[1][1])/mp.mpf('0.75')-1,40)))
    names=['code/check_rank_deficit_limit.py','results/susskind-memory-rank-blocks/product-rank-moments.json']
    target.write_text(json.dumps(dict(status='passed',uniform_theta=str(theta),
        uniform_scaled_delta_moment_limits=list(map(str,limits)),rational_degrees=degrees,records=rows,
        precisions=[80,120],environment=dict(python=platform.python_version(),sympy=s.__version__,mpmath=mp.__version__),
        scope='Uniform alpha,beta in [0,1], filling rho in a compact subset of (0,1), fixed positive lambda. Third delta moment bounds D^3, hence first two actual deficit moments converge. Not a proof of the complete Poisson law.',
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')


if __name__=='__main__':main()
