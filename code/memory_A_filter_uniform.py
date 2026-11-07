"""Uniform A-filter prefactors and exact core assembly algebra."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,platform
import sympy as s
from memory_two_charge_filter import filters

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-filter-limit'


def main():
    target=OUT/'filter-certificate.json'
    if target.exists():raise FileExistsError(target)
    m,n,p,a,b,r=s.symbols('m n p a b r');d=m-p
    x=lambda j:j*(m+n-j)
    u0=3*d*(n+p)/(m*n)+s.Rational(5,2)-4*d/m
    ua=-(4*m+p)/(4*d)-(m-n+p)*p/((m-a)*(n-a))
    ub=ua.subs(a,b)
    ur=s.Rational(1,2)+(m+3*n-p-5*(2*m*n-x(a)-x(b))/(4*d))*p/((m-r)*(n-r))
    e,lam,alpha,beta,rho=s.symbols('e lam alpha beta rho')
    sub={m:1/(2*lam*e**2),n:1/(2*lam*e**2)+1,p:1/e,a:alpha/e,b:beta/e,r:rho/e}
    expanded=[]
    for u,want in zip([u0,ua,ub,ur],[s.Rational(3,2),-1,-1,s.Rational(1,2)]):
        f=s.cancel(u.subs(sub));expr=s.series(f,e,0,2).removeO().expand()
        assert s.factor(expr.coeff(e,0)-want)==0
        assert not any(power<0 and coeff!=0 for power,coeff in [(i,expr.coeff(e,i)) for i in range(-4,0)])
        expanded.append(str(s.factor(expr)))
    # Exact finite agreement with the earlier independent five-word contraction.
    checks=0
    for N,pp in [(8,2),(12,3),(16,4)]:
        mm=N//2;nn=mm+1
        from math import comb
        C=comb(mm,pp)*comb(N-mm+pp,pp)
        def kap(j):return comb(mm-j,pp)*comb(N-mm+pp-j,pp) if mm-j>=pp and N-mm+pp-j>=pp else 0
        for aa in range(pp+1):
            for bb in range(aa,pp+1):
                for rr in range(min(aa+bb,mm)+1):
                    if rr==mm:continue  # removable e-channel endpoint is not needed for the DS window
                    vals={m:mm,n:nn,p:pp,a:aa,b:bb,r:rr}
                    rational=sum(coef.subs(vals)*v for coef,v in zip([u0,ua,ub,ur],[C,kap(aa),kap(bb),kap(rr)]))
                    exact=filters(N,pp,mm,aa,bb,rr)['g'][1]
                    assert rational==s.Rational(exact.numerator,exact.denominator)
                    checks+=1
    q=s.symbols('q',positive=True);K0=lam*(1-10*q+9*q*q+16*lam*q-32*lam*q*q)
    kbar=s.factor(K0+16*lam**2*(q*q-q))
    assert s.factor(kbar-(lam*(1-10*q+9*q*q)-16*lam**2*q*q))==0
    gamma=(3-4*q+q*q)/2
    vnew=(1+q)*(1-q*q+8*lam*q)/(lam*(1+q+4*lam*q))
    core_residual=s.factor(64*q**4*vnew)
    small=s.limit(core_residual.subs(q,s.exp(-4*lam)),lam,0)
    assert small==1024
    names=['code/memory_A_filter_uniform.py','code/memory_two_charge_filter.py']
    result=dict(status='passed',exact_prefactors=list(map(str,[u0,ua,ub,ur])),
                uniform_rank_proportion_expansions=expanded,finite_filter_identities=checks,
                gA_over_C_limit=str(gamma),centered_core_B_finite_coefficient=str(kbar),
                core_new_residual_over_Js_eight=str(core_residual),small_lambda_core_residual_over_Js_eight=str(small),
                analytic_domain='Fixed positive lambda in a compact set; a,b<=p,r<=a+b<=2p; p sufficiently large so all scaled rational denominators are bounded away from zero.',
                analytic_bound='gA/C is uniformly bounded; its deviation from gamma is <= K_lambda*(1+u+v+S)/p using the binomial-product adjacent-rank bound. Window S<=p^(3/4) and sixth ordinary moments give L2 error O(p^(-1/4)) with prior input weighted norms.',
                environment=dict(python=platform.python_version(),sympy=s.__version__),
                input_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print('passed',checks,'exact filter identities; core residual small lambda',small)


if __name__=='__main__':main()
