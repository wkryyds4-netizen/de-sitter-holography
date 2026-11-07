"""Symbolic cancellation of p^2 and p terms in the coherent A variance."""
from pathlib import Path
import hashlib,json,platform
import sympy as s

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-residual-norm'


def main():
    target=OUT/'symbolic-limit.json'
    if target.exists():raise FileExistsError(target)
    j,k,u,v,lam,e=s.symbols('j k u v lam e');p=1/e;m=p*p/(2*lam);N=2*m;n=m+1
    a=p-u;b=p-v;x=lambda a:a*(N+1-a)
    g=(m-p)**2/m;h=2*p-m;eta=p*p/m
    P=g*g*x(a)*x(b)/((m*n-x(a))*(m*n-x(b)))
    cf=(m-1)*n*(N+1)/N
    Ps=s.series(P*e**2,e,0,3).removeO().expand()
    Ds=s.series((g*g/n-P/n+2*h*g)*j/(m-j)*e**2,e,0,3).removeO().expand()
    unsigned=s.factor(s.series(eta**2+P*(m*n/cf-1),e,0,1).removeO())
    # Verify the per-step logarithmic expansion from the exact term ratio.
    factor=s.cancel(((-a+k)*(a-N-1+k)*(-b+k)*(b-N-1+k))/((-m+k)**2*(-N-1+k)*(-4*lam)))
    step=s.series(factor,e,0,3).removeO().expand()
    assert step.coeff(e,0)==1
    l1=-u-v-2*k-2*lam
    l2=-((u+k)**2+(v+k)**2)/2+lam*(u+v+1+3*k)-lam**2
    assert s.factor(step.coeff(e,1)-l1)==0
    assert s.factor(step.coeff(e,2)-l1*l1/2-l2)==0
    L1=-j*(u+v+2*lam)-j*(j-1)
    L2=-j*(u*u+v*v)/2-(u+v)*j*(j-1)/2-j*(j-1)*(2*j-1)/6+lam*j*(u+v+1)+s.Rational(3,2)*lam*j*(j-1)-lam*lam*j
    assert s.factor(s.summation(l1,(k,0,j-1))-L1)==0
    assert s.factor(s.summation(l2,(k,0,j-1))-L2)==0
    expansion=s.Poly(s.expand((-Ps+Ds+eta**2*e**2)*(1+L1*e+(L2+L1**2/2)*e**2)),e)
    z=s.symbols('z');touchard=[s.Integer(1)]
    for degree in range(10):touchard.append(s.expand(z*(touchard[-1]+s.diff(touchard[-1],z))))
    def expectation(poly):
        poly=s.Poly(s.expand(poly),j)
        return s.factor(sum(c*touchard[power[0]].subs(z,-4*lam) for power,c in poly.terms()))
    averaged=[expectation(expansion.coeff_monomial(e**power)) for power in range(3)]
    assert averaged==[0,0,4*lam**2]
    assert unsigned==4*lam*(lam+1)
    q,theta=s.symbols('q theta',positive=True)
    zz=1+q;bb=theta*(1-q);ss=theta**2*zz+theta*(1-q);aa=2*theta+theta**2*zz
    G=s.Matrix([[zz,0,bb,bb],[0,theta*zz,0,0],[bb,0,ss,ss],[bb,0,ss,aa]])
    vs=s.factor(ss-bb*bb/zz);delta=s.factor(aa-ss)
    assert delta==theta*zz
    assert s.factor(G.det()-theta**2*zz**3*vs)==0
    ids=[0,1,3];schur=s.factor(G[2,2]-(G.extract([2],ids)*G.extract(ids,ids).inv()*G.extract(ids,[2]))[0])
    assert s.factor(schur-vs*delta/(vs+delta))==0
    physical=s.factor(schur.subs(theta,2*lam)/lam**2)
    small=s.limit(physical.subs(q,s.exp(-4*lam)),lam,0)
    assert small==16
    names=['code/memory_norm_asymptotic.py']
    result=dict(status='passed',unsigned_limit=str(unsigned),log_terms=[str(L1),str(L2)],
                P_scaled_series=[str(s.factor(Ps.coeff(e,power))) for power in range(3)],
                summand_coefficients=[str(s.factor(expansion.coeff_monomial(e**power))) for power in range(3)],
                averaged_coefficients=list(map(str,averaged)),
                norm_over_C_squared='4*lambda + 4*lambda^2*(1+q)',
                gram_over_C_squared=[[str(x) for x in row] for row in G.tolist()],
                determinant=str(s.factor(G.det())),difference_norm_over_C_squared=str(delta),
                residual_over_Js_four=str(physical),small_lambda_residual_over_Js_four=str(small),
                scope='Symbolic coefficients, cancellations, and projection. The accompanying proof supplies growing logarithmic windows and remainder bounds; this script alone is not that bound.',
                environment=dict(python=platform.python_version(),sympy=s.__version__),
                input_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print('passed: divergent coefficients',averaged[:2],'; A norm',unsigned,'+',averaged[2],'q; residual small lambda',small)


if __name__=='__main__':main()
