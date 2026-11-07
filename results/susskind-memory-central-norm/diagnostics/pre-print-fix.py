"""Central-window symbolic cancellation for the general-charge A norm."""
from pathlib import Path
import hashlib,json,platform
import sympy as s

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-central-norm'


def main():
    OUT.mkdir(exist_ok=True)
    target=OUT/'symbolic-limit.json'
    if target.exists():raise FileExistsError(target)
    j,k,u,v,lam,e,x=s.symbols('j k u v lam e x');p=1/e;N=p*p/lam;m=N/2+x*p;n=N-m+1
    a=p-u;b=p-v;cas=lambda a:a*(N+1-a)
    g=(m-p)**2/m;h=2*p-m;eta=p*p/m
    P=g*g*cas(a)*cas(b)/((m*n-cas(a))*(m*n-cas(b)))
    cf=(m-1)*n*(N+1)/N
    Ps=s.series(s.cancel(P*e**2),e,0,3).removeO().expand()
    Ds=s.series(s.cancel(((g*g-P)*(n-m)/n+2*h*g)*j/(m-j)*e**2),e,0,3).removeO().expand()
    unsigned=s.factor(s.series(s.cancel(eta**2+P*(m*n/cf-1)),e,0,1).removeO())
    # Verify the per-step logarithmic expansion from the exact term ratio.
    factor=s.cancel(((-a+k)*(a-N-1+k)*(-b+k)*(b-N-1+k))/((-m+k)*(m-N+k)*(-N-1+k)*(-4*lam)))
    step=s.series(factor,e,0,3).removeO().expand()
    assert step.coeff(e,0)==1
    l1=-u-v-2*k-2*lam
    l2=-((u+k)**2+(v+k)**2)/2+lam*(u+v+1+3*k)-lam**2+4*lam**2*x**2
    assert s.factor(step.coeff(e,1)-l1)==0
    assert s.factor(step.coeff(e,2)-l1*l1/2-l2)==0
    L1=-j*(u+v+2*lam)-j*(j-1)
    L2=-j*(u*u+v*v)/2-(u+v)*j*(j-1)/2-j*(j-1)*(2*j-1)/6+lam*j*(u+v+1)+s.Rational(3,2)*lam*j*(j-1)-lam*lam*j+4*lam**2*x**2*j
    assert s.factor(s.summation(l1,(k,0,j-1))-L1)==0
    assert s.factor(s.summation(l2,(k,0,j-1))-L2)==0
    expansion=s.Poly(s.expand((-Ps+Ds+s.series(s.cancel(eta**2*e**2),e,0,3).removeO())*(1+L1*e+(L2+L1**2/2)*e**2)),e)
    z=s.symbols('z');touchard=[s.Integer(1)]
    for degree in range(10):touchard.append(s.expand(z*(touchard[-1]+s.diff(touchard[-1],z))))
    def expectation(poly):
        poly=s.Poly(s.expand(poly),j)
        return s.factor(sum(c*touchard[power[0]].subs(z,-4*lam) for power,c in poly.terms()))
    averaged=[expectation(expansion.coeff_monomial(e**power)) for power in range(3)]
    assert averaged==[0,0,4*lam**2]
    assert unsigned==4*lam*(lam+1)
    # Record only the newly established central norm, not the inherited
    # half-filled four-field Gram involving the still-open output field.
    nu=p-u
    rr=(nu/(N-nu+1))**2*(1-((nu-1)/(N-nu+2))**2)/(1-(nu/(N-nu+1))**2)
    rr*=(m-nu+1)/(p-nu+1)*(N-p-nu+1)/(N-m-nu+1)
    ratio_limit=s.factor(s.limit(s.cancel(rr),e,0))
    assert s.factor(ratio_limit-lam/(u+1))==0
    names=['code/memory_central_norm_asymptotic.py','code/memory_norm_asymptotic.py']
    result=dict(status='passed',unsigned_limit=str(unsigned),log_terms=[str(L1),str(L2)],
                P_scaled_series=[str(s.factor(Ps.coeff(e,power))) for power in range(3)],
                summand_coefficients=[str(s.factor(expansion.coeff_monomial(e**power))) for power in range(3)],
                averaged_coefficients=list(map(str,averaged)),
                norm_over_C_squared='4*lambda + 4*lambda^2*(1+q)',
                input_deficit_successive_ratio_limit=str(ratio_limit),
                scope='Central-window general-charge coherent A norm; exact polynomial coefficients and cancellations. The general-charge cross formula and uniform logarithmic windows are supplied by accompanying derivation and independent CAR checks. This is not the complete central memory Gram.',
                environment=dict(python=platform.python_version(),sympy=s.__version__),
                input_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print('central passed: divergent coefficients',averaged[:2],'; A norm',unsigned,'+',averaged[2],'q; residual small lambda',small)


if __name__=='__main__':main()
