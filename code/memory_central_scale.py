"""Central charge normalization drift and exact raw-filter asymptotics."""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import hashlib,json,platform
import sympy as s

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-central-window'


def main():
    target=OUT/'central-scale.json'
    if target.exists():raise FileExistsError(target)
    e,lam,x,y,t,v=s.symbols('e lam x y t v',real=True)
    # e=1/p, i=y*p, charge shift=x*p, N=p^2/lam.
    local=s.log(1+2*lam*x*e/(1-2*lam*y*e))+s.log(1-2*lam*x*e/(1+2*lam*(1-y)*e))
    local_series=s.series(local,e,0,3).removeO().expand()
    drift=s.factor(local_series.coeff(e,2))
    assert s.factor(drift-4*lam**2*x*(1-x))==0
    # rank j=t*p-v; paired factors use m and N-m+p.
    ranklog=s.log(1-(t/e-v)/(1/(2*lam*e**2)+(x-y)/e))
    ranklog+=s.log(1-(t/e-v)/(1/(2*lam*e**2)+(1-x-y)/e))
    rankseries=s.series(ranklog,e,0,3).removeO().expand()
    ranklead=s.factor(rankseries.coeff(e,1))
    ranknext=s.factor(s.integrate(rankseries.coeff(e,2),(y,0,1)))
    assert ranklead==-4*lam*t
    assert s.factor(ranknext-(4*lam*v-4*lam**2*t**2))==0
    # The e coefficient is independent of y, so the sum-to-integral endpoint
    # correction does not change ranknext.
    old=json.loads((ROOT/'results/susskind-memory-two-charge/symbolic-filter.json').read_text())
    m,n,p,a,b,r=s.symbols('m n p a b r')
    u,v,S,q=s.symbols('u v S q')
    sub={m:1/(2*lam*e**2)+x/e,n:1/(2*lam*e**2)-x/e+1,
         p:1/e,a:1/e-u,b:1/e-v,r:2/e-u-v-S}
    expansions=[]
    for val in old['H_squared_coefficients']:
        ex=s.series(s.cancel(s.sympify(val).subs(sub)),e,0,1).removeO().expand()
        expansions.append(s.factor(ex))
    kexp=[1,q*(1+e*(4*lam*u-4*lam**2)),q*(1+e*(4*lam*v-4*lam**2)),
          q*q*(1+e*(4*lam*(u+v+S)-16*lam**2))]
    combined=s.expand(sum(z*k for z,k in zip(expansions,kexp)))
    lead=s.factor(combined.coeff(e,-1));finite=s.factor(combined.coeff(e,0))
    assert s.factor(lead-2*(1-q)**2)==0
    half=s.sympify(old['gH_over_C_finite'],locals={'lam':lam,'u':u,'v':v,'S':S,'q':q})
    difference=s.factor(finite-half)
    pref=json.loads((ROOT/'results/susskind-memory-filter-limit/filter-certificate.json').read_text())
    aproportions=[]
    alpha,beta,rho=s.symbols('alpha beta rho')
    prop=dict(sub);prop.update({a:alpha/e,b:beta/e,r:rho/e})
    for val in pref['exact_prefactors']:
        ex=s.series(s.cancel(s.sympify(val).subs(prop)),e,0,2).removeO().expand()
        aproportions.append(s.factor(ex))
    records=[]
    for pp in [16,32,64,128,256]:
        nn=4*pp*pp;half_m=nn//2
        base=comb(half_m,pp)*comb(nn-half_m+pp,pp)
        for xx in [F(-2),F(-1),F(0),F(1,2),F(1),F(2)]:
            mm=half_m+int(xx*pp)
            rat=F(comb(mm,pp)*comb(nn-mm+pp,pp),base)
            reflected=nn+pp-mm
            assert comb(reflected,pp)*comb(nn-reflected+pp,pp)==rat*base
            value=pp*(rat-1);expected=F(1,4)*xx*(1-xx)
            records.append(dict(N=nn,p=pp,x=str(xx),ratio=str(rat),scaled_difference=float(value),
                                limit=str(expected),absolute_limit_deviation=float(value-expected)))
    names=['code/memory_central_scale.py','results/susskind-memory-two-charge/symbolic-filter.json',
           'results/susskind-memory-filter-limit/filter-certificate.json']
    data=dict(status='passed',local_log_ratio_coefficient=str(drift),
        rank_log_leading=str(ranklead),rank_log_first_correction=str(ranknext),
        H_prefactor_expansions=list(map(str,expansions)),gH_leading=str(lead),gH_finite=str(finite),
        gH_central_minus_half=str(difference),A_prefactor_expansions=list(map(str,aproportions)),
        exact_variance_checks=records,
        scope='Fixed compact x=(m-N/2)/p and positive lambda. Rational coefficient and binomial-product asymptotics; not a complete operator Gram convergence theorem.',
        environment=dict(python=platform.python_version(),sympy=s.__version__),
        input_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print('central variance drift',drift,'H finite shift',difference,'A prefactors',aproportions,flush=True)


if __name__=='__main__':main()
