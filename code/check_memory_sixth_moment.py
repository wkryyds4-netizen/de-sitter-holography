"""Independent finite-distribution and large-parameter checks through order six."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,platform
import mpmath as mp
from memory_sixth_moment import finite_deficit_moments
from rank_product_moments import product_distribution

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-filter-limit'


def touchard6(t):return t**6+15*t**5+65*t**4+90*t**3+31*t*t+t


def main():
    target=OUT/'sixth-moment-checks.json'
    if target.exists():raise FileExistsError(target)
    small=[]
    for N,m,a,b in [(14,6,2,3),(14,7,3,3),(16,7,3,4),(16,8,4,4)]:
        pi=product_distribution(N,m,a,b);mom=finite_deficit_moments(N,m,a,b)
        x=lambda r:r*(N+1-r)
        actual=[sum(v*F(x(a+b)-x(r),N)**k for r,v in enumerate(pi)) for k in range(7)]
        assert actual==mom,(N,m,a,b,actual,mom)
        sixthS=sum(v*(a+b-r)**6 for r,v in enumerate(pi))
        assert sixthS<=mom[6]*F(N,N+1-2*(a+b))**6
        small.append(dict(N=N,m=m,a=a,b=b,moments=list(map(str,mom)),sixth_actual_S=str(sixthS)))
    rows=[]
    for lam in [F(1,10),F(1,4),F(1,2)]:
        for p in [16,32,64,128]:
            N=int(p*p/lam)
            for shift in [0,-1]:
                m=N//2+shift
                for u,v in [(0,0),(1,2)]:
                    mom=finite_deficit_moments(N,m,p-u,p-v)
                    upper=mom[6]*F(N,N+1-4*p)**6
                    runs=[]
                    for precision in [80,120]:
                        mp.mp.dps=precision
                        cv=lambda z:mp.mpf(z.numerator)/z.denominator
                        runs.append([cv(mom[6]),cv(upper),cv(touchard6(2*lam))])
                    assert all(abs(a-b)<mp.mpf('1e-65') for a,b in zip(*runs))
                    vals=runs[-1]
                    rows.append(dict(N=N,p=p,charge_shift=shift,input_deficits=[u,v],lambda_value=str(lam),
                                     sixth_delta_casimir_over_N=mp.nstr(vals[0],65),
                                     rigorous_S_sixth_upper_bound_exact=str(upper),
                                     S_sixth_upper_bound_decimal=mp.nstr(vals[1],65),
                                     limit=mp.nstr(vals[2],65),
                                     proxy_absolute_limit_deviation=mp.nstr(vals[0]-vals[2],40)))
    names=['code/check_memory_sixth_moment.py','code/memory_sixth_moment.py',
           'code/rank_product_moments.py','code/cubic_linear_dynamics.py']
    result=dict(status='passed',finite_full_distribution_checks=small,
                finite_moment_equalities=28,records=rows,precisions=[80,120],
                scope='Exact Fraction finite moments before decimal evaluation; exact rational S-sixth upper bounds. Large-system proxy is X/N, not S itself. No full Poisson-law claim or decimal interval certificate.',
                environment=dict(python=platform.python_version(),mpmath=mp.__version__),
                input_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print('passed 28 finite moment identities and',len(rows),'large-parameter cases; last deviation',rows[-1]['proxy_absolute_limit_deviation'])


if __name__=='__main__':main()
