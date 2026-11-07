"""Central rank filters: independent finite formula checks and full-rank norms.

Never overwrites a prior result. Numerical refinement is not interval certification.
"""
from pathlib import Path
from fractions import Fraction
from math import comb
import hashlib,json,platform
import mpmath as mp
import sympy as sy
from cubic_linear_limit import distribution,scaled_coefficient
from cubic_linear_channels import coeff
from cubic_covariance_bridge import ratio

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-central-raw'

def calculate(lam,p,x,dps):
    with mp.workdps(dps):
        lam=mp.mpf(lam);N=int(p*p/lam);m=N//2+x*p
        w=distribution(N,p,m);q=mp.exp(-4*lam)
        t=[scaled_coefficient(N,p,m,p-u) for u in range(p+1)]
        mean=lambda vals:mp.fsum(a*b for a,b in zip(w,vals))
        center=lambda vals:[v-mean(vals) for v in vals]
        f=center(t);variance=mean([v*v for v in f])
        unit=[v/mp.sqrt(variance) for v in f]
        cv=lambda z:mp.mpf(z.numerator)/z.denominator
        r1=[cv(ratio(N,p,m,u,1)) for u in range(p+1)]
        r2=[cv(ratio(N,p,m,u,2)) for u in range(p+1)]
        proxy=[mp.sqrt(lam)*(r-1) for r in r1]
        ff=[lam*b+(1-2*lam)*a+lam for a,b in zip(r1,r2)]
        kap=[];sj=[]
        for u in range(p+1):
            a=p-u
            k=mp.fprod((1-mp.mpf(a)/(m-j))*(1-mp.mpf(a)/(N-m+p-j)) for j in range(p))
            J=mp.mpf(p)*(p*(N-m+1)-a*(N-a+1))/((m-a)*(N-m-a+1))
            kap.append(k);sj.append(k*J)
        bz=center(kap);bs=center(sj);mu=mean(range(p+1))
        # All coefficient errors are dimensionless, divided by C_m^3 in norm squared.
        slope_error=mean([(f[u]-8*lam*q*(u-mu))**2 for u in range(p+1)])
        measurements={
            'variance_over_C3':variance,'variance_limit_over_C3':64*lam**3*q*q,
            'variance_relative_limit_deviation':variance/(64*lam**3*q*q)-1,
            'F_slope_L2_error_squared_over_C3':slope_error,
            'normalized_F_proxy_L2_error_squared':mean([(a-b)**2 for a,b in zip(unit,proxy)]),
            'HF_covariance_L1_error':mean([abs(a-b) for a,b in zip(unit,proxy)]),
            'FF_covariance_L1_error':mean([abs(a*a-b) for a,b in zip(unit,ff)]),
            'pBZ_minus_Fhalf_norm_squared_over_C3':mean([(p*a-b/2)**2 for a,b in zip(bz,f)]),
            'pBS_minus_scaledF_norm_squared_over_C3':mean([(p*a-(1-2*lam)*b/2)**2 for a,b in zip(bs,f)]),
            'rank_deficit_mean':mu,'rank_probability_tail_ge_12':mp.fsum(w[12:]),
            'rank_weighted_slope_tail_ge_12':mp.fsum(w[u]*(f[u]-8*lam*q*(u-mu))**2 for u in range(12,p+1)),
        }
        for rr in [r1,r2]:assert abs(mean(rr)-1)<mp.mpf('1e-65')
        record=dict(N=N,p=p,m=m,x=x,lambda_value=str(lam),precision=dps,
                    **{k:mp.nstr(v,65) for k,v in measurements.items()})
        return record

def main():
    OUT.mkdir(exist_ok=True);target=OUT/'central-filter-checks.json'
    if target.exists():raise FileExistsError(target)
    prior=json.loads((OUT/'symbolic-descendants.json').read_text())
    for name,sha in prior['input_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha,name
    lam,x=sy.symbols('lam x')
    for row in prior['rows']:
        assert sy.factor(sy.sympify(row['central_B_shift'])-2*lam*x*sy.sympify(row['B_leading']))==0
    exact=[]
    with mp.workdps(100):
        for N in [8,10,12]:
            for p in [1,2,3]:
                for m in range(p+1,N-p):
                    C=comb(m,p)*comb(N-m+p,p)
                    for a in range(p+1):
                        expected=coeff(N,p,m,a)/C
                        actual=scaled_coefficient(N,p,m,a)
                        error=abs(actual-mp.mpf(expected.numerator)/expected.denominator)
                        assert error<mp.mpf('1e-90')
                        exact.append(dict(N=N,p=p,m=m,a=a,exact=str(expected)))
    records=[];changes=[]
    for lam in ['0.1','0.25','0.5']:
        for p in [16,32,64,128]:
            for x in [-2,0,2]:
                lo=calculate(lam,p,x,80);hi=calculate(lam,p,x,120)
                with mp.workdps(130):
                    errs=[abs(mp.mpf(lo[k])-mp.mpf(hi[k])) for k in hi if k not in ['N','p','m','x','lambda_value','precision']]
                    changes.append(max(errs));assert max(errs)<mp.mpf('1e-60')
                records.append(hi)
            print('central full-rank',lam,p,'passed',flush=True)
    names=['code/check_memory_central_raw.py','code/cubic_linear_limit.py','code/cubic_linear_channels.py',
           'code/cubic_covariance_bridge.py','results/susskind-memory-central-raw/symbolic-descendants.json']
    data=dict(status='passed',records=records,finite_exact_coefficient_reference=exact,
              central_shift_identities=5,precisions=[80,120],max_saved_precision_absolute_change=mp.nstr(max(changes),30),
              scope='All p+1 ranks; finite deviations retained. Precision comparisons use serialized 65-digit values; no interval or limit certificate.',
              environment=dict(python=platform.python_version(),mpmath=mp.__version__,sympy=sy.__version__),
              input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print('saved',len(exact),'finite reference coefficients;',len(records),'central sums',flush=True)

if __name__=='__main__':main()
