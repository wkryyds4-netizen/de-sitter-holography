"""Full-rank total-variation errors for the shifted-body covariance bridge."""
from pathlib import Path
from fractions import Fraction as F
from math import prod
import json,hashlib,platform
import mpmath as mp
from cubic_linear_dynamics import weights
from cubic_linear_limit import distribution,scaled_coefficient
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-linear-kernel'
def ratio(N,p,m,s,j):
    if s<j:return F(0)
    return F(prod((s-k)*(N-2*p+s+1+k)*(N-m+p-j+1+k) for k in range(j)),
             prod((m-p+1+k)*(p-k)**2 for k in range(j)))
def calc(lam,p,dps):
    mp.mp.dps=dps;lam=mp.mpf(lam);N=int(mp.nint(p*p/lam));m=N//2
    w=distribution(N,p,m)
    c=[scaled_coefficient(N,p,m,p-s) for s in range(p+1)]
    avg=mp.fsum(x*y for x,y in zip(w,c));delta=[x-avg for x in c]
    var=mp.fsum(x*y*y for x,y in zip(w,delta));f=[x/mp.sqrt(var) for x in delta]
    cv=lambda z:mp.mpf(z.numerator)/z.denominator
    r1=[cv(ratio(N,p,m,s,1)) for s in range(p+1)]
    r2=[cv(ratio(N,p,m,s,2)) for s in range(p+1)]
    fh=[mp.sqrt(lam)*(x-1) for x in r1]
    ff=[lam*y+(1-2*lam)*x+lam for x,y in zip(r1,r2)]
    eh=mp.fsum(w[s]*abs(f[s]-fh[s]) for s in range(p+1))
    ef=mp.fsum(w[s]*abs(f[s]**2-ff[s]) for s in range(p+1))
    assert abs(mp.fsum(w[s]*r1[s] for s in range(p+1))-1)<mp.mpf('1e-60')
    assert abs(mp.fsum(w[s]*r2[s] for s in range(p+1))-1)<mp.mpf('1e-60')
    return dict(N=N,p=p,m=m,lambda_value=str(lam),precision=dps,
      HF_covariance_L1_error=mp.nstr(eh,65),FF_covariance_L1_error=mp.nstr(ef,65),
      exact_HF_total_variation=mp.nstr(mp.fsum(w[s]*abs(f[s]) for s in range(p+1)),40),
      approximate_FF_total_variation=mp.nstr(mp.fsum(w[s]*abs(ff[s]) for s in range(p+1)),40),
      shifted_rank_probability_sums=[mp.nstr(mp.fsum(w[s]*rr[s] for s in range(p+1)),30) for rr in [r1,r2]],
      coefficient_total_variation_bound=mp.nstr(max(1,2*mp.sqrt(lam),2*lam+abs(1-2*lam)),30))
def main():
    OUT.mkdir(parents=True,exist_ok=True);path=OUT/'covariance-bridge.json'
    if path.exists():raise FileExistsError(path)
    checks=0
    for N in range(4,17):
        for m in range(2,N-1):
            for p in range(2,min(m,N-m,4)+1):
                w=weights(N,m,p)
                for j in range(3):
                    z=weights(N,m,p-j)
                    for a in range(p+1):
                        assert w[a]*ratio(N,p,m,p-a,j)==(z[a] if a<=p-j else 0)
                        checks+=1
    rows=[];changes=[]
    for lam in ['0.1','0.25','0.5','1.0']:
        for p in [16,32,64,128,256]:
            a=calc(lam,p,80);b=calc(lam,p,120)
            for name in ['HF_covariance_L1_error','FF_covariance_L1_error']:
                change=abs(mp.mpf(a[name])/mp.mpf(b[name])-1);assert change<mp.mpf('1e-60');changes.append(change)
            rows.append(b)
            print(lam,p,'L1',b['HF_covariance_L1_error'][:12],b['FF_covariance_L1_error'][:12],flush=True)
    names=['code/cubic_covariance_bridge.py','code/cubic_linear_dynamics.py','code/cubic_linear_channels.py','code/cubic_linear_limit.py']
    path.write_text(json.dumps(dict(status='passed',exact_ratio_checks=checks,records=rows,precisions=[80,120],
      max_relative_precision_change=mp.nstr(max(changes),30),
      scope='All p+1 ranks; L1 of signed scalar covariance-channel weights; 65-digit saved precision comparison, not intervals',
      environment=dict(python=platform.python_version(),mpmath=mp.__version__),
      input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')
if __name__=='__main__':main()
