"""Full rank-channel sums, stable without enormous binomial coefficients."""
from pathlib import Path
import json,hashlib
import mpmath as mp
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-cubic-linear'
def distribution(N,p,m):
    vals=[mp.mpf(1)]
    for s in range(p):
        nu=p-s;a=mp.mpf(nu)/(N-nu+1);b=mp.mpf(nu-1)/(N-nu+2)
        ratio=a*a*(1-b*b)/(1-a*a)*mp.mpf(m-nu+1)/(p-nu+1)*mp.mpf(N-p-nu+1)/(N-m-nu+1)
        vals.append(vals[-1]*ratio)
    total=mp.fsum(vals)
    return [v/total for v in vals]
def scaled_coefficient(N,p,m,nu):
    N,p,m,nu=map(mp.mpf,(N,p,m,nu))
    cp=(m-p)/m*(N-m+p+1)/(N-m+1)
    jc=(m-p)/m
    base=-(m+2*p)*cp+(m-3*p)+3*p*jc
    kap=mp.fprod((1-nu/(m-j))*(1-nu/(N-m+p-j)) for j in range(int(p)))
    A=m-nu;B=N-m-nu+1
    al=nu*(N-nu-p+1)/(N-2*nu+1)
    ratio=(1-p/A)*(1+p/B)
    E=(m-p)*ratio-m+3*(al*(1+p/B)+(p-al)*(1-p/A))
    return base+kap*E
def run(lam,p,dps):
    mp.mp.dps=dps;lam=mp.mpf(lam);N=int(mp.nint(p*p/lam));assert N%2==0 and mp.mpf(p*p)/N==lam
    m=N//2;prob=distribution(N,p,m)
    coeff=[scaled_coefficient(N,p,m,p-s) for s in range(p+1)]
    avg=mp.fsum(w*c for w,c in zip(prob,coeff))
    var=mp.fsum(w*(c-avg)**2 for w,c in zip(prob,coeff))
    mean_s=mp.fsum(s*w for s,w in enumerate(prob))
    target=64*lam**3*mp.exp(-8*lam)
    return dict(N=N,p=p,m=m,lambda_value=str(lam),precision=dps,
      mean_s=mp.nstr(mean_s,60),variance_coefficient=mp.nstr(var,60),limit_variance=mp.nstr(target,60),
      relative_limit_deviation=mp.nstr(var/target-1,40),H_projection_coefficient_over_C=mp.nstr(avg,60),
      s1_slope=mp.nstr(coeff[1]-coeff[0],60),limit_slope=mp.nstr(8*lam*mp.exp(-4*lam),60),
      rank_tail_s_ge_12=mp.nstr(mp.fsum(prob[12:]),40),
      probabilities=[mp.nstr(x,50) for x in prob],coefficients_over_C=[mp.nstr(x,60) for x in coeff])
def main():
    target=OUT/'double-scale.json'
    if target.exists():raise FileExistsError(target)
    records=[];maxdiff=mp.mpf(0)
    for lam in ['0.1','0.25','0.5']:
        for p in [16,32,64,128]:
            low=run(lam,p,80);high=run(lam,p,110)
            diff=abs(mp.mpf(low['variance_coefficient'])/mp.mpf(high['variance_coefficient'])-1)
            maxdiff=max(maxdiff,diff);assert diff<mp.mpf('1e-55')
            records.append(high)
            print(lam,p,'mean',high['mean_s'][:12],'variance ratio-1',high['relative_limit_deviation'][:14],flush=True)
    target.write_text(json.dumps(dict(status='passed',scope='all p+1 rank channels, no rank cutoff; full high-precision binomial-ratio sums at half filling; finite-p differences retained',
      records=records,precisions=[80,110],max_relative_precision_change=mp.nstr(maxdiff,30),
      environment=dict(mpmath=mp.__version__),input_sha256={'code/cubic_linear_limit.py':hashlib.sha256((ROOT/'code/cubic_linear_limit.py').read_bytes()).hexdigest()}),indent=2)+'\n')
if __name__=='__main__':main()
