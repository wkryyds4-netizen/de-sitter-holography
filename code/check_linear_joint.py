"""Finite-rank quartic law of h+z*y, retaining odd mixed contractions."""
from pathlib import Path
import hashlib,json
import mpmath as mp
from cubic_linear_limit import distribution,scaled_coefficient
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-linear-kernel'
def compute(lam,p,dps):
    mp.mp.dps=dps;lam=mp.mpf(lam);N=int(p*p/lam);m=N//2
    w=list(reversed(distribution(N,p,m)))
    c=[scaled_coefficient(N,p,m,a) for a in range(p+1)]
    mean=mp.fsum(w[a]*c[a] for a in range(p+1))
    var=mp.fsum(w[a]*(c[a]-mean)**2 for a in range(p+1))
    u=[(x-mean)/mp.sqrt(var) for x in c]
    z=[mp.mpf(1)]*(p+1);h=mp.mpf(1);terms=[[] for _ in range(5)]
    for k in range(p+1):
        A=[mp.fsum(w[a]*u[a]**j*z[a] for a in range(p+1)) for j in range(3)]
        vals=[A[0]**2,4*A[0]*A[1],2*A[0]*A[2]+4*A[1]**2,4*A[1]*A[2],A[2]**2]
        for j in range(5):terms[j].append(h*vals[j])
        if k<p:
            for a in range(p+1):z[a]*=(-a+k)*(a-N-1+k)
            h/=(-m+k)*(m-N+k)*(-N-1+k)*(k+1)
    val=[mp.fsum(terms[j])+[2,0,4,0,2][j] for j in range(5)]
    q=mp.exp(-4*lam);target=[2+q,0,2*(2+q),0,2+q]
    return dict(N=N,p=p,lambda_value=str(lam),precision=dps,
      polynomial_coefficients=[mp.nstr(x,65) for x in val],
      limit_coefficients=[mp.nstr(x,65) for x in target],
      max_absolute_limit_deviation=mp.nstr(max(abs(x-y) for x,y in zip(val,target)),30))
def main():
    path=OUT/'joint-quartic.json'
    if path.exists():raise FileExistsError(path)
    rows=[];errors=[]
    for lam in ['0.1','0.5','1.0']:
        for p in [16,64,256]:
            a=compute(lam,p,80);b=compute(lam,p,120)
            error=max(abs(mp.mpf(x)-mp.mpf(y)) for x,y in zip(a['polynomial_coefficients'],b['polynomial_coefficients']))
            assert error<mp.mpf('1e-60');errors.append(error);rows.append(b)
            print(lam,p,'quartic deviation',b['max_absolute_limit_deviation'][:14],flush=True)
    names=['code/check_linear_joint.py','code/cubic_linear_limit.py']
    path.write_text(json.dumps(dict(status='passed',records=rows,precisions=[80,120],
      scope='Exact rank-covariance formula for Etr(h+z*y)^4, all ranks and all terminating Racah terms; finite deviations retained',
      max_absolute_precision_change=mp.nstr(max(errors),30),
      input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')
if __name__=='__main__':main()
