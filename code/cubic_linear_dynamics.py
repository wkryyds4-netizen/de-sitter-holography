"""Rank-resolved covariance and same-disorder quadratic dynamics.

Exact Fraction formulas and a full-rank mpmath contraction; no rank cutoff.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import hashlib,json,platform
import mpmath as mp
from cubic_linear_channels import coeff
from cubic_linear_limit import distribution,scaled_coefficient

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-cubic-dynamics'
def b(n,k):return comb(n,k) if 0<=k<=n else 0
def weights(N,m,r):
    C=b(m,r)*b(N-m+r,r)
    return [F((b(N,a)**2-b(N,a-1)**2)*b(m-a,m-r)*b(N-r-a,m-r),b(N,m)*C) for a in range(r+1)]
def crossing(N,m,a,n):
    term=F(1);value=term
    for k in range(min(a,n)):
        term*=F((-a+k)*(a-N-1+k)*(-n+k)*(n-N-1+k),(-m+k)*(m-N+k)*(-N-1+k)*(k+1))
        value+=term
    return value
def finite(N,p,m):
    w=weights(N,m,p);C=b(m,p)*b(N-m+p,p)
    t=[coeff(N,p,m,a) for a in range(p+1)]
    avg=sum(x*y for x,y in zip(w,t));f=[x-avg for x in t]
    var=sum(w[a]*f[a]**2 for a in range(p+1))
    X=sum(w[a]*w[n]*f[a]*f[n]*crossing(N,m,a,n) for a in range(p+1) for n in range(p+1))
    U=sum(w[a]*w[n]*f[a]**2*crossing(N,m,a,n) for a in range(p+1) for n in range(p+1))
    direct=C*C*sum(w[a]*w[n]*(f[a]-f[n])**2*(1-crossing(N,m,a,n)) for a in range(p+1) for n in range(p+1))
    assert direct==2*C*C*(var+X-U)
    return dict(C=F(C),norm=C*var,M2=direct,B20=C*C*(var+X),B11=C*C*U,cross=C*C*X)

def ds(lam,p,dps):
    mp.mp.dps=dps;lam=mp.mpf(lam);N=int(mp.nint(p*p/lam));m=N//2
    assert N%2==0 and mp.mpf(p*p)/N==lam
    w=list(reversed(distribution(N,p,m)))
    c=[scaled_coefficient(N,p,m,a) for a in range(p+1)]
    avg=mp.fsum(x*y for x,y in zip(w,c));f=[x-avg for x in c]
    var=mp.fsum(w[a]*f[a]**2 for a in range(p+1))
    z=[mp.mpf(1)]*(p+1);h=mp.mpf(1);xt=[];ut=[]
    for k in range(p+1):
        a0=mp.fsum(w[a]*z[a] for a in range(p+1))
        a1=mp.fsum(w[a]*f[a]*z[a] for a in range(p+1))
        a2=mp.fsum(w[a]*f[a]**2*z[a] for a in range(p+1))
        xt.append(h*a1*a1);ut.append(h*a0*a2)
        if k<p:
            for a in range(p+1):z[a]*=(-a+k)*(a-N-1+k)
            h/=(-m+k)*(m-N+k)*(-N-1+k)*(k+1)
    X=mp.fsum(xt);U=mp.fsum(ut);m2=2*(var+X-U)
    q=mp.exp(-4*lam);target=2*(1-q)/lam
    mu2=m2/var/lam # Js=1, C -> 1/lambda; finite C fixed to 1/lambda.
    g=mp.fsum(ut)/var
    return dict(N=N,p=p,m=m,lambda_value=str(lam),precision=dps,
      convention='sector variance C=Js^2/lambda, Js=1; full p+1 rank and terminating series sums',
      variance_coefficient=mp.nstr(var,65),M2_over_C4=mp.nstr(m2,65),
      normalized_mu2=mp.nstr(mu2,65),limit_mu2=mp.nstr(target,65),
      relative_limit_deviation=mp.nstr(mu2/target-1,45),
      B20_over_Cnorm=mp.nstr(1+X/var,65),B11_over_Cnorm=mp.nstr(g,65),
      same_disorder_cross_over_Cnorm=mp.nstr(X/var,65),
      contraction_cancellation_factor=mp.nstr((var+mp.fsum(map(abs,xt))+mp.fsum(map(abs,ut)))/abs(var+X-U),30))
def main():
    OUT.mkdir(parents=True,exist_ok=True);path=OUT/'double-scale.json'
    if path.exists():raise FileExistsError(path)
    rows=[];changes=[]
    for lam in ['0.1','0.25','0.5']:
        for p in [16,32,64,128]:
            lo=ds(lam,p,80);hi=ds(lam,p,120)
            change=abs(mp.mpf(lo['normalized_mu2'])/mp.mpf(hi['normalized_mu2'])-1)
            assert change<mp.mpf('1e-60');changes.append(change);rows.append(hi)
            print(lam,p,'mu2',hi['normalized_mu2'][:15],'relative',hi['relative_limit_deviation'][:15],flush=True)
    names=['code/cubic_linear_dynamics.py','code/cubic_linear_channels.py','code/cubic_linear_limit.py']
    path.write_text(json.dumps(dict(status='passed',records=rows,precisions=[80,120],
      precision_scope='65-digit stored normalized moments; no interval rounding certificate',
      max_relative_precision_change=mp.nstr(max(changes),30),
      environment=dict(python=platform.python_version(),mpmath=mp.__version__),
      input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')
if __name__=='__main__':main()
