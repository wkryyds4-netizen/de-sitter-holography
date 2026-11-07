"""Full-rank L2 estimates for the S/Z linear counterterms."""
from pathlib import Path
import hashlib,json
import mpmath as mp
from cubic_linear_limit import distribution,scaled_coefficient
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-cubic-cross'
def calc(lam,p,dps):
    mp.mp.dps=dps;lam=mp.mpf(lam);N=int(p*p/lam);m=N//2
    w=list(reversed(distribution(N,p,m)));z=[];ss=[];tt=[]
    cp=mp.mpf(m-p)/m*mp.mpf(N-m+p+1)/(N-m+1)
    a0=2*p-m+m*cp;K0=a0+mp.mpf(p*p)/m
    for a in range(p+1):
        kap=mp.fprod((1-mp.mpf(a)/(m-j))*(1-mp.mpf(a)/(N-m+p-j)) for j in range(p))
        J=mp.mpf(p)*(p*(N-m+1)-a*(N-a+1))/((m-a)*(N-m-a+1))
        z.append(2+kap);ss.append(K0+kap*J);tt.append(scaled_coefficient(N,p,m,a))
    avg=lambda v:mp.fsum(x*y for x,y in zip(w,v))
    center=lambda v:[x-avg(v) for x in v]
    bz=center(z);bs=center(ss);f=center(tt);norm=avg([x*x for x in f])
    c=(avg(ss)-a0)/(avg(z)-1);clim=2*lam*(1-mp.exp(-4*lam))/(1+mp.exp(-4*lam))
    ez=avg([(p*x-y/2)**2 for x,y in zip(bz,f)])/norm
    es=avg([(p*x-(1-2*lam)*y/2)**2 for x,y in zip(bs,f)])/norm
    ex=avg([(p*(x-c*zv)-(1-2*lam-clim)*y/2)**2 for x,zv,y in zip(bs,bz,f)])/norm
    return dict(N=N,p=p,m=m,lambda_value=str(lam),precision=dps,
      Z_scaled_L2_error_squared_over_Fnorm=mp.nstr(ez,65),
      S_scaled_L2_error_squared_over_Fnorm=mp.nstr(es,65),
      Jordan_scaled_L2_error_squared_over_Fnorm=mp.nstr(ex,65),
      BZ_norm_squared_over_Fnorm=mp.nstr(avg([x*x for x in bz])/norm,65),
      BS_norm_squared_over_Fnorm=mp.nstr(avg([x*x for x in bs])/norm,65),
      finite_quadratic_projection=mp.nstr(c,65),limit_projection=mp.nstr(clim,65))
def main():
    path=OUT/'counterterm-limits.json'
    if path.exists():raise FileExistsError(path)
    rows=[];errs=[]
    fields=['Z_scaled_L2_error_squared_over_Fnorm','S_scaled_L2_error_squared_over_Fnorm','Jordan_scaled_L2_error_squared_over_Fnorm']
    for lam in ['0.1','0.25','0.5']:
        for p in [16,32,64,128,256]:
            a=calc(lam,p,80);b=calc(lam,p,120)
            e=max(abs(mp.mpf(a[k])-mp.mpf(b[k])) for k in fields);assert e<mp.mpf('1e-60');errs.append(e);rows.append(b)
            print(lam,p,'scaled errors squared',*[b[k][:12] for k in fields],flush=True)
    names=['code/run_cubic_cross_limits.py','code/cubic_linear_limit.py']
    path.write_text(json.dumps(dict(status='passed',records=rows,precisions=[80,120],max_absolute_precision_change=mp.nstr(max(errs),30),
      scope='Full rank sums, no cutoff; errors are squared L2 errors divided by F norm, not moment errors',
      input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')
if __name__=='__main__':main()
