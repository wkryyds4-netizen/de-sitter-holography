"""All-rank high-precision A norm and finite Casimir-proxy projection."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,platform,time
import mpmath as mp
from cubic_linear_limit import distribution
from memory_residual_norm import exact_norm

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-residual-norm'


def evaluate(N,p,precision):
    mp.mp.dps=precision;m=N//2;n=m+1
    lam=mp.mpf(p*p)/N;eta=mp.mpf(p*p)/m
    omega=list(reversed(distribution(N,p,m)))
    xa=[a*(N+1-a) for a in range(p+1)]
    fa=[mp.mpf(x)/(m*n-x) for x in xa]
    ma=mp.fsum(a*omega[a] for a in range(p+1));mx=mp.fsum(xa[a]*omega[a] for a in range(p+1))
    mf=mp.fsum(fa[a]*omega[a] for a in range(p+1))
    cf=mp.mpf(m*m)*(N+1)/N;cfm=mp.mpf((m-1)*n)*(N+1)/N
    g=mp.mpf((m-p)**2)/m;h=2*p-m
    z=[mp.mpf(1)]*(p+1);hh=mp.mpf(1);sums=[[] for _ in range(7)]
    for j in range(p+1):
        z0=mp.fsum(omega[a]*z[a] for a in range(p+1))
        z1=mp.fsum(omega[a]*a*z[a] for a in range(p+1))
        z2=mp.fsum(omega[a]*a*a*z[a] for a in range(p+1))
        zx=mp.fsum(omega[a]*xa[a]*z[a] for a in range(p+1))
        zf=mp.fsum(omega[a]*fa[a]*z[a] for a in range(p+1))
        zd=mp.fsum(omega[a]*a*a*(m+1-a)*(z[a-1]-z[a]) for a in range(1,p+1))
        hd=hh*mp.mpf(j)/(m-j)
        for array,val in zip(sums,[hh*z0*z0,hh*z1*z1,hh*zf*zf,hd*z0*z0,hd*z1*z1,hd*zf*zf,
                                   hh*(4*z2*z2-2*z1*zx+zd*z0)]):array.append(val)
        if j<p:
            for a in range(p+1):z[a]*=(-a+j)*(a-N-1+j)
            hh/=(-m+j)**2*(-N-1+j)*(j+1)
    chi,abchi,chiff,delta,abdelta,deltaff,signedXX=map(mp.fsum,sums)
    parts=[eta*eta*(1+chi),g*g*mf*mf*(m*n/cfm-1),-g*g*chiff,
           (g*g/n+2*h*g)*delta,-g*g/n*deltaff]
    avar=mp.fsum(parts)
    bvar=1+chi;bA=eta*(1+chi)+g*delta
    ux=-2*ma*ma+mx*mx/cf
    ud=mx*mx/(mp.mpf(m)*(m-1)*(N+1))
    xA=(eta*(ux-2*abchi)+g*(ud+m*delta-2*abdelta))/N
    bx=(ux-2*abchi)/N
    a21=-mp.mpf(2*(N-1))/(m*m*(N+1))
    a22=mp.mpf((N-1)*(N-2))/(m*m*(m-1)**2*N*(N+1))
    chi2=[1+a21*x+a22*x*(x-N) for x in xa]
    unsigned=[]
    for a in range(p+1):
        for b in range(p+1):
            er=xa[a]+xa[b]-mp.mpf(xa[a])*xa[b]/cf
            falling2=(chi2[a]*chi2[b]-1-a21*er)/a22
            top=(a+b)*(N+1-a-b)
            unsigned.append(omega[a]*omega[b]*(top*top-2*top*er+falling2+N*er))
    xvar=(mp.fsum(unsigned)+signedXX)/(N*N)
    diff=avar+xvar-2*xA
    det=bvar*avar-bA*bA
    assert det>0
    residual=xvar-(avar*bx*bx-2*bA*bx*xA+bvar*xA*xA)/det
    assert avar>0 and diff>0 and residual>0
    q=mp.exp(-4*lam);vs=2*lam*(1-q)+16*lam*lam*q/(1+q);dd=2*lam*(1+q)
    limits=[4*lam+4*lam*lam*(1+q),dd,vs*dd/(vs+dd)]
    vals=[avar,diff,residual]
    return dict(N=N,p=p,lambda_value=str(lam),precision=precision,
                values=[mp.nstr(v,65) for v in vals],limits=[mp.nstr(v,65) for v in limits],
                absolute_deviations=[mp.nstr(v-t,40) for v,t in zip(vals,limits)],
                gram_B_X_A=[[mp.nstr(v,65) for v in row] for row in [[bvar,bx,bA],[bx,xvar,xA],[bA,xA,avar]]],
                A_cancellation_factor=mp.nstr(mp.fsum(map(abs,parts))/abs(avar),30))


def main():
    target=OUT/'limit-checks.json'
    if target.exists():raise FileExistsError(target)
    start=time.time();old=ROOT/'results/susskind-memory-coherent-cross/car-cross.json'
    finite=json.loads(old.read_text());checks=[]
    for row in finite['rows']:
        N,p,m=row['N'],row['p'],row['m'];C=__import__('math').comb(m,p)*__import__('math').comb(N-m+p,p)
        out=evaluate(N,p,120)
        ids=[0,2,3]
        for i,a in enumerate(ids):
            for j,b in enumerate(ids):
                exact=F(row['exact_gram'][a][b])/(C*C)
                assert abs(mp.mpf(out['gram_B_X_A'][i][j])-mp.mpf(exact.numerator)/exact.denominator)<mp.mpf('1e-60')
        checks.append(dict(N=N,p=p,gram_entries=9))
    # A further exact all-pair formula comparison outside the saved CAR systems.
    for N,p in [(16,4),(24,4)]:
        out=evaluate(N,p,120);exact=exact_norm(N,p)
        assert abs(mp.mpf(out['values'][0])-mp.mpf(exact.numerator)/exact.denominator)<mp.mpf('1e-60')
    rows=[];changes=[]
    for lam in [F(1,10),F(1,4),F(1,2)]:
        for p in [16,32,64,128]:
            N=int(p*p/lam)
            low=evaluate(N,p,80);high=evaluate(N,p,120)
            change=max(abs(mp.mpf(v)-mp.mpf(w)) for v,w in zip(low['values'],high['values']))
            assert change<mp.mpf('1e-60');changes.append(change);rows.append(high)
            print('norm',lam,p,high['absolute_deviations'],flush=True)
    names=['code/check_memory_norm_limit.py','code/memory_residual_norm.py','code/cubic_linear_limit.py',
           'code/cubic_linear_dynamics.py','code/memory_coherent_cross.py',
           'results/susskind-memory-coherent-cross/car-cross.json']
    data=dict(status='passed',finite_gram_checks=checks,additional_exact_norms=2,
              records=rows,precisions=[80,120],max_absolute_precision_change=mp.nstr(max(changes),30),
              quantities=['||W2 A||^2/C^2','||W2 A-X_N/N||^2/C^2','X_N/N residual after B,A divided by C^2'],
              scope='Complete p+1 rank sums and terminating series. The last two finite quantities use the Casimir proxy X_N/N, not S_N. Analytic L2 convergence justifies their limiting targets. Decimal precision refinement is not an interval certificate.',
              elapsed_seconds=time.time()-start,
              environment=dict(python=platform.python_version(),mpmath=mp.__version__),
              input_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names})
    target.write_text(json.dumps(data,indent=2)+'\n')


if __name__=='__main__':main()
