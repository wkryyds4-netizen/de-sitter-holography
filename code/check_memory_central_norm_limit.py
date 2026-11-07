"""Full-rank central-window B,L,O Gram and weighted A norm sums."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,platform,time
import mpmath as mp
from cubic_linear_limit import distribution
from memory_central_norm import exact_gram

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-central-norm'


def evaluate(N,p,m,dps):
    mp.mp.dps=dps;n=N-m+1;lam=mp.mpf(p*p)/N
    omega=list(reversed(distribution(N,p,m)))
    u=[p-a for a in range(p+1)];mean=mp.fsum(w*v for w,v in zip(omega,u))
    var=mp.fsum(w*(v-mean)**2 for w,v in zip(omega,u))
    fa=[mp.mpf(a*(N+1-a))/(m*n-a*(N+1-a)) for a in range(p+1)]
    g=mp.mpf((m-p)**2)/m;h=2*p-m;eta=mp.mpf(p*p)/m
    cf=mp.mpf((m-1)*n)*(N+1)/N;alpha=mp.mpf(n-m)/n
    def moments(z):return [mp.fsum(omega[a]*z[a]*u[a]**k for a in range(p+1)) for k in range(3)]
    def bilinear(z,weighted):
        z0,z1,z2=z
        return z0*z0+4*z0*z1+2*z0*z2+2*z1*z1 if weighted else z0*z0
    un=moments([mp.mpf(1)]*(p+1));ff=moments(fa)
    series=[[] for _ in range(10)];zz=[mp.mpf(1)]*(p+1);hh=mp.mpf(1)
    for j in range(p+1):
        z=moments(zz);zf=moments([v*f for v,f in zip(zz,fa)])
        hd=hh*mp.mpf(j)/(m-j)
        values=[]
        for weighted in [False,True]:
            values.extend([hh*bilinear(z,weighted),hh*bilinear(zf,weighted),
                           hd*bilinear(z,weighted),hd*bilinear(zf,weighted)])
        centered=z[1]-mean*z[0];centered2=z[2]-2*mean*z[1]+mean*mean*z[0]
        values.extend([2*hh*z[0]*centered,hh*(2*z[0]*centered2+2*centered*centered)])
        for arr,v in zip(series,values):arr.append(v)
        if j<p:
            for a in range(p+1):zz[a]*=(-a+j)*(a-N-1+j)
            hh/=(-m+j)*(m-N+j)*(-N-1+j)*(j+1)
    sums=list(map(mp.fsum,series));norms=[];conditions=[]
    for weighted,offset in [(False,0),(True,4)]:
        chi,chiff,delta,deltaff=sums[offset:offset+4]
        parts=[eta*eta*(bilinear(un,weighted)+chi),
               g*g*bilinear(ff,weighted)*(m*n/cf-1),-g*g*chiff,
               (g*g*alpha+2*h*g)*delta,-g*g*alpha*deltaff]
        val=mp.fsum(parts);assert val>0;norms.append(val)
        conditions.append(mp.fsum(map(abs,parts))/abs(val))
    chi,chiff,delta,deltaff=sums[:4];BL=sums[8];LL=2*var+sums[9]
    # Compute L--O by inserting the centered input weight in delta.
    zz=[mp.mpf(1)]*(p+1);hh=mp.mpf(1);dL=[]
    for j in range(p+1):
        z=moments(zz)
        dL.append(2*hh*mp.mpf(j)/(m-j)*z[0]*(z[1]-mean*z[0]))
        if j<p:
            for a in range(p+1):zz[a]*=(-a+j)*(a-N-1+j)
            hh/=(-m+j)*(m-N+j)*(-N-1+j)*(j+1)
    BB=1+chi;BO=eta*BB+g*delta;LO=eta*BL+g*mp.fsum(dL)
    gram=[[BB,BL,BO],[BL,LL,LO],[BO,LO,norms[0]]]
    assert mp.det(mp.matrix(gram))>0
    q=mp.exp(-4*lam);aa=4*lam+4*lam*lam*(1+q)
    target=[[1+q,mp.mpf(0),2*lam*(1-q)],
            [mp.mpf(0),2*lam*(1+q),mp.mpf(0)],
            [2*lam*(1-q),mp.mpf(0),aa]]
    return dict(N=N,p=p,m=m,x=str(F(2*m-N,2*p)),lambda_value=str(lam),precision=dps,
        gram=[[mp.nstr(v,65) for v in row] for row in gram],
        limit_gram=[[mp.nstr(v,65) for v in row] for row in target],
        max_absolute_gram_deviation=mp.nstr(max(abs(v-w) for row,tar in zip(gram,target) for v,w in zip(row,tar)),40),
        weighted_A_norm=mp.nstr(norms[1],65),weighted_A_limit=mp.nstr(aa*((1+2*lam)**2+2*lam),65),
        A_relative_limit_error=mp.nstr(norms[0]/aa-1,40),
        weighted_A_relative_limit_error=mp.nstr(norms[1]/(aa*((1+2*lam)**2+2*lam))-1,40),
        cancellation_factors=list(map(lambda v:mp.nstr(v,30),conditions)),
        input_mean=mp.nstr(mean,50),input_variance=mp.nstr(var,50))


def main():
    target=OUT/'limit-checks.json'
    if target.exists():raise FileExistsError(target)
    start=time.time();checks=[]
    for N,p,m in [(5,2,3),(6,2,4),(16,3,7),(16,3,9)]:
        exact=exact_gram(N,p,m);data=evaluate(N,p,m,120)
        for i,ii in enumerate([0,2]):
            for j,jj in enumerate([0,2]):
                z=exact[i][j]
                assert abs(mp.mpf(data['gram'][ii][jj])-mp.mpf(z.numerator)/z.denominator)<mp.mpf('1e-60')
        checks.append(dict(N=N,p=p,m=m,exact_entries=4))
    rows=[];changes=[]
    for lam in [F(1,10),F(1,4),F(1,2)]:
        for p in [16,32,64,128]:
            N=int(p*p/lam)
            for x in [-2,0,2]:
                m=N//2+x*p;low=evaluate(N,p,m,80);high=evaluate(N,p,m,120)
                lv=[v for row in low['gram'] for v in row]+[low['weighted_A_norm']]
                hv=[v for row in high['gram'] for v in row]+[high['weighted_A_norm']]
                change=max(abs(mp.mpf(a)-mp.mpf(b)) for a,b in zip(lv,hv))
                assert change<mp.mpf('1e-55');changes.append(change);rows.append(high)
            print('all rank central Gram',str(lam),p,'passed',flush=True)
    names=['code/check_memory_central_norm_limit.py','code/memory_central_norm.py',
           'code/cubic_linear_limit.py','code/cubic_linear_dynamics.py']
    target.write_text(json.dumps(dict(status='passed',records=rows,finite_checks=checks,precisions=[80,120],
        max_absolute_precision_change=mp.nstr(max(changes),30),elapsed_seconds=time.time()-start,
        scope='All p+1 input ranks and full terminating Racah series, no numerical rank cutoff. B,L,O block and (1+u+v)^2-weighted positive A norm. Does not check S--O or the full memory spectrum. Floating evaluations are not interval certificates.',
        environment=dict(python=platform.python_version(),mpmath=mp.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')


if __name__=='__main__':main()
