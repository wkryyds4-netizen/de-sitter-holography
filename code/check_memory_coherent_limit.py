"""Stable coherent cancellations and complete input-rank sums at two precisions."""
from pathlib import Path
from fractions import Fraction as F
import hashlib, json, platform
import mpmath as mp
import sympy as s
from memory_coherent_cross import charged, exchanged, mixed_moments, cross_BX_A
from cubic_linear_dynamics import crossing, weights
from cubic_linear_limit import distribution

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-coherent-cross'


def full_sum(lam,p,precision):
    mp.mp.dps=precision
    lam=mp.mpf(lam);N=int(mp.nint(p*p/lam));m=N//2
    assert N%2==0 and mp.mpf(p*p)/N==lam
    omega=list(reversed(distribution(N,p,m)))
    ma=mp.fsum(a*omega[a] for a in range(p+1))
    mx=mp.fsum(a*(N+1-a)*omega[a] for a in range(p+1))
    cf=mp.mpf(m*m)*(N+1)/N;gain=mp.mpf((m-p)**2)/m;theta=mp.mpf(p*p)/m
    z=[mp.mpf(1)]*(p+1);h=mp.mpf(1)
    terms=[[] for _ in range(6)]
    for k in range(p+1):
        z0=mp.fsum(omega[a]*z[a] for a in range(p+1))
        z1=mp.fsum(omega[a]*a*z[a] for a in range(p+1))
        hd=h*mp.mpf(k)/(m-k)
        ll=2*(ma*z0-z1)*z0
        for arr,value in zip(terms,[h*z0*z0,h*z1*z1,hd*z0*z0,hd*z1*z1,h*ll,hd*ll]): arr.append(value)
        if k<p:
            for a in range(p+1): z[a]*=(-a+k)*(a-N-1+k)
            h/=(-m+k)**2*(-N-1+k)*(k+1)
    chi,abchi,diff,abdiff,lchi,ldiff=map(mp.fsum,terms)
    ux=-2*ma*ma+mx*mx/cf
    ud=mx*mx/(mp.mpf(m)*(m-1)*(N+1))
    bparts=[theta,theta*chi,gain*diff]
    xparts=[theta*ux,-2*theta*abchi,gain*ud,gain*m*diff,-2*gain*abdiff]
    vals=[mp.fsum(bparts),mp.fsum(xparts)/N,theta*lchi+gain*ldiff]
    q=mp.exp(-4*lam);targets=[2*lam*(1-q),4*lam*lam*(1+q)+2*lam*(1-q),mp.mpf(0)]
    cancellation=max(mp.fsum(map(abs,bparts))/abs(vals[0]),mp.fsum(map(abs,xparts))/abs(vals[1]*N))
    return dict(N=N,p=p,lambda_value=str(lam),precision=precision,
                values=[mp.nstr(x,65) for x in vals],limits=[mp.nstr(x,65) for x in targets],
                absolute_deviations=[mp.nstr(x-y,40) for x,y in zip(vals,targets)],
                cancellation_factor=mp.nstr(cancellation,30))


def main():
    target=OUT/'limit-checks.json'
    if target.exists(): raise FileExistsError(target)
    # Contiguous identity reduces the first signed moment to two neutral-type sums.
    checks=0
    for N in range(4,22,2):
        m=N//2
        for a in range(m):
            for b in range(m):
                xi=exchanged(N,a,b);chi=crossing(N,m,a,b)
                direct=mixed_moments(N,a,b)[1][1]
                assert direct == -2*a*b*xi+m*(xi-chi)
                checks+=1
    m=s.symbols('m',positive=True);N=2*m
    cf=m*m*(N+1)/N;ct=m*(m-1)*(N+1)/(N-1)
    assert s.factor(1/ct-1/cf-1/(m*(m-1)*(N+1)))==0
    k,b=s.symbols('k b');ratio=(b-k)*(b-N-1)/(b*(b-N+k-1))
    T_over_xi=(N+1-b-k)/(N+1-b)
    assert s.factor(b*(N+1-b)/(N+1-2*b)*(ratio-1)*T_over_xi+k)==0
    # Full-sum implementation versus exact finite pair summation.
    for NN,pp in [(8,2),(12,3),(16,4)]:
        ww=weights(NN,NN//2,pp);mean=sum(a*ww[a] for a in range(pp+1))
        expected=[F(0)]*3
        for aa in range(pp+1):
            for bb in range(pp+1):
                bv,xv=cross_BX_A(NN,pp,aa,bb)
                expected[0]+=ww[aa]*ww[bb]*bv
                expected[1]+=ww[aa]*ww[bb]*xv
                expected[2]+=ww[aa]*ww[bb]*(2*mean-aa-bb)*bv
        out=full_sum(str(F(pp*pp,NN).numerator/F(pp*pp,NN).denominator),pp,120)
        assert all(abs(mp.mpf(v)-mp.mpf(e.numerator)/e.denominator)<mp.mpf('1e-60') for v,e in zip(out['values'],expected))
    rows=[];maxchange=mp.mpf(0)
    for lam in ['0.1','0.25','0.5']:
        for p in [16,32,64,128]:
            low=full_sum(lam,p,80);high=full_sum(lam,p,120)
            change=max(abs(mp.mpf(x)-mp.mpf(y)) for x,y in zip(low['values'],high['values']))
            assert change<mp.mpf('1e-60');maxchange=max(maxchange,change)
            rows.append(high)
            print('full rank cross',lam,p,high['absolute_deviations'],flush=True)
    names=['code/check_memory_coherent_limit.py','code/memory_coherent_cross.py',
           'code/cubic_linear_dynamics.py','code/cubic_linear_limit.py']
    result=dict(status='passed',contiguous_identities=checks,symbolic_identities=2,
                independent_exact_full_sums=3,records=rows,precisions=[80,120],
                max_absolute_precision_change=mp.nstr(maxchange,30),
                quantities=['<B,W2 A>/C^2','<X_N/N,W2 A>/C^2','<L_input,W2 A>/C^2'],
                scope='All p+1 input ranks, all terminating-series terms; no cutoff. X_N/N is a finite Casimir proxy, not S_N. Limit deviations and cancellation retained. No interval certificate.',
                environment=dict(python=platform.python_version(),mpmath=mp.__version__,sympy=s.__version__),
                input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__': main()
