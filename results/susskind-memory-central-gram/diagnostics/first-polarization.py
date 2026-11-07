"""Neutral-moment drift for the adjacent-charge polarization proof."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,platform,time
import sympy as s
import mpmath as mp
from cubic_linear_dynamics import weights,crossing
from cubic_linear_limit import distribution

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-central-gram'

def neutral_pair(N,m,a,b):
    cf=F(m*(N-m)*(N+1),N)
    return (-2*a*b*(1+crossing(N,m,a,b))+F(a*(N+1-a)*b*(N+1-b),cf))/N

def neutral_moment(N,p,m):
    w=list(reversed(distribution(N,p,m)))
    ma=mp.fsum(w[a]*a for a in range(p+1))
    mx=mp.fsum(w[a]*a*(N+1-a) for a in range(p+1))
    z=[mp.mpf(1)]*(p+1);h=mp.mpf(1);parts=[]
    for j in range(p+1):
        za=mp.fsum(w[a]*a*z[a] for a in range(p+1))
        parts.append(h*za*za)
        if j<p:
            z=[v*(-a+j)*(a-N-1+j) for a,v in enumerate(z)]
            h/=(-m+j)*(m-N+j)*(-N-1+j)*(j+1)
    cf=mp.mpf(m*(N-m)*(N+1))/N
    return (-2*(ma*ma+mp.fsum(parts))+mx*mx/cf)/N

def evaluate(N,p,m,dps):
    mp.mp.dps=dps
    upper=neutral_moment(N,p,m);lower=neutral_moment(N,p,m-1)
    drift=lower-upper
    return [upper,lower,drift,N*p*drift,(m-2*p)*drift/2]

def main():
    OUT.mkdir(exist_ok=True);target=OUT/'neutral-drift.json'
    if target.exists():raise FileExistsError(target)
    start=time.time();e,lam,x,u,k=s.symbols('e lam x u k')
    p=1/e;N=p*p/lam;m=N/2+x*p;n=N-m+1;a=p-u;xa=a*(N+1-a)
    zeta=(m-p)*n*n/((n+p)*(m*n-xa))
    numerator=(n+p)*xa-p*n*(N+1)
    assert s.factor(zeta-1-numerator/((n+p)*(m*n-xa)))==0
    score=s.factor(s.limit(N*(zeta-1),e,0))
    assert score==4*(lam-u)
    r=2*p-k;ratio=(m-p)*(n+p)/(m*n);ell2=m*(m-r)*(n-r)/n
    desired=m-2*p+k-2*lam
    correction=s.factor(s.limit(p*(ratio**2*ell2-desired**2)/(2*m),e,0))
    exact_ratios=0;finite=[]
    for NN in range(6,18):
        for pp in range(1,min(4,NN//2)):
            for mm in range(pp+1,NN-pp+1):
                w=weights(NN,mm,pp);wm=weights(NN,mm-1,pp);nn=NN-mm+1
                for aa in range(pp+1):
                    zz=F((mm-pp)*nn*nn,(nn+pp)*(mm*nn-aa*(NN+1-aa)))
                    assert wm[aa]==w[aa]*zz;exact_ratios+=1
    for NN,pp,mm in [(6,1,3),(8,2,4),(10,2,6),(12,3,5)]:
        w=weights(NN,mm,pp)
        ex=sum(w[a]*w[b]*neutral_pair(NN,mm,a,b) for a in range(pp+1) for b in range(pp+1))
        mp.mp.dps=120;got=neutral_moment(NN,pp,mm)
        assert abs(got-mp.mpf(ex.numerator)/ex.denominator)<mp.mpf('1e-100')
        finite.append(dict(N=NN,p=pp,m=mm,exact_neutral_moment=str(ex)))
    rows=[];refinements=[]
    for ll in [F(1,10),F(1,4),F(1,2)]:
        for pp in [16,32,64,128]:
            NN=int(pp*pp/ll)
            for xx in [-2,0,2]:
                mm=NN//2+xx*pp;lo=evaluate(NN,pp,mm,80);hi=evaluate(NN,pp,mm,120)
                change=max(abs(aa-bb) for aa,bb in zip(lo,hi))
                assert change<mp.mpf('1e-60');refinements.append(change)
                rows.append(dict(N=NN,p=pp,m=mm,x=xx,lambda_value=str(ll),
                    values=[mp.nstr(v,75) for v in hi],absolute_precision_change=mp.nstr(change,30)))
            print('neutral drift',str(ll),pp,'passed',flush=True)
    names=['code/memory_central_polarization.py','code/cubic_linear_dynamics.py','code/cubic_linear_limit.py']
    data=dict(status='passed',symbolic=dict(weight_ratio_numerator=str(numerator),weight_score_limit=str(score),
        first_transfer_correction=str(correction)),exact_weight_ratios=exact_ratios,finite_checks=finite,
        labels=['F_m','F_m_minus_1','adjacent_difference','N_p_times_difference','half_h_times_difference'],
        records=rows,precisions=[80,120],max_precision_change=mp.nstr(max(refinements),30),
        scope='Neutral Racah and rank-weight identities only; no candidate charged Casimir formula. Every input rank and terminating-series term retained. Uniform O(1/(N p)) proof is in the companion note; floating arithmetic is not interval certified.',
        environment=dict(python=platform.python_version(),sympy=s.__version__,mpmath=mp.__version__),elapsed_seconds=time.time()-start,
        input_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names})
    target.write_text(json.dumps(data,indent=2)+'\n')

if __name__=='__main__':main()
