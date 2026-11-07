"""Central fixed-charge occupation ratio with an explicit log remainder."""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import hashlib,json,platform,time
import mpmath as mp

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-central-spectrum'

def evaluate(N,m,L,t,dps):
    mp.mp.dps=dps;s=L-t
    assert max(t,s,L)<min(m,N-m,N)
    falling=lambda a,k:mp.loggamma(a+1)-mp.loggamma(a-k+1)
    logR=L*mp.log(2)+falling(m,t)+falling(N-m,s)-falling(N,L)
    approx=t*mp.log(mp.mpf(2*m)/N)+s*mp.log(mp.mpf(2*(N-m))/N)
    approx-=mp.mpf(t*(t-1))/(2*m)+mp.mpf(s*(s-1))/(2*(N-m))
    approx+=mp.mpf(L*(L-1))/(2*N)
    def remainder(a,k):
        return mp.mpf(k*(k-1)*(2*k-1))/(12*a*a*(1-mp.mpf(max(k-1,0))/a)) if k else mp.mpf(0)
    error=remainder(m,t)+remainder(N-m,s)+remainder(N,L)
    delta=mp.mpf(2*m-N)/(2*N)
    upper=L*mp.log(1+2*abs(delta))+mp.mpf(L*(L-1))/(2*(N-L+1))
    tol=mp.mpf(10)**(-dps+12)
    assert abs(logR-approx)<=error+tol
    assert logR<=upper+tol
    ratio=mp.exp(logR)
    assert abs(ratio-1)<=mp.exp(upper)*abs(logR)+tol
    return [ratio,logR,approx,error,upper,abs(logR-approx)]

def main():
    OUT.mkdir(exist_ok=True);target=OUT/'trace-ratio.json'
    if target.exists():raise FileExistsError(target)
    start=time.time();finite=[]
    for N,m,L in [(12,5,3),(14,8,4),(20,8,6),(24,13,7)]:
        for t in range(L+1):
            exact=F(2**L*comb(N-L,m-t),comb(N,m))
            vals=evaluate(N,m,L,t,120)
            assert abs(vals[0]-mp.mpf(exact.numerator)/exact.denominator)<mp.mpf('1e-100')
            finite.append(dict(N=N,m=m,L=L,t=t,exact_ratio=str(exact)))
    rows=[];changes=[]
    for lam in [F(1,10),F(1,4),F(1,2)]:
        for p in [16,32,64,128]:
            N=int(p*p/lam)
            for x in [-2,0,2]:
                m=N//2+x*p
                for slots in [2,4,8]:
                    L=slots*p
                    for imbalance in [0,2,6]:
                        t=(L+imbalance)//2;lo=evaluate(N,m,L,t,90);hi=evaluate(N,m,L,t,130)
                        change=max(abs(a-b) for a,b in zip(lo,hi));assert change<mp.mpf('1e-75');changes.append(change)
                        rows.append(dict(N=N,p=p,m=m,x=x,lambda_value=str(lam),L=L,t=t,imbalance=imbalance,
                            values=[mp.nstr(v,75) for v in hi],p_times_ratio_minus_one=mp.nstr(p*(hi[0]-1),35)))
            print('central trace ratio',str(lam),p,'passed',flush=True)
    names=['code/memory_central_trace.py']
    data=dict(status='passed',finite_exact_checks=finite,records=rows,
        labels=['ratio','log_ratio','first_falling_factorial_approximation','analytic_log_remainder_bound','global_log_upper_bound','observed_log_remainder'],
        precisions=[90,130],max_precision_change=mp.nstr(max(changes),30),elapsed_seconds=time.time()-start,
        scope='Exact combinatorial occupation ratio and explicit analytic falling-factorial remainder. 324 balanced central configurations. Uniform marked-graph transfer follows from the note and prior exponential overlap moments; high precision is not interval arithmetic.',
        environment=dict(python=platform.python_version(),mpmath=mp.__version__),
        input_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names})
    target.write_text(json.dumps(data,indent=2)+'\n')

if __name__=='__main__':main()
