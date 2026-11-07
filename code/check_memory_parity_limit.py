"""High-precision finite Casimir-deficit moments, including the parity sign."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,platform
import mpmath as mp
from cubic_linear_dynamics import crossing
from memory_parity_moments import signed_dx

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-parity-gram'

def ordinary_dx(N,a,b):
    m=N//2;x=lambda j:j*(N+1-j)
    Cf=F(m*(N-m)*(N+1),N)
    M1=x(a)+x(b)-F(x(a)*x(b),1)/Cf
    a21=F(-2*(N-1),m*(N-m)*(N+1))
    a22=F(2*(N-1)*(N-2),m*(m-1)*(N-m)*(N-m-1)*N*(N+1)*2)
    # The coefficient is (-2)_2(1-N)_2 / [(-m)_2(m-N)_2(-N-1)_2 2!].
    F2=(crossing(N,m,a,2)*crossing(N,m,b,2)-1-a21*M1)/a22
    EX2=F2+N*M1
    return [F(1),x(a+b)-M1,x(a+b)**2-2*x(a+b)*M1+EX2]

def main():
    target=OUT/'limit-checks.json'
    if target.exists():raise FileExistsError(target)
    # Independent finite full-distribution check of the unsigned recurrence.
    from memory_parity_moments import product_table
    for N,a,b in [(6,1,2),(8,2,3),(10,3,3)]:
        _,_,pi=product_table(N);x=lambda j:j*(N+1-j)
        assert ordinary_dx(N,a,b)==[sum(pr*(x(a+b)-x(r))**k for r,pr in enumerate(pi[a,b])) for k in range(3)]
    rows=[]
    for lam in [F(1,10),F(1,4),F(1,2)]:
        for p in [16,32,64,128]:
            N=int(p*p/lam);assert F(p*p,N)==lam and N%2==0
            for u,v in [(0,0),(1,2)]:
                a,b=p-u,p-v
                signed=signed_dx(N,a,b);unsigned=ordinary_dx(N,a,b)
                runs=[]
                for precision in [80,120]:
                    mp.mp.dps=precision
                    cv=lambda f:mp.mpf(f.numerator)/f.denominator
                    theta=2*cv(lam);q=mp.exp(-2*theta)
                    values=[cv(f/F(N**k)) for k,f in enumerate(signed)]
                    evens=[cv((x+y)/F(N**k)) for k,(x,y) in enumerate(zip(signed,unsigned))]
                    target_signed=[q,-theta*q,(theta**2-theta)*q]
                    target_even=[1+q,theta*(1-q),theta**2*(1+q)+theta*(1-q)]
                    runs.append(values+evens)
                assert all(abs(x-y)<mp.mpf('1e-70') for x,y in zip(*runs))
                rows.append(dict(N=N,p=p,lambda_value=str(lam),input_deficits=[u,v],
                    signed_scaled_delta_casimir_moments=list(map(lambda x:mp.nstr(x,60),values)),
                    signed_limits=list(map(lambda x:mp.nstr(x,60),target_signed)),
                    signed_absolute_deviations=[mp.nstr(x-y,40) for x,y in zip(values,target_signed)],
                    even_weighted_proxy_moments=list(map(lambda x:mp.nstr(x,60),evens)),
                    even_weighted_limits=list(map(lambda x:mp.nstr(x,60),target_even))))
    names=['code/check_memory_parity_limit.py','code/memory_parity_moments.py','code/cubic_linear_dynamics.py']
    result=dict(status='passed',records=rows,precisions=[80,120],
        finite_quantity='Delta x/N, not S itself. Their first two signed moments have the same limit by the earlier third-moment tail bound.',
        scope='Exact finite Fraction formulas followed by decimal evaluation; all deviations retained, no interval certificate.',
        environment=dict(python=platform.python_version(),mpmath=mp.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print('passed',len(rows),'finite cases; last signed deviations',rows[-1]['signed_absolute_deviations'])

if __name__=='__main__':main()
