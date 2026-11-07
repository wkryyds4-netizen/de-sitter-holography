"""High-precision Jacobi expansion checks outside its auxiliary spectral interval."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,platform
import mpmath as mp
from memory_wick_rank import mu

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-uniform-bound'


def evaluate(N,p,m,s,x,dps):
    mp.mp.dps=dps;s=mp.mpf(s);x=mp.mpf(x)
    degree=64;order=150
    coeff=[F(0)]+[sum((mu(N,p,m,j*p) for j in range(n)),F(0)) for n in range(1,order+3)]
    a=[mp.mpf(v.numerator)/v.denominator for v in coeff]
    values=[mp.mpf(1),x]
    for n in range(1,degree):values.append(x*values[-1]-a[n]*values[-2])
    w=[mp.mpf(0)]*(order+2);w[0]=1
    c=[mp.mpc(0)]*(degree+1);factor=mp.mpc(1)
    for ell in range(order+1):
        for n in range(min(ell,degree)+1):c[n]+=factor*w[n]
        w=[(w[n-1] if n else 0)+a[n+1]*w[n+1] for n in range(order+1)]+[mp.mpf(0)]
        factor*=1j*s/(ell+1)
    actual=mp.fsum(c[n]*values[n] for n in range(degree+1))
    q=mu(N,p,m,p);A=1/(1-mp.mpf(q.numerator)/q.denominator)
    return dict(N=N,p=p,m=m,time=str(s),spectral_argument=str(x),C=1,
                auxiliary_spectral_bound=str(2*mp.sqrt(A)),outside_spectral_bound=abs(x)>2*mp.sqrt(A),
                retained_polynomials=degree+1,Jacobi_time_series_order=order,precision=dps,
                real=str(actual.real),imag=str(actual.imag),absolute_error=str(abs(actual-mp.exp(1j*s*x))))


if __name__=='__main__':
    dest=OUT/'jacobi-entire-checks.json'
    if dest.exists():raise FileExistsError(dest)
    exact=0
    for N,p,m in [(16,1,8),(32,2,16),(64,4,32),(24,2,6),(32,3,24)]:
        q=mu(N,p,m,p)
        for j in range(21):
            assert 0<=mu(N,p,m,j*p)<=q**j
            exact+=1
    rows=[]
    for N,p,m in [(16,1,8),(32,2,16),(64,4,32)]:
        for x in ('-6','6'):
            lo=evaluate(N,p,m,'.6',x,80);hi=evaluate(N,p,m,'.6',x,120)
            difference=abs(mp.mpc(lo['real'],lo['imag'])-mp.mpc(hi['real'],hi['imag']))
            assert hi['outside_spectral_bound']
            assert difference<mp.mpf('1e-65')
            assert mp.mpf(hi['absolute_error'])<mp.mpf('1e-40')
            hi['precision_difference']=str(difference);rows.append(hi)
    names=['code/check_memory_jacobi_entire.py','code/memory_wick_rank.py',
           'code/cubic_linear_dynamics.py','code/complex_syk_dictionary.py']
    result=dict(exact_geometric_cases=exact,records=rows,
                scope='Auxiliary finite-parameter Jacobi expansion evaluated at x=+-6 outside its compact spectral bound, with 80/120 digit repeat. This checks cancellation and whole-plane evaluation, not a substitute for the normal-convergence/ODE proof; no interval-rounding certification.',
                environment=dict(python=platform.python_version(),mpmath=mp.__version__),
                input_sha256={x:hashlib.sha256((ROOT/x).read_bytes()).hexdigest() for x in names})
    with dest.open('x',encoding='utf-8') as f:f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='passed',exact_geometric_cases=exact,whole_plane_rows=len(rows))))
