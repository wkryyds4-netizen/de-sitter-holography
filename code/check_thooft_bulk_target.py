"""Independent finite-part integrals, angular Gram entries and spectral arithmetic."""
from pathlib import Path
import hashlib, json, platform
import mpmath as mp

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-bulk-target'

def chebyu(n,z):
    a,b=mp.mpf(1),2*z
    if n==0:return a
    for k in range(1,n):a,b=b,2*z*b-a
    return b

def main():
    target=OUT/'checks.json'
    if target.exists(): raise FileExistsError(target)
    mp.mp.dps=65;rows=[];gram=[]
    for n in range(1,5):
        f=lambda y:mp.sqrt(1-(2*y-1)**2)*chebyu(n-1,2*y-1)
        for sx in ['.2','.5','.8']:
            x=mp.mpf(sx);v=f(x);der=[mp.diff(f,x,k) for k in range(1,5)]
            def regular(y):
                z=y-x
                if abs(z)<mp.mpf('1e-15'):
                    return der[1]/2+der[2]*z/6+der[3]*z*z/24
                return (f(y)-v-der[0]*z)/(z*z)
            got=-(mp.quad(regular,[0,x,1])+der[0]*mp.log((1-x)/x)-v/(x*(1-x)))
            expected=2*mp.pi*n*chebyu(n-1,2*x-1)
            error=abs(got-expected)
            assert error<mp.mpf('1e-38')
            rows.append(dict(n=n,x=sx,finite_part=str(got),expected=str(expected),absolute_error=str(error)))
    for m in range(1,7):
        for n in range(1,7):
            got=mp.quad(lambda t:mp.sin(m*t)*mp.sin(n*t)*mp.sin(t),[0,mp.pi])
            exact=mp.mpf(0) if (m+n)%2 else mp.mpf(1)/(1-(m-n)**2)-mp.mpf(1)/(1-(m+n)**2)
            assert abs(got-exact)<mp.mpf('1e-60')
            gram.append(dict(m=m,n=n,absolute_error=str(abs(got-exact))))
    moments=[]
    for d in [2,3]:
        for k in [1,2]:
            got=mp.quad(lambda p:2**(2*d)*p**(2*d-1-k)*mp.exp(-2*p)/mp.gamma(2*d),[0,1,mp.inf])
            exact=2**k*mp.gamma(2*d-k)/mp.gamma(2*d)
            assert abs(got-exact)<mp.mpf('1e-60')
            moments.append(dict(d=d,k=k,inverse_moment=str(got)))
    spec=json.loads((OUT/'spectrum.json').read_text())
    a=spec['chebyshev'][-1]['twice_lambda'];b=spec['fredholm'][-1]['twice_lambda']
    diff=max(abs(mp.mpf(str(x))/mp.mpf(str(y))-1) for x,y in zip(a,b))
    assert diff<mp.mpf('2e-9')
    masses=[mp.pi*mp.mpf(str(v)) for v in a]
    assert 2<masses[0]<3*mp.pi/4
    files=['code/check_thooft_bulk_target.py','code/thooft_bulk_target.py','results/susskind-bulk-target/spectrum.json']
    payload=dict(status='passed',finite_part=rows,angular_gram=gram,inverse_moments=moments,
        cross_representation_max_relative_error=str(diff),mass_squared_over_tension=[str(v) for v in masses],
        mass_ratios_to_ground=[str(mp.sqrt(v/masses[0])) for v in masses],
        dps=mp.mp.dps,environment=dict(python=platform.python_version(),mpmath=mp.__version__),
        error_scope='Empirical arbitrary-precision integral checks, not interval enclosures. Near-diagonal Taylor branch threshold 1e-15, degree 2 in the regularized integrand.',
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files})
    target.write_text(json.dumps(payload,indent=2)+'\n')
    print('12 finite-part integrals, 36 Gram integrals, 4 inverse moments; max two-method relative difference',diff)

if __name__=='__main__':main()
