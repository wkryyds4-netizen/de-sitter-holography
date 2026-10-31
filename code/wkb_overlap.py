"""CSTY reduced recoil kernel: finite-nu confluence and prefactor diagnostics.

No finite-G amplitude, physical mass packet, or interval certificate is computed.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, platform, sys
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
mp.mp.dps=60

def scaled_kernel(Q,nu=None):
    Q=mp.mpf(Q)
    # theta=u/sqrt(Q); factor exp(-2Q) keeps the sampled scale moderate.
    end=mp.pi*mp.sqrt(Q)/2
    points=sorted(set([mp.mpf(0),end]+[mp.mpf(v) for v in [1,3,8,16] if v<end]))
    def integrand(u):
        z=mp.cos(u/mp.sqrt(Q))
        if nu is None:
            phase=2*Q*z; amp=1
        else:
            r=Q*z/nu
            phase=2*nu*mp.asin(r); amp=1/mp.sqrt(1-r*r)
        return (mp.exp(phase-2*Q)-mp.exp(-phase-2*Q))*amp/(4*Q*mp.sqrt(Q))
    return mp.quad(integrand,points)

def upper_log(Q,nu):
    r=mp.mpf(Q)/nu
    return mp.log(mp.asin(r)/r)-mp.log(1-r*r)/2+2*nu*(mp.asin(r)-r)

def reduced_I(A,t):
    # y=A*sinh(t)+sqrt(A*cosh(t))*u, with saddle exponent subtracted.
    A=mp.mpf(A); t=mp.mpf(t); c=mp.cosh(t); s=mp.sinh(t)
    def f(u):
        z=s+mp.sqrt(c/A)*u
        phase=2*A*(-z*mp.asinh(z)+mp.sqrt(1+z*z)-1+z*t-(c-1))
        return mp.exp(phase)
    return mp.quad(f,[-mp.inf,-8,-3,0,3,8,mp.inf])*mp.sqrt(A*c)

def run():
    checks=[]; samples=[]; critical=[]; prefactors=[]
    def check(name,error,tol):
        checks.append(dict(name=name,error=str(error),tolerance=str(tol),passed=bool(error<=tol)))
    for nu,Q in [(2,.1),(2,1),(5,2),(10,5),(30,10),(100,20)]:
        nu=mp.mpf(nu);Q=mp.mpf(str(Q))
        F=scaled_kernel(Q,nu);H=scaled_kernel(Q)
        ref=mp.hyper([1,1-1j*nu,1+1j*nu],[mp.mpf('1.5')]*2,(Q/nu)**2)*mp.exp(-2*Q)
        href=mp.hyper([1],[mp.mpf('1.5')]*2,Q**2)*mp.exp(-2*Q)
        check(f'finite-nu-integral-{nu}-{Q}',abs(F/ref-1),mp.mpf('1e-48'))
        check(f'limit-integral-{Q}',abs(H/href-1),mp.mpf('1e-48'))
        lr=mp.log(F/H); ub=upper_log(Q,nu)
        check(f'positive-bound-{nu}-{Q}',max(-lr,lr-ub,mp.mpf(0)),mp.mpf('1e-48'))
        samples.append(dict(nu=str(nu),Q=str(Q),log_ratio=str(lr),upper_log=str(ub)))
    for nu in [1000,10000,1000000]:
        for c in [mp.mpf('.5'),mp.mpf(1),mp.mpf(2)]:
            Q=c*mp.mpf(nu)**(mp.mpf(2)/3)
            F=scaled_kernel(Q,nu); H=scaled_kernel(Q)
            with mp.workdps(85):
                high=mp.log(scaled_kernel(Q,nu)/scaled_kernel(Q))
            lr=mp.log(F/H); r=Q/nu
            saddle=2*nu*(mp.asin(r)-r)-mp.log(1-r*r)/4
            check(f'critical-precision-{nu}-{c}',abs(lr-high),mp.mpf('1e-45'))
            check(f'critical-upper-{nu}-{c}',max(lr-upper_log(Q,nu),-lr,mp.mpf(0)),mp.mpf('1e-45'))
            critical.append(dict(nu=nu,c=str(c),Q=str(Q),log_ratio=str(lr),critical_limit=str(c**3/3),saddle_log_ratio=str(saddle),saddle_residual=str(lr-saddle)))
    for A in [10,40,160]:
        i0=reduced_I(A,0); h0=scaled_kernel(A)
        for t in [mp.mpf('.5'),mp.mpf(1)]:
            c=mp.cosh(t); it=reduced_I(A,t); ht=scaled_kernel(A*c)
            # common exponential cancels in this normalized comparison.
            ratio=(it/i0)/(ht/h0)
            with mp.workdps(80):
                hi=(reduced_I(A,t)/reduced_I(A,0))/(scaled_kernel(A*mp.cosh(t))/scaled_kernel(A))
            check(f'prefactor-precision-{A}-{t}',abs(ratio/hi-1),mp.mpf('1e-40'))
            prefactors.append(dict(A=A,t=str(t),normalized_I_over_H=str(ratio),leading_cosh_squared=str(c*c),relative_residual=str(ratio/(c*c)-1)))
    files=['code/wkb_overlap.py','sources/chen-stanford-tang-yang-2607.14042v1.html','sources/chandra-liska-post-2610.30247v1.html']
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='passed' if all(c['passed'] for c in checks) else 'failed',scope='reduced pure-recoil kernels and action-only inverse integral; empirical arbitrary-precision checks, not finite-G or interval certification',environment=dict(python=sys.version,platform=platform.platform(),mpmath=mp.__version__,dps=60,crosscheck_dps=[80,85]),checks=checks,samples=samples,critical_scaling=critical,prefactors=prefactors,hashes={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files})
    out=ROOT/'results/stage4/wkb-overlap';out.mkdir(parents=True,exist_ok=True)
    (out/'diagnostics.json').write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=data['status'],checks=len(checks),failed=sum(not c['passed'] for c in checks))))
    return data

if __name__=='__main__': run()
