"""Joint endpoint estimates and central absolute-conditioning diagnostics.

Finite-point checks, not interval certificates or a signed-amplitude calculation.
The old kernel implementation is imported without modification.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, platform, sys
import mpmath as mp
from recoil_bessel_limit import log_q

ROOT = Path(__file__).resolve().parents[1]


def log_u(t, K, alpha):
    return (-2j*mp.pi*t + 2*t*mp.log(K/2)
            + mp.loggamma(mp.mpf('.5')-t+alpha/2)
            + mp.loggamma(mp.mpf('.5')-t-alpha/2))


def envelope(K, alpha):
    return K*(mp.besselk(alpha, K)+mp.pi*mp.besseli(alpha, K))


def F(y):
    return mp.sqrt(1+y*y)-y*mp.asinh(y)


def B(y):
    return y*mp.acosh(y)-mp.sqrt(y*y-1)


def h(y, r):
    return ((r+y/2)*mp.log(2*(r+y/2))
            +(r-y/2)*mp.log(2*(r-y/2))-2*r)


def vertical_log_ratio(N, s, K, nu, alpha):
    """Independent exact real-product expression; no log_q evaluation."""
    v=nu/2; M=nu*nu/(2*K); bm=M-nu; cp=M+nu
    def lc(x):
        return mp.log(mp.cosh(mp.pi*x))
    # Four negative-integer Gamma factors, including their w^-1 terms.
    out=-mp.pi*s+2*(mp.log(v/(v+s))
        +mp.log(mp.sinh(mp.pi*v))-mp.log(mp.sinh(mp.pi*(v+s))))
    out+=2*mp.fsum(mp.log((j*j+v*v)/(j*j+(v+s)**2)) for j in range(1,N+1))
    # Positive half-integer Gamma and reciprocal negative half-integer Gamma.
    out+=(lc(bm)-lc(bm-s)+lc(cp+s)-lc(cp))/2
    out+=mp.fsum(mp.log(((j+mp.mpf('.5'))**2+(bm-s)**2)
                         /((j+mp.mpf('.5'))**2+bm*bm))
                 +mp.log(((j+mp.mpf('.5'))**2+(cp+s)**2)
                         /((j+mp.mpf('.5'))**2+cp*cp))
                 for j in range(N))/2
    for a in [N+mp.mpf('.5')+alpha/2,N+mp.mpf('.5')-alpha/2]:
        out+=mp.re(mp.loggamma(a+1j*s))-mp.loggamma(a)
    return out


def window(K, order, dps):
    """Integrate |U| and |Q| at y=2, scaled by E; unit-or-shorter cells."""
    with mp.workdps(dps):
        K=mp.mpf(K); alpha=2*K; nu=32*K**mp.mpf('1.5')
        xstar=-mp.sqrt(3)*K/2; half=mp.sqrt(K)/4
        lo=xstar-half; hi=xstar+half
        cuts=[lo]+[mp.mpf(j) for j in range(int(mp.floor(lo))+1,int(mp.ceil(hi)))]+[hi]
        nodes,weights=mp.gauss_quadrature(order,'legendre')
        loge=mp.log(envelope(K,alpha)); vals=[mp.mpf(0),mp.mpf(0)]
        for a,b in zip(cuts[:-1],cuts[1:]):
            for z,w in zip(nodes,weights):
                t=(a+b)/2+(b-a)*z/2+1j
                vals[0]+=w*(b-a)/2*mp.exp(mp.re(log_u(t,K,alpha))-loge)
                vals[1]+=w*(b-a)/2*mp.exp(mp.re(log_q(t,K,nu,alpha))-loge)
        return vals,dict(K=str(K),y='2',alpha=str(alpha),nu=str(nu),rho='0.0009765625',
                         height='1',lo=str(lo),hi=str(hi),cells=len(cuts)-1,order=order,dps=dps,
                         U_absolute_over_E=str(vals[0]),Q_absolute_over_E=str(vals[1]),
                         U_log_rate=str(mp.log(vals[0])/K),Q_log_rate=str(mp.log(vals[1])/K))


def run():
    mp.mp.dps=75
    checks=[]; vertical=[]; bases=[]; positive=[]; central=[]; local=[]
    def check(name,error,tolerance):
        checks.append(dict(name=name,error=str(abs(error)),tolerance=str(tolerance),
                           passed=bool(abs(error)<=tolerance)))
    for ki in [8,32,128,512]:
        K=mp.mpf(ki); nu=32*K**mp.mpf('1.5'); M=nu*nu/(2*K)
        for yi in [0,1,2]:
            y=mp.mpf(yi); alpha=K*y; N=int(mp.ceil(K*mp.sqrt(1+y*y)/2))
            lam=2*mp.pi-N*(1/(M-nu)+1/(M+nu)); rho=K**3/nu**2
            logbase=mp.re(log_q(-N,K,nu,alpha)); loge=mp.log(envelope(K,alpha))
            check(f'lambda-{ki}-{yi}',max(0,mp.pi-lam),0)
            samples=[mp.mpf(0),mp.mpf('.1'),mp.mpf(1),mp.mpf(10),K,nu,M-nu,M,M+nu,2*M]
            for s in samples:
                exact=mp.re(log_q(-N-1j*s,K,nu,alpha))-logbase
                independent=vertical_log_ratio(N,s,K,nu,alpha)
                check(f'vertical-identity-{ki}-{yi}-{s}',exact-independent,mp.mpf('1e-55'))
                check(f'vertical-bound-{ki}-{yi}-{s}',max(0,exact+lam*s),mp.mpf('1e-55'))
                vertical.append(dict(K=str(K),y=str(y),alpha=str(alpha),nu=str(nu),N=N,
                                     s=str(s),log_ratio=str(exact),product_log_ratio=str(independent),
                                     lambda_value=str(lam),log_upper=str(-lam*s)))
            bases.append(dict(K=str(K),y=str(y),alpha=str(alpha),nu=str(nu),N=N,rho=str(rho),
                              scaled_base=str(mp.exp(logbase-loge)*mp.sqrt(K)),
                              exponent_removed_base=str(mp.exp(logbase+K*F(y))),
                              U_exponent_removed=str(mp.exp(mp.re(log_u(-N,K,alpha))+K*F(y))),
                              multiplier_real_over_rho=str((logbase-mp.re(log_u(-N,K,alpha)))/rho)))
            # Y=2 gives a common R=12K. Samples are NOT a supremum certificate.
            for dx in [mp.mpf(0),mp.mpf('.5'),mp.mpf(1)]:
                x=12*K+dx; t=x+1j; lq=mp.re(log_q(t,K,nu,alpha)); lu=mp.re(log_u(t,K,alpha))
                reflected=(-2j*mp.pi*t+2*t*mp.log(K/2)+2*mp.log(mp.pi)
                    -mp.log(mp.cos(mp.pi*(t-alpha/2)))-mp.log(mp.cos(mp.pi*(t+alpha/2)))
                    -mp.loggamma(mp.mpf('.5')+t-alpha/2)-mp.loggamma(mp.mpf('.5')+t+alpha/2))
                check(f'reflection-{ki}-{yi}-{dx}',mp.exp(log_u(t,K,alpha)-reflected)-1,mp.mpf('1e-55'))
                positive.append(dict(K=str(K),y=str(y),nu=str(nu),x=str(x),
                                     Q_exponent_removed_log=str(lq+K*h(y,x/K)),
                                     U_exponent_removed_log=str(lu+K*h(y,x/K)),
                                     log_Q_over_E=str(lq-loge)))
            # Check the phase-preserving expansion modulo 2*pi*i.
            for t in [-N,mp.mpc(-K*mp.sqrt(3)/2,1),mp.mpc(12*K,1)]:
                residual=mp.exp(log_q(t,K,nu,alpha)-log_u(t,K,alpha)-4j*(t*t-t)/nu)
                local.append(dict(K=str(K),y=str(y),nu=str(nu),t=str(t),rho=str(rho),
                                  phase_removed_multiplier=str(residual),
                                  principal_log_over_rho=str(mp.log(residual)/rho)))
        print(f'endpoint samples K={ki}',flush=True)
    delta=B(mp.mpf(2))+F(mp.mpf(2))
    for K in [8,16,32,64,128]:
        low,_=window(K,20,55); high,row=window(K,32,75)
        changes=[abs(a/b-1) for a,b in zip(low,high)]
        for idx,name in enumerate(['U','Q']):
            check(f'window-order-{name}-{K}',changes[idx],mp.mpf('1e-16'))
        row['quadrature_changes']=[str(x) for x in changes]
        if K==128:
            precise,_=window(K,32,95)
            pc=[abs(a/b-1) for a,b in zip(precise,high)]
            for idx,name in enumerate(['U','Q']):
                check(f'window-precision-{name}-{K}',pc[idx],mp.mpf('1e-55'))
            row['precision_changes']=[str(x) for x in pc]
        central.append(row)
        print(f'central absolute window K={K}',flush=True)
    inputs=['code/recoil_bessel_endpoints.py','code/recoil_bessel_limit.py',
            'notes/stage4-recoil-bessel-cancellation.md','sources/dlmf-5.11.html']
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
              status='passed' if all(c['passed'] for c in checks) else 'failed',
              scope='endpoint identities and finite-point inequalities; central absolute norms only, not signed full amplitudes',
              environment=dict(python=sys.version,platform=platform.platform(),mpmath=mp.__version__,
                               dps=[55,75,95],orders=[20,32]),
              hierarchy='nu=32*K^(3/2), rho=1/1024 fixed: bound test, not rho -> 0 convergence',
              delta_y2=str(delta),checks=checks,vertical_rows=vertical,base_rows=bases,
              positive_base_samples=positive,central_rows=central,local_multiplier_samples=local,
              hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    out=ROOT/'results/stage4/recoil-bessel-endpoints';out.mkdir(parents=True,exist_ok=True)
    (out/'diagnostics.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=data['status'],checks=len(checks),delta=str(delta),
                         failed=[c['name'] for c in checks if not c['passed']])))
    if data['status']!='passed':raise SystemExit(1)


if __name__=='__main__':run()
