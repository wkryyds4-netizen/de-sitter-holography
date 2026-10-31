"""Pure-recoil combined vertex diagnostics, not a complete 6j or physical OTOC."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, platform, sys
import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
mp.mp.dps = 65

def zs(ms):
    S = sum(ms)
    return [(1+1j*S)/2] + [(1-1j*(S-2*m))/2 for m in ms]

def logF(ms):
    S = sum(ms)
    return sum(mp.loggamma(1-1j*m) for m in ms)-2*sum(
        mp.loggamma((1-1j*(S-2*m))/2) for m in ms)

def logT(ms):
    S = sum(ms)
    return (sum(mp.log(mp.pi*m/mp.sinh(mp.pi*m)) for m in ms)/2
            +sum(mp.log(mp.cosh(mp.pi*(S-2*m)/2)/mp.pi) for m in ms))

def logE(ms):
    """Exponentially small factor after removing c sqrt(m1 m2 m3)."""
    S = sum(ms)
    return sum(mp.log1p(mp.exp(-mp.pi*(S-2*m))) for m in ms)-sum(
        mp.log1p(-mp.exp(-2*mp.pi*m)) for m in ms)/2

def ebound(ms):
    S=sum(ms)
    eta=mp.exp(-mp.pi*min(mp.re(S-2*m) for m in ms))
    delta=mp.exp(-2*mp.pi*min(mp.re(m) for m in ms))
    assert eta < 1 and delta < 1
    return 3*eta/(1-eta)+mp.mpf('1.5')*delta/(1-delta)

def logB(ms,k):
    """Ordinary Gamma ratio from two exact inverse-b shifts."""
    L=-1j*k
    return 2*mp.loggamma(L+zs(ms)[0])-sum(mp.loggamma(L+1+1j*m) for m in ms)

def triples(M,nu,w,u,v):
    m=[M+w/2+u,M-w/2+v,M+w/2-u,M-w/2-v]
    return [(m[i],m[(i+1)%4],nu) for i in range(4)],m

def run():
    checks=[]; gamma=[]; shift=[]; boxes=[]; lattice=[]; observer=[]
    def check(name,error,tol):
        checks.append(dict(name=name,error=str(error),tolerance=str(tol),passed=bool(abs(error)<=tol)))
    a,b,c=sp.symbols('a b c',real=True); ms=[a,b,c]; S=sum(ms)
    z=[(1+sp.I*S)/2]+[(1-sp.I*(S-2*m))/2 for m in ms]
    leg=[1+s*sp.I*m for m in ms for s in [-1,1]]
    for j,want in enumerate([-2,2,4]):
        err=sp.simplify(sum(x**j for x in leg)-2*sum(x**j for x in z)-want)
        check(f'balanced-moment-{j}',int(err),0)
    for sample in [(4,5,6),(20,21,4),(100,102,20),(1000,1001,40)]:
        ms=list(map(mp.mpf,sample))
        err=abs(mp.exp(mp.re(logF(ms))-logT(ms))-1)
        check(f'Gamma-hyperbolic-{sample}',err,mp.mpf('1e-55'))
        e=logE(ms); bound=ebound(ms)
        check(f'vertex-exponential-bound-{sample}',abs(e),bound)
        gamma.append(dict(masses=list(map(str,ms)),relative_identity_error=str(err),log_remainder=str(e),bound=str(bound)))
    ref=list(map(mp.mpf,[4,5,6]))
    for k in [1000,10000,100000]:
        for sample in [(6,7,8),(20,21,4)]:
            ms=list(map(mp.mpf,sample)); v=mp.exp(logB(ms,k)-logB(ref,k))
            with mp.workdps(95): high=mp.exp(logB(ms,k)-logB(ref,k))
            check(f'large-Gamma-precision-{k}-{sample}',abs(v/high-1),mp.mpf('1e-52'))
            shift.append(dict(kappa=k,masses=list(map(str,ms)),ratio=str(v),relative_to_one=str(abs(v-1))))
    # Sum over n analytically using polygamma; m tail has an explicit bound.
    for k in [100,1000]:
        for z in [mp.mpc('.5','3'),mp.mpc('1','-8')]:
            N=40
            partial=sum(mp.polygamma(2,z-1j*k*m) for m in range(1,N+1))
            tail=16/(mp.mpf(k)**3)*mp.zeta(3,N+1)+8/(mp.mpf(k)**2)*mp.zeta(2,N+1)
            totalbound=16*mp.zeta(3)/k**3+8*mp.zeta(2)/k**2
            check(f'lattice-derivative-bound-{k}-{z}',abs(partial)+tail,totalbound)
            lattice.append(dict(kappa=k,z=str(z),m_cutoff=N,partial=str(partial),tail_bound=str(tail),total_bound=str(totalbound)))
    for nu0 in [100,1000,10000]:
        nu=mp.mpf(nu0); M=nu**mp.mpf('1.6'); A=nu**2/(2*M)
        for y in [0,1,2]:
            w=1j*A*y;u=mp.sqrt(M)/2;v=-mp.sqrt(M)/3
            ts,m=triples(M,nu,w,u,v); rs,_=triples(M,nu,0,0,0)
            poly=mp.fprod(m)/M**4
            exact=mp.exp(sum(logT(t) for t in ts)-sum(logT(t) for t in rs))
            correction=sum(logE(t) for t in ts)-sum(logE(t) for t in rs)
            bound=sum(ebound(t) for t in ts+rs)
            check(f'graph-factorization-{nu0}-{y}',abs(exact/(poly*mp.exp(correction))-1),mp.mpf('1e-52'))
            check(f'graph-exponential-bound-{nu0}-{y}',abs(correction),bound)
            boxes.append(dict(nu=str(nu),M=str(M),y=y,u=str(u),v=str(v),ratio=str(exact),polynomial=str(poly),log_correction=str(correction),log_bound=str(bound),relative_to_one=str(abs(exact-1))))
    for k in [10000,100000,1000000]:
        M=mp.mpf(100); ts,m=triples(M,20,mp.mpf('.4'),mp.mpf('2'),mp.mpf('-3'))
        # Log-sinh avoids the enormous common entropy factor.
        lz=lambda x: 2*mp.log(mp.sin(mp.pi*x/k))+2*mp.log(mp.sinh(mp.pi*(k-x)))
        ratio=mp.exp(sum(lz(x)-lz(M) for x in m)); target=(mp.fprod(m)/M**4)**2
        with mp.workdps(95): high=mp.exp(sum(lz(x)-lz(M) for x in m))
        check(f'observer-precision-{k}',abs(ratio/high-1),mp.mpf('1e-50'))
        observer.append(dict(kappa=k,masses=list(map(str,m)),ratio=str(ratio),target=str(target),relative_error=str(abs(ratio/target-1))))
    inputs=['code/recoil_vertex.py','sources/chandra-liska-post-2610.30247v1.html','sources/dlmf-5.11.html']
    out=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='passed' if all(c['passed'] for c in checks) else 'failed',
             scope='exact cancellation moments, limiting vertex and local box checks; not full double-Gamma/6j numerical evaluation or finite-G OTOC',
             environment=dict(python=sys.version,platform=platform.platform(),mpmath=mp.__version__,sympy=sp.__version__,dps=65,crosscheck_dps=95),
             checks=checks,gamma_rows=gamma,shift_rows=shift,lattice_rows=lattice,box_rows=boxes,observer_rows=observer,
             hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    dest=ROOT/'results/stage4/recoil-vertex';dest.mkdir(parents=True,exist_ok=True)
    (dest/'diagnostics.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=out['status'],checks=len(checks),failed=sum(not c['passed'] for c in checks))))
    if out['status']!='passed':raise SystemExit(1)

if __name__=='__main__':run()
