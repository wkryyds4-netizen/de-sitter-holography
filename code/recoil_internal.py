"""Continued double-sine internal integral: exact shifts and pure-recoil Barnes limit."""
from pathlib import Path
from datetime import datetime, timezone
import sys,platform,json,hashlib
import mpmath as mp
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
VERT=[(0,1,2),(0,4,5),(1,3,5),(2,3,4)]
QUAD=[(0,1,3,4),(0,2,3,5),(1,2,4,5)]
mp.mp.dps=45

def parameters(ms):
    return ([sum(ms[i] for i in v)/2 for v in VERT],
            [sum(ms[i] for i in q)/2 for q in QUAD])

def barnes(z,a,r):
    return mp.exp(-1j*mp.pi*z+sum(mp.loggamma(z+1j*x) for x in a)
                  +sum(mp.loggamma(mp.mpf('.5')-z-1j*x) for x in r)
                  -mp.loggamma(z+mp.mpf('.5')))

def sb_scaled(w,k):
    """CLP A.5 at x=b*w, using q-Pochhammer products including constant."""
    q=mp.exp(-2*mp.pi/k);t=mp.exp(-2*mp.pi*k);b2=1j/k
    Q2=b2+2+1/b2
    return (mp.qp(t*mp.exp(2j*mp.pi*w),t,maxterms=100000)
            /mp.qp(q**w,q,maxterms=100000)
            *mp.exp(1j*mp.pi*w/2+mp.pi*(w*w-w)/(2*k)-1j*mp.pi*(Q2+1)/12))

def direct_normalized(z,a,r,k):
    L=-1j*k;q=mp.exp(-2*mp.pi/k)
    value=mp.fprod(sb_scaled(z+1j*x,k) for x in a)/(
        mp.fprod(sb_scaled(z+mp.mpf('.5')+1j*x+L,k) for x in r)
        *sb_scaled(z+mp.mpf('.5')-L,k))
    return -8*mp.pi**3*mp.exp(-mp.pi*k)*(1-q)**2*value

def normalized(z,a,r,k):
    q=mp.exp(-2*mp.pi/k);t=mp.exp(-2*mp.pi*k)
    num=[z+1j*x for x in a];den=[z+mp.mpf('.5')+1j*x for x in r]+[z+mp.mpf('.5')]
    delta=sum(w*w-w for w in num)-sum(w*w-w for w in den)
    g=mp.fprod(mp.qgamma(w,q,maxterms=100000) for w in num)/mp.fprod(
        mp.qgamma(w,q,maxterms=100000) for w in den)
    dual=mp.fprod(mp.qp(t*mp.exp(2j*mp.pi*w),t,maxterms=100000) for w in num)/mp.fprod(
        mp.qp(t*mp.exp(2j*mp.pi*w),t,maxterms=100000) for w in den)
    return (mp.pi**3*mp.exp(-1j*mp.pi*z+mp.pi*delta/(2*k))
            *(1+t*mp.exp(2j*mp.pi*z))*g*dual
            /mp.fprod(mp.sin(mp.pi*w) for w in den[:3]))

def vertical(a,r,dps=45):
    with mp.workdps(dps):
        shifts=sorted(set([mp.mpf(0)]+[-x for x in a]+[-x for x in r]))
        intervals=[-mp.inf]+shifts+[mp.inf]
        return mp.quad(lambda y:barnes(mp.mpf('.25')+1j*y,a,r),intervals)/(2*mp.pi)

def residue(a,r,dps=45):
    with mp.workdps(dps):
        terms=[]
        for i,x in enumerate(a):
            others=[a[j] for j in range(4) if j!=i]
            coeff=mp.exp(-mp.pi*x)*mp.fprod(mp.gamma(1j*(v-x)) for v in others)*mp.fprod(
                mp.gamma(mp.mpf('.5')+1j*(x-v)) for v in r)/mp.gamma(mp.mpf('.5')-1j*x)
            upper=[mp.mpf('.5')+1j*x]+[mp.mpf('.5')+1j*(x-v) for v in r]
            lower=[1+1j*(x-v) for v in others]
            terms.append(coeff*mp.hyper(upper,lower,1))
        return sum(terms),sum(abs(t) for t in terms),terms

def run():
    checks=[]; point=[]; integrals=[]; tails=[]
    def check(name,err,tol):
        checks.append(dict(name=name,error=str(err),tolerance=str(tol),passed=bool(abs(err)<=tol)))
    es=sp.symbols('e1:7',real=True);aa,rr=parameters(es)
    check('mass-sum-balance',int(sp.simplify(sum(aa)-sum(rr))),0)
    z=sp.symbols('z');num=[z+sp.I*x for x in aa];den=[z+sp.Rational(1,2)+sp.I*x for x in rr]+[z+sp.Rational(1,2)]
    check('small-argument-exponent',int(sp.simplify(sum(num)-sum(den)+2)),0)
    check('quadratic-z-coefficient',int(sp.simplify(sp.diff(sum(w*w-w for w in num)-sum(w*w-w for w in den),z)+4)),0)
    M,nu,w,u,v=sp.symbols('M nu w u v')
    cyc=[M+w/2+u,M-w/2+v,M+w/2-u,M-w/2-v]
    ag,rg=parameters([cyc[0],nu,cyc[1],cyc[2],nu,cyc[3]])
    check('symmetric-left-independent-omega',int(sp.simplify(sum(sp.diff(x,w)**2 for x in ag))),0)
    ms=list(map(mp.mpf,['2','1.1','2.3','2.1','1.2','2.4']));a,r=parameters(ms)
    c=(max(a)+min(r))/2
    for k in [20,60,180]:
        for x in ['-1.2','.25','1.7']:
            z=mp.mpc(x,-c)
            ex=normalized(z,a,r,k);dr=direct_normalized(z,a,r,k);lim=barnes(z,a,r)
            check(f'exact-shift-{k}-{x}',abs(ex/dr-1),mp.mpf('1e-34'))
            with mp.workdps(65):hi=normalized(z,a,r,k)
            check(f'point-precision-{k}-{x}',abs(ex/hi-1),mp.mpf('1e-32'))
            point.append(dict(kappa=k,z=str(z),exact=str(ex),limit=str(lim),relative_to_limit=str(abs(ex/lim-1))))
    for ms0 in [['2','1.1','2.3','2.1','1.2','2.4'],['2','1','2','2','1','2'],
                ['2','1.1','2.3','2.07','1.27','2.41']]:
        ms=list(map(mp.mpf,ms0));a,r=parameters(ms)
        val=vertical(a,r);hi=vertical(a,r,65)
        check('vertical-precision-'+','.join(ms0),abs(val/hi-1),mp.mpf('1e-32'))
        row=dict(masses=ms0,a=list(map(str,a)),r=list(map(str,r)),vertical=str(val),high=str(hi))
        if len(set(a))==4:
            res,absum,terms=residue(a,r)
            check('four-residue-hypergeometric',abs(res/hi-1),mp.mpf('1e-30'))
            row.update(residue=str(res),cancellation_ratio=str(absum/abs(res)),residue_terms=list(map(str,terms)))
        integrals.append(row)
    # Finite-kappa tails are tested on both sides; theorem supplies the uniform bound.
    ms=list(map(mp.mpf,['2','1.1','2.3','2.1','1.2','2.4']));a,r=parameters(ms);c=(max(a)+min(r))/2
    for k in [20,60]:
        for x in [-120,-40,40,120]:
            z=mp.mpc(x,-c);ex=normalized(z,a,r,k)
            with mp.workdps(65):hi=normalized(z,a,r,k)
            check(f'tail-precision-{k}-{x}',abs(ex/hi-1),mp.mpf('1e-30'))
            tails.append(dict(kappa=k,x=x,value=str(ex),scaled_x2=str(abs(ex)*x*x)))
    inputs=['code/recoil_internal.py','sources/chandra-liska-post-2610.30247v1.html','sources/dlmf-5.11.html']
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='passed' if all(c['passed'] for c in checks) else 'failed',
              scope='internal fixed-mass recoil integral only; exact local identities, limiting quadrature and tail samples, no external physical OTOC',
              environment=dict(python=sys.version,platform=platform.platform(),mpmath=mp.__version__,sympy=sp.__version__,dps=45,crosscheck_dps=65),
              checks=checks,point_rows=point,integral_rows=integrals,tail_rows=tails,
              hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    out=ROOT/'results/stage4/recoil-internal';out.mkdir(parents=True,exist_ok=True)
    (out/'diagnostics.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=data['status'],checks=len(checks),failed=sum(not c['passed'] for c in checks))))
    if data['status']!='passed':raise SystemExit(1)
if __name__=='__main__':run()
