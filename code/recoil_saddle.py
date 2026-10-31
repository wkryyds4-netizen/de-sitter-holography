"""Reflected Barnes saddle: algebra and local diagnostics, not a full cycle proof."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, platform, sys
import mpmath as mp
import sympy as sp
from recoil_internal import parameters, barnes
from gram_prefactor import recoil_action, transverse_h

ROOT = Path(__file__).resolve().parents[1]
mp.mp.dps = 80

def labels(M, nu, w=0, u=0, v=0):
    return parameters([M+w/2+u, nu, M-w/2+v, M+w/2-u, nu, M-w/2-v])

def coefficients(a, r):
    e2 = lambda t: sum(t[i]*t[j] for i in range(len(t)) for j in range(i+1,len(t)))
    e3 = lambda t: sum(t[i]*t[j]*t[k] for i in range(len(t)) for j in range(i+1,len(t)) for k in range(j+1,len(t)))
    return e2(a)-e2(r), -e3(a)+e3(r), mp.fprod(a)

def saddle(a, r, sign=1):
    A,B,C = coefficients(a,r)
    delta = B*B-4*A*C
    # Transport the positive square root of -Delta; sqrt(Delta) lies on a cut.
    return -1j*(-B+sign*1j*mp.sqrt(-delta))/(2*A), delta

def phase(z,a,r):
    f = lambda t: t*(mp.log(t)-1)
    return 2j*mp.pi*z-mp.pi*sum(r)+sum(f(z+1j*t) for t in a)-sum(f(z+1j*t) for t in r)-f(z)

def phase_prime(z,a,r):
    return 2j*mp.pi+sum(mp.log(z+1j*t) for t in a)-sum(mp.log(z+1j*t) for t in r)-mp.log(z)

def phase_second(z,a,r):
    return sum(1/(z+1j*t) for t in a)-sum(1/(z+1j*t) for t in r)-1/z

def correction(z,r):
    return mp.fprod(1/(1+mp.exp(2j*mp.pi*(z+1j*t))) for t in r)

def reflected(z,a,r):
    lg = sum(mp.loggamma(z+1j*t) for t in a)-sum(mp.loggamma(z+mp.mpf('.5')+1j*t) for t in r)-mp.loggamma(z+mp.mpf('.5'))
    return 8*mp.pi**3*mp.exp(2j*mp.pi*z-mp.pi*sum(r)+lg)*correction(z,r)

def local_integrand(z,a,r):
    return 8*mp.pi**3*mp.exp(phase(z,a,r)-sum(mp.log(z+1j*t) for t in a)/2)

def phase_action(M,nu,w=0,u=0,v=0):
    a,r=labels(M,nu,w,u,v);z,_=saddle(a,r)
    return -2*mp.re(phase(z,a,r))

def run():
    checks=[]; rows=[]; scaling=[]; reflection_rows=[]
    def check(name,err,tol):
        checks.append(dict(name=name,error=str(abs(err)),tolerance=str(tol),passed=bool(abs(err)<=tol)))
    M,n,w,u,v,x=sp.symbols('M n w u v x')
    aa=[M+n/2+(p*u+q*v)/2 for p,q in [(1,1),(1,-1),(-1,1),(-1,-1)]]
    rr=[M+n+w/2,2*M,M+n-w/2]
    poly=sp.Poly(sp.prod(x-t for t in aa)-x*sp.prod(x-t for t in rr),x)
    check('quadratic-degree',poly.degree()-2,0)
    disc=sp.discriminant(poly.as_expr(),x)
    check('symmetric-discriminant',int(sp.expand(disc.subs({u:0,v:0})+(n**4-4*M*M*w*w)*(4*M*M-2*n*n+w*w)/16)!=0),0)
    ms=[M+w/2+u,M-w/2+v,M+w/2-u,M-w/2-v]
    faces=sp.Matrix([[0,ms[1],n,ms[0]],[ms[1],0,ms[2],n],[n,ms[2],0,ms[3]],[ms[0],n,ms[3],0]])
    cm=sp.ones(5);cm[0,0]=0
    for i in range(4):
        for j in range(4):
            if i==j: cm[i+1,j+1]=0
            else:
                k,l=[t for t in range(4) if t not in (i,j)]
                cm[i+1,j+1]=sp.expand(faces[k,l]**2)
    check('dual-Cayley-Menger',int(sp.expand(disc+cm.det(method='domain-ge')/32)!=0),0)
    A,B,C=poly.all_coeffs()
    c=-B/(2*A)
    check('symmetric-pole-height',int(sp.factor((M+n-c).subs({u:0,v:0})-n*(n*n+w*w)/(4*M*M+2*n*n+w*w))!=0),0)
    for ww,uu,vv in [('0','0','0'),('.7','1','2'),('1.2','-.7','1.3'),('1j','1','2'),('2j','.4','-.3')]:
        W=mp.mpc(0,ww[:-1]) if ww.endswith('j') else mp.mpf(ww)
        U,V=mp.mpf(uu),mp.mpf(vv);MM,NN=mp.mpf(100),mp.mpf(20)
        a,r=labels(MM,NN,W,U,V);z,D=saddle(a,r)
        check('stationary-'+ww,abs(phase_prime(z,a,r)),mp.mpf('1e-50'))
        N=mp.fprod(z+1j*t for t in a)
        check('determinant-'+ww,phase_second(z,a,r)*N+mp.sqrt(-D),mp.mpf('1e-48'))
        I=recoil_action(MM,W,U,V,NN);Ip=phase_action(MM,NN,W,U,V)
        check('classical-action-'+ww,Ip-I,mp.mpf('1e-48'))
        exact=reflected(z,a,r);direct=barnes(z,a,r)
        check('exact-reflection-'+ww,exact/direct-1,mp.mpf('1e-48'))
        with mp.workdps(110):
            ah,rh=labels(MM,NN,W,U,V);zh,Dh=saddle(ah,rh)
            hi=reflected(zh,ah,rh)
        check('precision-'+ww,exact/hi-1,mp.mpf('1e-46'))
        rows.append(dict(w=ww,u=uu,v=vv,z=str(z),discriminant=str(D),action=str(I),phase_action=str(Ip),
                         exact_reflection_factor=str(correction(z,r)),minimum_pole_height=str(min(mp.im(z+1j*t) for t in r))))
    for ww in [mp.mpf(0),mp.mpf('.7'),1j]:
        fn=lambda uu,vv:phase_action(mp.mpf(100),mp.mpf(20),ww,uu,vv)
        h=transverse_h(100,ww,20)
        for ij in [(2,0),(0,2),(1,1)]:
            got=mp.diff(fn,(mp.mpf(0),mp.mpf(0)),ij)
            check('transverse-'+str(ww)+'-'+str(ij),got-(0 if ij==(1,1) else h),mp.mpf('1e-45'))
    # Three distinct asymptotic sectors; small reflection factors are not assumed.
    for exponent in [mp.mpf('1.4'),mp.mpf('1.5'),mp.mpf('1.7')]:
        for NN in [mp.mpf('1e3'),mp.mpf('1e5'),mp.mpf('1e7')]:
            MM=NN**exponent;w0=1j*NN*NN/(2*MM)
            a,r=labels(MM,NN,w0);z,D=saddle(a,r)
            delta=mp.im(z+1j*r[0]);qmax=mp.exp(-2*mp.pi*delta)
            E=correction(z,r);bound=(1-qmax)**(-3)-1
            check('reflection-bound-'+str(exponent)+'-'+str(NN),max(mp.mpf(0),abs(E-1)-bound),mp.mpf('1e-48'))
            ratio=reflected(z,a,r)/local_integrand(z,a,r)
            # Remove exact pole factor to isolate Stirling, rather than conflate errors.
            stirling_ratio=ratio/E
            with mp.workdps(110):
                ah,rh=labels(MM,NN,w0);zh,_=saddle(ah,rh)
                eh=correction(zh,rh)
                rhigh=reflected(zh,ah,rh)/local_integrand(zh,ah,rh)
            check('large-mass-precision-'+str(exponent)+'-'+str(NN),ratio/rhigh-1,mp.mpf('1e-40'))
            scaling.append(dict(nu=str(NN),M=str(MM),exponent=str(exponent),w=str(w0),z=str(z),
                                pole_height=str(delta),qmax=str(qmax),reflection_factor=str(E),
                                reflection_error=str(abs(E-1)),bound=str(bound),stirling_relative_error=str(abs(stirling_ratio-1)),
                                full_local_relative_error=str(abs(ratio-1))))
    # Explicitly record the failed principal-log shortcut, not as a valid asymptotic.
    a,r=labels(mp.mpf(100),mp.mpf(20),mp.mpf('.7'),mp.mpf(1),mp.mpf(2));z,_=saddle(a,r)
    naive=-1j*mp.pi+sum(mp.log(z+1j*t) for t in a)-sum(mp.log(-z-1j*t) for t in r)-mp.log(z)
    check('principal-log-failure-is-minus-two-pi-i',naive+2j*mp.pi,mp.mpf('1e-48'))
    inputs=['code/recoil_saddle.py','code/recoil_internal.py','code/gram_prefactor.py',
            'sources/chen-stanford-tang-yang-2607.14042v1.html','sources/chandra-liska-post-2610.30247v1.html']
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='passed' if all(t['passed'] for t in checks) else 'failed',
              scope='symbolic saddle/volume identities and finite-point local asymptotics; no complete Barnes contour asymptotic',
              environment=dict(python=sys.version,platform=platform.platform(),mpmath=mp.__version__,sympy=sp.__version__,dps=80,crosscheck_dps=110),
              checks=checks,saddle_rows=rows,scaling_rows=scaling,naive_log_derivative=str(naive),
              hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    out=ROOT/'results/stage4/recoil-saddle';out.mkdir(parents=True,exist_ok=True)
    (out/'diagnostics.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=data['status'],checks=len(checks),failed=[t['name'] for t in checks if not t['passed']])))
    if data['status']!='passed':raise SystemExit(1)

if __name__=='__main__':run()
