"""Exact angle-Gram geometry and conditional one-loop diagnostics.

No full Virasoro amplitude or physical external energy contour is computed.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib,json,sys,platform
import mpmath as mp
import sympy as sp
from wkb_overlap import scaled_kernel

ROOT=Path(__file__).resolve().parents[1]
mp.mp.dps=65
PAIRS=[(0,1,4),(0,2,5),(0,3,6),(1,2,3),(1,3,2),(2,3,1)]
SETS=[(1,4),(2,5),(3,6),(1,2,6),(1,3,5),(2,3,4),(4,5,6),(1,2,3,4,5,6)]

def determinant(theta):
    c=[0]+[mp.cos(x) for x in theta];G=mp.eye(4)
    for i,j,k in PAIRS: G[i,j]=G[j,i]=-c[k]
    return mp.det(G)

def discriminant(theta):
    x=[0]+[mp.exp(-1j*t) for t in theta]
    q2=sum(mp.fprod(x[i] for i in a) for a in SETS)
    q0=sum(mp.fprod(1/x[i] for i in a) for a in SETS)
    q1=-sum((1/x[i]-x[i])*(1/x[i+3]-x[i+3]) for i in [1,2,3])
    return (q1*q1-4*q0*q2)/16

def angles(k,H,nu,nup,w):
    d=mp.pi*w/(2*k)
    return [H+d,mp.pi-mp.pi*nup/k,H-d,H+d,mp.pi-mp.pi*nu/k,H-d]

def factorized(k,H,nu,nup,w):
    d=mp.pi*w/(2*k); b=mp.pi*nu/k;bp=mp.pi*nup/k
    return 16*(mp.cos(b/2)**2*mp.cos(bp/2)**2-mp.cos(H)**2*mp.cos(d)**2)*(mp.sin(b/2)**2*mp.sin(bp/2)**2-mp.sin(H)**2*mp.sin(d)**2)

def weighted_I(A,t):
    A=mp.mpf(A);t=mp.mpf(t);c=mp.cosh(t);s=mp.sinh(t)
    def f(u):
        z=s+mp.sqrt(c/A)*u
        e=2*A*(-z*mp.asinh(z)+mp.sqrt(1+z*z)-1+z*t-(c-1))
        return mp.exp(e)/mp.sqrt(1+z*z)
    return mp.sqrt(A*c)*mp.quad(f,[-mp.inf,-8,-3,0,3,8,mp.inf])

def area(a,b,c):
    return mp.sqrt((a+b+c)*(a+b-c)*(a+c-b)*(b+c-a))/4

def recoil_action(M,w,u,v,nu):
    """CSTY 2.54 and B.6; branch continued from nondegenerate real w=u=v=0."""
    ms=[M+w/2+u,M-w/2+v,M+w/2-u,M-w/2-v]
    m=mp.matrix([[0,ms[1],nu,ms[0]],[ms[1],0,ms[2],nu],[nu,ms[2],0,ms[3]],[ms[0],nu,ms[3],0]])
    action=0
    for i in range(4):
        for j in range(i+1,4):
            k,l=[x for x in range(4) if x not in [i,j]]
            a,b,c,d,e,f=m[i,j],m[i,k],m[i,l],m[j,k],m[j,l],m[k,l]
            co=((a*a+b*b-c*c)*(a*a+d*d-e*e)-2*a*a*(b*b+d*d-f*f))/(16*area(a,b,c)*area(a,d,e))
            action+=a*mp.acos(co)
    return action

def transverse_h(M,w,nu):
    return 4*mp.sqrt((nu**4-4*M*M*w*w)*(4*M*M-2*nu*nu+w*w))/((4*M*M-nu*nu)*(nu*nu-w*w))

def run():
    checks=[]; geometry=[]; fixed=[]; recoil=[]; integrals=[]; transverse=[]
    def check(name,e,tol): checks.append(dict(name=name,error=str(e),tolerance=str(tol),passed=bool(e<=tol)))
    a,b,u,v=sp.symbols('a b u v')
    G=sp.Matrix([[1,-a,-v,-b],[-a,1,-b,-u],[-v,-b,1,-a],[-b,-u,-a,1]])
    fact=((1-u)*(1-v)-(a+b)**2)*((1+u)*(1+v)-(a-b)**2)
    check('symbolic-determinant-factorization',int(sp.expand(G.det()-fact)!=0),0)
    J=sp.Matrix([[1,sp.Rational(1,2),1,0],[1,-sp.Rational(1,2),0,1],[1,sp.Rational(1,2),-1,0],[1,-sp.Rational(1,2),0,-1]])
    check('mass-jacobian',abs(int(J.det()))-4,0)
    M,w,U,V,t=sp.symbols('M w U V t',real=True)
    masses=J*sp.Matrix([M,w,U,V]); gaps=[sp.pi/2-sp.I*t,sp.pi/2+sp.I*t]*2
    phase=sp.expand(sum(m*g for m,g in zip(masses,gaps))-2*sp.pi*M+2*sp.I*w*t)
    check('antipodal-phase',int(phase!=0),0)
    for idx in range(5):
        th=[mp.mpf('1.1')+mp.mpf('.1')*j+mp.mpf('.013')*idx for j in range(6)]
        if idx==4: th[2]+=mp.mpc(0,'.07')
        d=determinant(th);disc=discriminant(th)
        check(f'discriminant-{idx}',abs(disc/d-1),mp.mpf('1e-55'))
        geometry.append(dict(theta=[str(x) for x in th],det=str(d),discriminant_over_16=str(disc)))
    for k in [1000,10000,100000]:
        k=mp.mpf(k);H=2*mp.pi/3;nu=k**mp.mpf('.65');nup=mp.mpf('.8')*nu
        A=mp.pi*nu*nup/(2*k*mp.sin(H));d0=factorized(k,H,nu,nup,0)
        for z in [mp.mpf('.4'),mp.mpc(0,'.7'),mp.mpc(0,'1.5')]:
            om=A*z;f=factorized(k,H,nu,nup,om);d=determinant(angles(k,H,nu,nup,om))
            check(f'exact-factor-{k}-{z}',abs(d/f-1),mp.mpf('1e-50'))
            with mp.workdps(95): hi=factorized(k,H,nu,nup,om)/factorized(k,H,nu,nup,0)
            ratio=f/d0
            check(f'fixed-angle-precision-{k}-{z}',abs(ratio/hi-1),mp.mpf('1e-50'))
            fixed.append(dict(kappa=str(k),thetaH=str(H),nu=str(nu),nup=str(nup),A=str(A),u=str(z),det_ratio=str(ratio),target=str(1-z*z),relative_error=str(abs(ratio/(1-z*z)-1)),small_parameter=str((nu*nu+nup*nup)/(k*k))))
    m=mp.mpf(100);nu=mp.mpf(20);nup=mp.mpf(30);A0=nu*nup/(2*m);B=4*m*m-nu*nu-nup*nup
    for k in [10000,100000,1000000]:
        k=mp.mpf(k);H=mp.pi*(1-m/k)
        for z in [mp.mpf('.4'),mp.mpc(0,1)]:
            om=A0*z
            ratio=factorized(k,H,nu,nup,om)/factorized(k,H,nu,nup,0)
            target=(1+om*om/B)*(1-om*om/(A0*A0))
            with mp.workdps(95): hi=factorized(k,H,nu,nup,om)/factorized(k,H,nu,nup,0)
            check(f'recoil-precision-{k}-{z}',abs(ratio/hi-1),mp.mpf('1e-47'))
            recoil.append(dict(kappa=str(k),m=str(m),nu=str(nu),nup=str(nup),u=str(z),ratio=str(ratio),exact_recoil_limit=str(target),relative_error=str(abs(ratio/target-1)),simplified_target=str(1-z*z)))
    for A in [20,80,320]:
        for t in [mp.mpf('.5'),mp.mpf(1)]:
            ratio=(weighted_I(A,t)/weighted_I(A,0))/(scaled_kernel(A*mp.cosh(t))/scaled_kernel(A))
            with mp.workdps(90): hi=(weighted_I(A,t)/weighted_I(A,0))/(scaled_kernel(A*mp.cosh(t))/scaled_kernel(A))
            check(f'gram-inverse-precision-{A}-{t}',abs(ratio/hi-1),mp.mpf('1e-45'))
            integrals.append(dict(A=A,t=str(t),normalized_gram_I_over_H=str(ratio),target_cosh=str(mp.cosh(t)),relative_error=str(ratio/mp.cosh(t)-1)))
    for M in [100,1000,10000]:
        M=mp.mpf(M);nu=mp.mpf(20);A0=nu*nu/(2*M);B=4*M*M-2*nu*nu
        for z in [mp.mpf('.4'),mp.mpc(0,1),mp.mpc(0,2)]:
            om=A0*z;fn=lambda u,v:recoil_action(M,om,u,v,nu)
            exact=transverse_h(M,om,nu)
            h11=mp.diff(fn,(0,0),(2,0));h22=mp.diff(fn,(0,0),(0,2));h12=mp.diff(fn,(0,0),(1,1))
            check(f'transverse-h11-{M}-{z}',abs(h11/exact-1),mp.mpf('1e-48'))
            check(f'transverse-h22-{M}-{z}',abs(h22/exact-1),mp.mpf('1e-48'))
            check(f'transverse-h12-{M}-{z}',abs(h12/exact),mp.mpf('1e-48'))
            check(f'transverse-stationary-{M}-{z}',max(abs(mp.diff(fn,(0,0),(1,0))),abs(mp.diff(fn,(0,0),(0,1)))),mp.mpf('1e-48'))
            rg=1/mp.sqrt((1+om*om/B)*(1-om*om/(A0*A0)))
            rh=transverse_h(M,0,nu)/exact
            combined=(1-om*om/(nu*nu))/((1+om*om/B)*(1-om*om/(A0*A0)))
            check(f'combined-factor-{M}-{z}',abs(rg*rh/combined-1),mp.mpf('1e-52'))
            transverse.append(dict(M=str(M),nu=str(nu),A0=str(A0),u=str(z),h11=str(h11),h22=str(h22),h12=str(h12),analytic_h=str(exact),gram_ratio=str(rg),transverse_ratio=str(rh),combined_ratio=str(combined),target=str(1/(1-z*z)),relative_error=str(abs(combined*(1-z*z)-1))))
    files=['code/gram_prefactor.py','code/wkb_overlap.py','sources/chandra-liska-post-2610.30247v1.html','sources/chen-stanford-tang-yang-2607.14042v1.html']
    out=ROOT/'results/stage4/gram-prefactor';out.mkdir(parents=True,exist_ok=True)
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='passed' if all(x['passed'] for x in checks) else 'failed',scope='exact Gram algebra, coordinate pushforward, classical transverse Hessian and conditional Gram-weighted action integral; not full 6j asymptotics or finite-G certification',environment=dict(python=sys.version,platform=platform.platform(),mpmath=mp.__version__,sympy=sp.__version__,dps=65,crosscheck_dps=[90,95]),checks=checks,general_geometry=geometry,fixed_angle=fixed,pure_recoil=recoil,weighted_integrals=integrals,transverse_hessian=transverse,hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files})
    (out/'diagnostics.json').write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=data['status'],checks=len(checks),failed=sum(not x['passed'] for x in checks))))

if __name__=='__main__':run()
