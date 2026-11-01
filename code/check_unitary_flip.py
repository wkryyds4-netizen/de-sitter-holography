from pathlib import Path
import json,hashlib
import sympy as s
import numpy as np
import mpmath as mp
from unitary_flip import phase
root=Path(__file__).resolve().parents[1];base=root/'results/stage4/unitary-flip'
out=base/'diagnostics.json'
if out.exists():raise RuntimeError('Do not overwrite diagnostics')
checks=[]
def ck(name,ok,**kw):checks.append(dict(name=name,passed=bool(ok),**kw))
U=s.Matrix([[1,0],[0,1],[-1,-1],[-1,-1],[0,1],[1,0]])
J=s.Matrix([[0,1],[-1,0]])
edges=[(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]
V=s.Matrix([[int(v in e) for v in range(1,5)] for e in edges])
E=s.zeros(6)
for a,b,c in [(1,3,2),(1,2,4),(1,4,3),(2,3,4)]:
    es=[edges.index(tuple(sorted(e))) for e in [(a,b),(b,c),(c,a)]]
    for i,j in zip(es,es[1:]+es[:1]):E[i,j]+=1;E[j,i]-=1
R=V*(V.T*V).inv()
for name,ok in [('oriented tetrahedron factorization',E==U*J*U.T),('exchange rank',E.rank()==2),
    ('puncture rank',V.rank()==4),('puncture kernel',E*V==s.zeros(6,4)),
    ('coordinate constraints',V.T*U==s.zeros(4,2)),('central section',V.T*R==s.eye(4)),
    ('CCR sign',-U*(U*J).T==E)]:ck(name,ok)
x,y=s.symbols('x y',real=True);f=s.Function('f')(x,y);M=s.exp(s.I*x*y/s.pi)
for coord,added in [(x,y),(y,x)]:
    lhs=s.simplify(s.conjugate(M)*(-s.pi*s.I)*s.diff(M*f,coord))
    rhs=-s.pi*s.I*s.diff(f,coord)+added*f
    ck('chirp conjugation '+str(coord),s.simplify(lhs-rhs)==0)
mp.mp.dps=55;scalar=[]
for h in [.25,.5,1.,2.]:
    for xx,yy in [(.2,.3),(2.,3.),(2.,-4.)]:
        z=mp.mpc(xx,h*yy);q=mp.exp(-mp.pi*h);qt=mp.exp(-mp.pi/h)
        ph=mp.qp(-qt*mp.exp(z/(1j*h)),qt*qt)/mp.qp(-q*mp.exp(z),q*q)
        val,tail=phase(np.array(xx),np.array(yy),h)
        error=float(abs(complex(ph/mp.conj(ph))-val))
        inv=float(abs(val*phase(np.array(-xx),np.array(-yy),h)[0]-np.exp(-1j*xx*yy/np.pi)))
        scalar.append(dict(h=h,x=xx,y=yy,error=error,involutivity_error=inv,**tail))
        ck('independent compact products',error<2e-13,h=h,x=xx,y=yy)
        ck('involutivity',inv<2e-13,h=h,x=xx,y=yy)
barnes=[];rot=mp.exp(-1j*mp.pi/4);eps=mp.mpf('.2')
for xx,yy in [('.2','.3'),('1','1'),('-1','.7')]:
    z=mp.mpc(xx,yy)
    def integrand(t):
        p=rot*(t+1j*eps)
        return mp.exp(-1j*p*z)*rot/(mp.sinh(mp.pi*p)*mp.sinh(mp.pi*1j*p)*p)
    logphi=-mp.quad(integrand,[-30,-4,-1,0,1,4,30])/4
    q=mp.exp(-mp.pi);ph=mp.qp(-q*mp.exp(z/1j),q*q)/mp.qp(-q*mp.exp(z),q*q)
    error=abs(mp.exp(logphi)/ph-1)
    barnes.append(dict(x=xx,y=yy,relative_error=mp.nstr(error,30)))
    ck('slanted Barnes vs compact',error<mp.mpf('1e-35'),x=xx,y=yy)
result=dict(status='passed' if all(r['passed'] for r in checks) else 'failed',checks=checks,
    matrices={k:str(v) for k,v in dict(U=U,J=J,V=V,E=E,R=R).items()},scalar=scalar,barnes=barnes,
    precision=55,barnes_contour='p=exp(-i*pi/4)*(t+0.2i), -30<t<30; empirical truncation, not interval certified',
    scope='Scalar and algebra checks only. Does not validate the pentagon reproduction.',
    hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in ['code/unitary_flip.py','code/check_unitary_flip.py']})
out.write_text(json.dumps(result,indent=2));print(result['status'],len(checks))
