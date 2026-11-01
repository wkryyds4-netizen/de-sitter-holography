"""Exact local geometry and numerical diagnostics of weak continuation."""
from pathlib import Path
import json,hashlib
import sympy as s
import mpmath as m
root=Path(__file__).resolve().parents[1];base=root/'results/stage4/weak-continuation'
out=base/'diagnostics.json'
if out.exists():raise RuntimeError('Do not overwrite')
checks=[]
def ck(name,ok,**kw):checks.append(dict(name=name,passed=bool(ok),**kw))
a,b=s.symbols('a b',positive=True,real=True);H=a+s.I*b
x=s.pi*(a*a+a+b*b)/b;y=-s.pi*(a+1)/b;w=y+s.I*s.pi
J=s.det(s.Matrix([[s.diff(x,a),s.diff(x,b)],[s.diff(y,a),s.diff(y,b)]]))
ck('pole-map Jacobian',s.simplify(J-s.pi**2*((a+1)**2+b*b)/b**3)==0)
ck('Jacobian equals absolute-w squared over b',s.simplify(J-w*s.conjugate(w)/b)==0)
ck('affine divisor factorization',s.simplify(x+s.I*s.pi+H*w)==0)
ck('inverse pole parameter',s.simplify(-(x+s.I*s.pi)/w-H)==0)
ck('density Jacobian factor',s.simplify(J/w-s.conjugate(w)/b)==0)
ck('w in parameter coordinates',s.simplify(w+s.pi*(s.conjugate(H)+1)/b)==0)
z,zb,r=s.symbols('z zb r',nonzero=True)
C=s.pi*r*r/s.Integer(9)*(1-(1-z*zb/r**2)**9)/z
ck('radial Cauchy dbar formula',s.simplify(s.diff(C,zb)-s.pi*(1-z*zb/r**2)**8)==0)
ck('meromorphic and literal difference',s.simplify(C-s.pi*r*r/(9*z)+s.pi*r*r/(9*z)*(1-z*zb/r**2)**9)==0)
ck('other-pole separation in Q',s.Rational(37,11)>s.Rational(7,3))
ck('chosen source and loop within Q',s.Rational(6,100)<s.Rational(1,10))
ck('source within loop',s.Rational(4,100)<s.Rational(6,100))
m.mp.dps=50
def dec(z):return m.mpc(*z)
study=json.loads((base/'study.json').read_text())
for row in study['rows']:
    if row['kind']=='integral':
        ck('independent local integral',m.mpf(row['direct_difference'])<m.mpf('1e-18'),offset=row['offset'])
        ck('smooth quadrature doubling',m.mpf(row['quadrature_difference'])<m.mpf('1e-25'),offset=row['offset'])
cr=[row for row in study['rows'] if row['kind']=='wirtinger']
for prev,nxt in zip(cr[:-1],cr[1:]):
    ratio=m.mpf(prev['relative_error'])/m.mpf(nxt['relative_error'])
    ck('second-order Wirtinger approach',m.mpf('3.8')<ratio<m.mpf('4.2'),delta=nxt['delta'],ratio=m.nstr(ratio,25))
ck('nonzero derivative approached',m.mpf(cr[-1]['relative_error'])<m.mpf('.001'))
for row in cr:
    if 'direct_difference' in row:ck('independent derivative quadrature',m.mpf(row['direct_difference'])<m.mpf('1e-25'),delta=row['delta'])
loops=[row for row in study['rows'] if row['kind']=='morera']
ck('nonzero Morera witness',m.mpf(loops[-1]['relative_error'])<m.mpf('1e-18'))
ck('Morera refinement improves',all(m.mpf(a['relative_error'])>m.mpf(b['relative_error']) for a,b in zip(loops,loops[1:])))
positive=json.loads((base/'positive-form.json').read_text())
ck('positive state normalized',abs(m.mpf(positive['norm'])-1)<m.mpf('1e-35'))
ck('positive-state derivative nonzero',m.mpf(positive['rows'][-1]['relative_error'])<m.mpf('.001'))
ck('positive-state refinement improves',all(m.mpf(a['relative_error'])>m.mpf(b['relative_error']) for a,b in zip(positive['rows'],positive['rows'][1:])))
data=dict(status='passed' if all(c['passed'] for c in checks) else 'failed',checks=checks,
    scope='Exact symbolic identities plus finite quadrature checks. Weak non-holomorphy is proved analytically, not certified by finite sampling.',
    versions=dict(sympy=s.__version__,mpmath=m.__version__),
    hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in ['code/check_weak_continuation.py','results/stage4/weak-continuation/study.json','results/stage4/weak-continuation/positive-form.json']})
out.write_text(json.dumps(data,indent=2));print(data['status'],len(checks))
