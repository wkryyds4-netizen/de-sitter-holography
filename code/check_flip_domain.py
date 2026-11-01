"""Exact Weyl/domain algebra plus independent high-precision pole controls."""
from pathlib import Path
import json,hashlib,platform
import sympy as s
import mpmath as m
root=Path(__file__).resolve().parents[1];base=root/'results/stage4/flip-domain'
p=base/'diagnostics.json'
if p.exists():raise RuntimeError('Do not overwrite diagnostics')
checks=[]
def ck(name,ok,**kw):checks.append(dict(name=name,passed=bool(ok),**kw))
h=s.symbols('h',positive=True,real=True);a,b=s.symbols('a b',positive=True,real=True)
J=s.zeros(4);J[:2,2:]=s.eye(2);J[2:,:2]=-s.eye(2)
zb=s.Matrix([[1,-s.I*h,0,0]]);w=s.Matrix([[0,0,s.I*h,1]])
ck('formal cross-chiral commutator zero',s.simplify((zb*J*w.T)[0])==0)
ck('spectral real-imaginary pair nonzero',(s.Matrix([[1,0,0,0]])*J*s.Matrix([0,0,h,0]))[0]==h)
ck('bounded Weyl anticommutation phase',s.simplify(s.exp(-s.I*s.pi*h*(1/h)))==-1)
xp=s.pi*(a*a+a+b*b)/b;yp=-s.pi*(a+1)/b
ck('moving pole location',s.simplify(xp+(a+s.I*b)*yp+s.I*s.pi*(1+a+s.I*b))==0)
ainv=(-xp*yp-s.pi**2)/(s.pi**2+yp**2);binv=s.pi*(xp-yp)/(s.pi**2+yp**2)
ck('inverse pole map real',s.simplify(ainv-a)==0)
ck('inverse pole map imaginary',s.simplify(binv-b)==0)
ck('pole sweep wedge',s.simplify(xp*yp+s.pi**2)==-s.pi**2*a*((a+1)**2+b*b)/b**2)
m.mp.dps=50
def phi(z,hh):
    q=m.exp(m.pi*1j*hh);qt=m.exp(-m.pi*1j/hh)
    return m.qp(-qt*m.exp(z/hh),qt*qt)/m.qp(-q*m.exp(z),q*q)
def residue(hh):
    q=m.exp(m.pi*1j*hh);qt=m.exp(-m.pi*1j/hh)
    return -m.qp(qt*qt,qt*qt)/m.qp(q*q,q*q)
pole_rows=[]
for hh in [m.mpc(0,'.5'),m.j,m.mpc(0,2),m.mpc('.4','.7')]:
    z0=-m.pi*1j*(1+hh);R=residue(hh)
    for theta in [m.mpf('.2'),m.mpf('1.1')]:
        eps=m.mpf('1e-12')*m.exp(1j*theta)
        ratio=eps*phi(z0+eps,hh)/R
        err=abs(ratio-1)
        pole_rows.append(dict(parameter=[str(hh.real),str(hh.imag)],theta=str(theta),relative_error=m.nstr(err,25),residue=[m.nstr(R.real,30),m.nstr(R.imag,30)]))
        ck('first pole residue',err<m.mpf('1e-9'),parameter=str(hh),theta=str(theta))
shells=[]
for hh in [m.mpf('.5'),m.mpf(1),m.mpf(2)]:
    R=residue(1j*hh)
    C=2/(m.pi*hh)*abs(R)**2*m.exp(-m.pi*(hh*hh+1/(hh*hh)))
    row=dict(h=str(hh),residue=m.nstr(R.real,35),log_coefficient=m.nstr(C,35),shells=[])
    for radius in ['.2','.1','.05','.025']:
        r=m.mpf(radius);n=64;total=m.mpf(0)
        for k in range(n):
            theta=2*m.pi*k/n;x=m.pi*hh+r*m.cos(theta);y=-m.pi/hh+r*m.sin(theta)/hh
            z=m.pi*hh-1j*m.pi+r*m.exp(1j*theta)
            total+=abs(phi(z,1j*hh))**2*m.exp(-(x*x+y*y)/m.pi)/m.pi**2
        S=r*r/hh*2*m.pi/n*total
        row['shells'].append(dict(radius=radius,angular_nodes=n,shell=m.nstr(S,35),relative_to_limit=m.nstr(S/C-1,30)))
    ck('positive logarithmic coefficient',C>0,h=str(hh))
    ck('small-shell empirical check',abs(m.mpf(row['shells'][-1]['relative_to_limit']))<m.mpf('.01'),h=str(hh))
    shells.append(row)
ck('h=1 residue exact compact-product equality',residue(m.j)==-1)
data=dict(status='passed' if all(c['passed'] for c in checks) else 'failed',checks=checks,
    pole_rows=pole_rows,shells=shells,precision=50,
    scope='Exact canonical/pole-map identities and high-precision scalar controls. Pole/domain proof is analytic; shell quadrature is empirical. No pentagon disproof.',
    versions=dict(python=platform.python_version(),mpmath=m.__version__,sympy=s.__version__),
    hashes={'code/check_flip_domain.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
p.write_text(json.dumps(data,indent=2));print(data['status'],len(checks))
