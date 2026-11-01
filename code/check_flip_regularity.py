"""Local phase-vortex controls, independent of the floating-point FFT helper.

The non-H1 conclusion follows from the Laurent expansion, not these samples.
Finite-radius angular quadrature is empirical; this is not a pentagon test.
"""
from pathlib import Path
import json,hashlib,platform
import mpmath as m
import sympy as s

root=Path(__file__).resolve().parents[1]
out=root/'results/stage4/flip-domain/regularity.json'
if out.exists():raise RuntimeError('Do not overwrite regularity controls')
checks=[]
def ck(name,value,**kw):checks.append(dict(name=name,passed=bool(value),**kw))
h=s.symbols('h',positive=True);theta=s.symbols('theta',real=True)
angular=s.integrate(4*(s.sin(theta)**2+h*h*s.cos(theta)**2)/h,(theta,0,2*s.pi))
ck('exact elliptic angular coefficient',s.simplify(angular-4*s.pi*(1+h*h)/h)==0)
m.mp.dps=50
def phi(z,h):
    q=m.exp(-m.pi*h);t=m.exp(-m.pi/h)
    return m.qp(-t*m.exp(z/(1j*h)),t*t)/m.qp(-q*m.exp(z),q*q)
rows=[]
for hv in ['.5','1','2']:
    hh=m.mpf(hv);z0=m.pi*hh-m.pi*1j
    limit=4*(1+hh*hh)/(m.pi*hh)*m.exp(-m.pi*(hh*hh+hh**-2))
    for radius in ['.0125','.00625']:
        r=m.mpf(radius);values=[]
        for n in [32,64]:
            total=m.mpf(0);phases=[]
            for k in range(n):
                t=2*m.pi*k/n;z=z0+r*m.exp(1j*t)
                v=phi(z,hh);L=m.diff(lambda zz:phi(zz,hh),z)/v
                x=z.real;y=z.imag/hh
                rho=m.exp(-(x*x+y*y)/m.pi)/m.pi**2
                total+=rho*((x/m.pi)**2+(y/m.pi)**2+4*L.imag**2+4*hh*hh*L.real**2)
                phases.append(v/m.conj(v))
            value=r*r/hh*2*m.pi/n*total;values.append(value)
            winding=sum(m.arg(phases[(k+1)%n]/phases[k]) for k in range(n))/(2*m.pi)
            row=dict(h=hv,radius=radius,angular_nodes=n,shell=m.nstr(value,40),
                     limit=m.nstr(limit,40),relative_to_limit=m.nstr(value/limit-1,30),
                     winding=m.nstr(winding,30))
            rows.append(row)
            ck('small-radius gradient shell',abs(value/limit-1)<m.mpf('.01'),h=hv,radius=radius,n=n)
            ck('sampled phase winding',abs(winding+2)<m.mpf('1e-40'),h=hv,radius=radius,n=n)
        ck('angular doubling',abs(values[1]/values[0]-1)<m.mpf('1e-30'),h=hv,radius=radius)
data=dict(status='passed' if all(c['passed'] for c in checks) else 'failed',checks=checks,rows=rows,
          precision=50,versions=dict(python=platform.python_version(),mpmath=m.__version__,sympy=s.__version__),
          scope='Local Laurent/vortex controls only. Empirical angular quadrature, no continuum pentagon certification.',
          hashes={'code/check_flip_regularity.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
out.write_text(json.dumps(data,indent=2));print(data['status'],len(checks),flush=True)
