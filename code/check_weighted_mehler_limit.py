"""Fresh modular, expansion, and Fourier checks; no previous driver reruns."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,platform
import mpmath as mp
import sympy as sp
from weighted_mehler_limit import (product_density,modular_density,leading_density,
    leading_gap,coefficient,image_quotient,image_tail,symmetric_weight,gap,response)

ROOT=Path(__file__).resolve().parents[1]
DIR=ROOT/'results/susskind-weighted-mehler-limit'
OUT=DIR/'modular-checks.json'
if OUT.exists():raise FileExistsError(OUT)
DIR.mkdir(exist_ok=True)
mp.mp.dps=70
rows=[]
for d in (1,2,3):
    for h in ('0.4','0.1','0.04'):
        for y,z in [('0','0'),('.7','.4'),('-1.2','-1.5')]:
            v=product_density(h,y,z,d);w=modular_density(h,y,z,d)
            error=abs(v-w)/abs(v)
            assert error<mp.mpf('1e-55'),(d,h,y,z,error)
            rows.append(dict(d=d,h=h,y=y,z=z,relative_error=str(error),value=str(v)))
# Boundary-layer points: both a spectral endpoint and the diagonal are close.
edges=[]
for h in (mp.mpf('.1'),mp.mpf('.4')):
    for theta,phi in [(mp.mpf('1e-5'),mp.mpf('3e-5')),
                       (mp.pi-mp.mpf('1e-5'),mp.pi-mp.mpf('2e-5')),
                       (mp.mpf('.1'),mp.pi-mp.mpf('.1'))]:
        y=(theta+phi-mp.pi)/mp.sqrt(h);z=(theta-phi)/h
        v=product_density(h,y,z);w=modular_density(h,y,z)
        assert abs(v-w)/abs(v)<mp.mpf('1e-50')
        edges.append(dict(h=str(h),theta=str(theta),phi=str(phi),relative_error=str(abs(v-w)/abs(v)),
                          ratio_to_majorant_shape=str(w/leading_density(y,z,2))))

expansions=[]
for d in (1,2,3):
    for y,z in [(mp.mpf('.7'),mp.mpf('.4')),(mp.mpf('1.2'),mp.mpf('1.5'))]:
        for h in (mp.mpf('.02'),mp.mpf('.01'),mp.mpf('.005')):
            ratio=modular_density(h,y,z,d)/leading_density(y,z,d)
            remainder=(ratio-1-h*coefficient(y,z,d))/h**2
            assert abs(remainder)<100
            expansions.append(dict(d=d,y=str(y),z=str(z),h=str(h),density_second_remainder=str(remainder),
                weight_second_remainder=str((symmetric_weight(h,y,z)-(y*y-1)-h*(z*z+y*y/2-y**4/12))/h**2),
                gap_second_remainder=str((gap(h,y,z)+2*z*(1+h*(mp.mpf(1)/4-y*y/8)))/h**2)))

# Symbolic Gaussian integration of the complete first-order weighted integrand.
y,z,u,d=sp.symbols('y z u d',real=True)
A=(4*d+1)/8-(2*d+1)*y*y/8-z*z/2
W0=y*y-1;W1=z*z+y*y/2-y**4/12
phase=-2*sp.I*u*z*(sp.Rational(1,4)-y*y/8)
poly=sp.Poly(sp.expand(W1+W0*(A+phase-sp.Rational(1,2))),y)
integrated=0
for (n,),c in poly.terms():
    integrated+=c*(sp.factorial2(n-1) if n%2==0 else 0)
assert sp.simplify(integrated-(z*z-d/2+sp.I*u*z/2))==0
g=sp.sech(u)**(2*d)
assert sp.simplify((-sp.diff(g,u,2)/4-d*g/2-u*sp.diff(g,u)/4
    -d*g*(u*sp.tanh(u)-(1+2*d)*sp.tanh(u)**2)/2).rewrite(sp.cosh))==0

# Independent real Fourier integral, with explicit exponential-tail estimate.
fourier=[]
Z=mp.mpf(32)
for dd in (1,2,3):
    for uu in (mp.mpf(0),mp.mpf('.5'),mp.mpf(1),mp.mpf(2)):
        def integrand(zz):
            return 2*leading_gap(zz,dd)*((zz*zz-mp.mpf(dd)/2)*mp.cos(2*uu*zz)
                                           +uu*zz/2*mp.sin(2*uu*zz))
        values=[]
        for precision in (50,70):
            with mp.workdps(precision):values.append(mp.quad(integrand,[0,1,2,4,8,16,Z]))
        # On z>=Z>=1, prod(j²+z²)<=prod(1+j²) z^(2d-2).
        pref=mp.mpf(2)**(2*dd-1)/mp.factorial(2*dd-1)*2*mp.fprod(1+j*j for j in range(1,dd))/(1-mp.exp(-2*mp.pi*Z))
        tail=2*pref*(mp.gammainc(2*dd+2,mp.pi*Z,mp.inf)/mp.pi**(2*dd+2)
            +mp.mpf(dd)/2*mp.gammainc(2*dd,mp.pi*Z,mp.inf)/mp.pi**(2*dd)
            +abs(uu)/2*mp.gammainc(2*dd+1,mp.pi*Z,mp.inf)/mp.pi**(2*dd+1))
        err=abs(values[-1]-response(uu,dd))
        assert err<tail+mp.mpf('1e-45')
        fourier.append(dict(d=dd,u=str(uu),value=str(values[-1]),formula=str(response(uu,dd)),
            absolute_difference=str(err),precision_difference=str(abs(values[-1]-values[0])),
            analytic_tail_bound=str(tail)))
result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),
    product_modular=rows,boundary_layers=edges,local_expansions=expansions,fourier=fourier,
    symbolic_gaussian_integral=str(sp.simplify(integrated)),
    image_quotient_terms=8,image_tail_bound_at_h04=str(image_tail(mp.mpf('.4'),8)),
    environment=dict(python=platform.python_version(),mpmath=mp.__version__,sympy=sp.__version__),
    input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in [
        'code/weighted_mehler_limit.py','code/check_weighted_mehler_limit.py',
        'results/susskind-selected-rindler-response/response-checks.json']},
    scope='Exact symbolic Gaussian coefficient and independent high-precision representations. Finite sample checks complement the analytic global majorant and remainder proof; quadrature is not interval-certified.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',product_rows=len(rows),boundary_rows=len(edges),expansion_rows=len(expansions),fourier_rows=len(fourier),sha256=hashlib.sha256(OUT.read_bytes()).hexdigest())))
