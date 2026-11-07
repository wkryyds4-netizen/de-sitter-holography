"""Exact differential check of the standard Jacobi--Laguerre tensor map."""
from pathlib import Path
import json,hashlib
import sympy as s
ROOT=Path(__file__).resolve().parent.parent
target=ROOT/'results/susskind-two-matter/tensor-check.json'
if target.exists():raise FileExistsError(target)
P,x=s.symbols('P x',positive=True)
a,b=s.symbols('a b',positive=True)
g=s.Function('g')(P,x)
def dl(f):return s.diff(f,P)+(1-x)/P*s.diff(f,x)
def dr(f):return s.diff(f,P)-x/P*s.diff(f,x)
f=g/s.sqrt(P)
direct=s.sqrt(P)*(-P*x*dl(dl(f))-dl(f)-P*(1-x)*dr(dr(f))-dr(f)+a*a/(P*x)*f+b*b/(P*(1-x))*f)
relative=-s.diff(x*(1-x)*s.diff(g,x),x)+(a*a/x+b*b/(1-x)+s.Rational(1,4))*g
expected=-P*s.diff(g,P,2)-s.diff(g,P)+relative/P
assert s.simplify(direct-expected)==0
# Conjugate the relative Sturm--Liouville expression by its Jacobi weight.
y=s.Function('y')(x);w=x**a*(1-x)**b
T=-s.diff(x*(1-x)*s.diff(w*y,x),x)+(a*a/x+b*b/(1-x)+s.Rational(1,4))*w*y
L=-x*(1-x)*s.diff(y,x,2)-(2*a+1-(2*a+2*b+2)*x)*s.diff(y,x)
assert s.simplify(s.expand(T/w-L-(a+b+s.Rational(1,2))**2*y))==0
boost=s.sqrt(P)*s.I*(P*x*dl(f)+P*(1-x)*dr(f)+f)
assert s.simplify(boost-s.I*(P*s.diff(g,P)+g/2))==0
data=dict(status='passed',checks=['full_two_variable_S_transformation','relative_Jacobi_conjugation_general_parameters','total_boost_transformation'],sympy=s.__version__,input_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
target.write_text(json.dumps(data,indent=2)+'\n')
print(data['checks'])
