"""Fresh exact trace coefficients and read-only analysis of saved packets."""
from pathlib import Path
from datetime import datetime,timezone
from fractions import Fraction
from math import comb,factorial,cosh
import hashlib,json,platform
import sympy as sp
from selected_rindler_response import energy_pair,add,scale,exponential_coefficient,response

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/susskind-weighted-mehler-limit/trace-mixture-checks.json'
if OUT.exists():raise FileExistsError(OUT)
u=sp.symbols('u')
g=sp.series(1/sp.cosh(u)**4,u,0,16).removeO().expand()
R=sp.series((u*sp.tanh(u)-5*sp.tanh(u)**2)/sp.cosh(u)**4,u,0,16).removeO().expand()
exact=[]
for n in range(0,16,2):
    poly={};k=n//2
    for a in range(n+1):poly=add(poly,scale(energy_pair(a,n-a),(-1)**(n-a)*comb(n,a)))
    for j in range(k):assert exponential_coefficient(poly,j,Fraction(0))==0
    lead=(-1)**k*exponential_coefficient(poly,k,Fraction(0))/factorial(n)
    correction=(-1)**k*exponential_coefficient(poly,k+1,Fraction(0))/factorial(n)
    assert lead==Fraction(str(g.coeff(u,n)))
    assert correction==-Fraction(str(R.coeff(u,n)))/2
    exact.append(dict(time_degree=n,leading=str(lead),first_correction=str(correction)))
parent='results/susskind-selected-rindler-response/response-checks.json'
assert hashlib.sha256((ROOT/parent).read_bytes()).hexdigest()=='814039c5cc5e3366db8d6957bfef66e80efee02adf409c7e1df0ac59e587fee0'
data=json.loads((ROOT/parent).read_text());rows=[]
for row in data['rows']:
    h=row['h'];values=[]
    for point in row['values']:
        u0=point['u'];base=cosh(u0)**-4;r=response(u0)
        mix=(point['trace']+point['selected'])/2
        values.append(dict(u=u0,trace_second_remainder=(point['trace']-base+h*r/2)/h**2,
                           selected_second_remainder=(point['selected']-base-h*r/2)/h**2,
                           mixture_error=mix-base,mixture_error_over_h2=(mix-base)/h**2))
    rows.append(dict(h=h,values=values))
result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),exact_taylor=exact,
    saved_packet_analysis=rows,environment=dict(python=platform.python_version(),sympy=sp.__version__),
    input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in [parent,
        'code/selected_rindler_response.py','code/check_mehler_trace_mixture.py']},
    scope='Exact finite trace Taylor checks and analysis of existing packet values; no old scientific driver rerun. The compact-time error order is established by the separate analytic density proof.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',exact_coefficients=16,saved_time_points=28,
    max_mixture_scaled_error=max(abs(p['mixture_error_over_h2']) for row in rows for p in row['values']),
    sha256=hashlib.sha256(OUT.read_bytes()).hexdigest())))
