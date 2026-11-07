"""Exact rational sign certificate from Mehler moments and a cosine tail."""
from pathlib import Path
from fractions import Fraction as F
from math import factorial,comb
from decimal import Decimal,localcontext
from datetime import datetime,timezone
import hashlib,json,platform

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-slater-selected-state/sign-certificate.json'
if OUT.exists():raise FileExistsError(OUT)
q=F(1,4);s=F(1,2);r=q*q;t=F(2);K=20
# Monic continuous q-Hermite basis: E H_n=H_(n+1)+[n]_q H_(n-1).
qn=[(1-q**n)/(1-q) for n in range(2*K+4)]
norm=[F(1)]
for n in range(1,len(qn)):norm.append(norm[-1]*qn[n])
powers=[[F(1)]]
for n in range(2*K+2):
    new=[F(0)]*(n+2)
    for j,c in enumerate(powers[-1]):
        new[j+1]+=c
        if j:new[j-1]+=c*qn[j]
    powers.append(new)
def two(a,b):
    return sum((powers[a][n]*powers[b][n]*norm[n]*r**n for n in range(min(a,b)+1)),F(0))
def delta(m,selected):
    return sum((F((-1)**(m-a)*comb(m,a))*(
        (1-s)*two(a,m-a)+s*two(a+2,m-a) if selected else two(a,m-a))
        for a in range(m+1)),F(0))
def dec(x):
    with localcontext() as ctx:
        ctx.prec=38
        return str(Decimal(x.numerator)/Decimal(x.denominator))
rows=[]
for selected in (False,True):
    values=[delta(2*k,selected) for k in range(K+1)]
    assert values[0]==1 and all(x>=0 for x in values)
    expected=2*(1-r)+(s*(1+q)*(1-r)**2 if selected else 0)
    assert values[1]==expected
    assert all(delta(n,selected)==0 for n in (1,3,5,7))
    partial=sum(((-1)**k*t**(2*k)*values[k]/factorial(2*k) for k in range(K+1)),F(0))
    # |E-E'| <= 4/sqrt(1-q) = 8/sqrt(3) < 5, hence |t Delta| < 10.
    # Positive unit-mass spectral measure makes the cosine remainder bound
    # uniform before integration. Every bound here is rational.
    remainder=F(10**(2*K+2),factorial(2*K+2))
    lower,upper=partial-remainder,partial+remainder
    assert upper<0 if selected else lower>0
    rows.append(dict(state='one-selected' if selected else 'trace',
        even_delta_moments=[str(x) for x in values],partial_sum=str(partial),
        error_bound=str(remainder),lower=str(lower),upper=str(upper),
        decimal_lower_approx=dec(lower),decimal_upper_approx=dec(upper),
        sign='strictly negative' if selected else 'strictly positive'))
    print(rows[-1]['state'],dec(lower),dec(upper),flush=True)
inputs=[Path(__file__),ROOT/'results/susskind-slater-selected-state/exact-checks.json']
result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),
    parameters=dict(q=str(q),sqrt_q=str(s),probe_r=str(r),time=str(t),cosine_order=2*K,
                    spectral_difference_upper_bound='5',physical_lambda='log(2)/2'),
    environment=dict(python=platform.python_version(),arithmetic='fractions.Fraction; Decimal only for displayed approximations'),
    rows=rows,
    input_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
    scope='Rational strict-sign certification in the limiting two-energy measure. Remainder uses its proven positivity, unit mass, compact support and exact q-Hermite moments. Decimal endpoints are display approximations; rational endpoints are authoritative. No floating quadrature or eigensolver enters the sign test.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print('passed',hashlib.sha256(OUT.read_bytes()).hexdigest(),flush=True)
