"""Analytic exact-arithmetic error bounds for already saved spectral rows."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,math,platform
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1]
STAGE='susskind-neutral-selected-response'
OUT=ROOT/f'results/{STAGE}/spectral-bounds.json'
if OUT.exists():raise FileExistsError(OUT)
src=f'results/{STAGE}/response-checks.json'
assert hashlib.sha256((ROOT/src).read_bytes()).hexdigest()=='475535bab4df2947c1de0d886ae42672e478becd612ffb66338afc14b09af843'
data=json.loads((ROOT/src).read_text());mp.mp.dps=60;rows=[]
for row in data['rows']:
    eps=mp.mpf(str(row['epsilon']));q=mp.exp(-eps)
    for L in row['dimensions']:
        m=(2*L-1)//3
        assert m+2<L and 3*m<=2*L-1 and 2*m+2<=2*L-1
        for point in row['values']:
            u,v=map(lambda x:mp.mpf(str(x)),(point['u'],point['v']))
            z=2*max(abs(u),abs(v),abs(u-v))/mp.sqrt(eps*(1-q))
            eta=z**(m+1)/mp.factorial(m+1)
            chain=6*eta+6*eta**2+2*eta**3
            F=1-mp.sqrt(q)+4*mp.sqrt(q)/(1-q)
            tilted=F*(4*eta+2*eta**2)
            selected=tilted+mp.sqrt(q)*eps/4*(1+q)*chain
            rows.append(dict(epsilon=str(eps),L=L,degree=m,u=str(u),v=str(v),z=str(z),uniform_taylor_error=str(eta),chain_bound=str(chain),selected_bound=str(selected),refined=L==row['dimensions'][-1]))
products=[]
for q in [mp.mpf('.25'),mp.mpf('.64')]:
    for T in [128,192]:
        mu=6*q**(T+1)/((1-q)*(1-q**(T+1)))
        kernel=(q*q+4*q)*q**T/((1-q)*(1-q**(T+1)))
        products.append(dict(q=str(q),terms=T,triple_measure_log_error=str(3*mu+2*kernel)))
result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),rows=rows,product_tails=products,
    formula='m=floor((2L-1)/3); eta=z^(m+1)/(m+1)!; B error <=6eta+6eta^2+2eta^3, tilted pair error <= F(4eta+2eta^2). Real exponential Taylor remainder, q-OU L2 contraction and Gauss polynomial exactness. No scientific driver rerun.',
    scope='Analytic exact-arithmetic truncation bounds; displayed decimal evaluation uses mpmath and does not enclose floating roundoff or quadrature discretization.',
    environment=dict(python=platform.python_version(),mpmath=mp.__version__,digits=60),
    input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in [src,'code/record_neutral_spectral_bounds.py','code/neutral_selected_response.py']})
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',bounds=len(rows),max_refined_selected_bound=str(max(mp.mpf(r['selected_bound']) for r in rows if r['refined'])),sha256=hashlib.sha256(OUT.read_bytes()).hexdigest())))
