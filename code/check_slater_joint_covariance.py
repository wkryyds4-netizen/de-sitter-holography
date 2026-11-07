"""Read-only comparison of sealed CAR counts with the EGUE rank formula."""
from pathlib import Path
from math import comb
from fractions import Fraction as F
from datetime import datetime, timezone
import hashlib,json,platform
ROOT=Path(__file__).resolve().parent.parent
SOURCE=ROOT/'results/susskind-slater-joint-fourth/exact-checks.json'
OUT=ROOT/'results/susskind-slater-joint-fourth/covariance-checks.json'
def c(n,k):return comb(n,k) if 0<=k<=n else 0
def lam(n,m,r,k):return c(m-r,k)*c(n-m+k-r,k)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
rows=[]
for x in json.loads(SOURCE.read_text())['fixtures']:
    n,p=x['N'],x['p'];m=n//2
    terms=[lam(n,m,r,m-p)*lam(n,m,r,p)*(c(n,r)**2-c(n,r-1)**2) for r in range(min(p,m-p)+1)]
    K=F(sum(terms),c(n,m));T=lam(n,m,p,p);D=lam(n,m,0,p)
    assert K==x['H_crossed_K'] and T==x['crossed_eigenvalue_T'] and D==x['D']
    rows.append(dict(N=n,p=p,K=str(K),T=T,D=D,rank_sum_terms=terms))
if OUT.exists():raise FileExistsError(OUT)
result=dict(status='passed',created_utc=datetime.now(timezone.utc).isoformat(),checks=12,rows=rows,
    source='Kota nucl-th/0410069 equations (17), (25)-(28), unit Gaussian covariance',
    environment=dict(python=platform.python_version(),dependencies='standard library'),
    input_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in (Path(__file__),SOURCE)},
    scope='Independent scalar literature formula vs saved finite CAR enumeration; no rerun of scientific generator and no general K=L^2 assertion.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',checks=12,sha256=sha(OUT))))
