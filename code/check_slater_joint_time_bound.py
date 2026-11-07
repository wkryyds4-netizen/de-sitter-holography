"""Exact coefficient row/column sums and diagnostics of the uniform time tail."""
from pathlib import Path
from itertools import combinations,product
from math import comb
from datetime import datetime,timezone
import hashlib,json,platform
import mpmath as mp
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-slater-joint-words/time-bound-checks.json'
def choose(n,p):return comb(n,p) if n>=p>=0 else 0
def act(s,P,Q):
    for i in Q:
        if not s>>i&1:return None
        s^=1<<i
    for i in reversed(P):
        if s>>i&1:return None
        s^=1<<i
    return s
rows=[];count=0
for n,p in [(4,1),(6,1),(6,2),(8,2)]:
    m=n//2;sets=list(combinations(range(n),p));labels=list(product(sets,repeat=2))
    selected=[(P,Q) for P,Q in labels if min(P)>=m and max(Q)<m]
    L=comb(m,p)**2;D=comb(m,p)*comb(m+p,p)
    for st in range(1<<n):
        if st.bit_count()!=m:continue
        x=(st>>m).bit_count()
        full=sum(act(st,P,Q) is not None for P,Q in labels)
        creation=sum(act(st,P,Q) is not None for P,Q in selected)
        annihilation=sum(act(st,Q,P) is not None for P,Q in selected)
        assert full==D and creation==choose(m-x,p)**2 and annihilation==choose(x,p)**2
        assert max(creation,annihilation)<=L
        count+=3;rows.append(dict(N=n,p=p,state=st,x=x,D=D,L=L,creation=creation,annihilation=annihilation))
mp.mp.dps=70;tails=[]
for T,d in product([mp.mpf('0.5'),mp.mpf(1),mp.mpf(3)],[0,2,4]):
    R=1
    def ratio(n):return T*mp.sqrt(2*(n+d)+1)/(n+1)
    while ratio(R)>mp.mpf('0.5'):R+=1
    def term(n):return T**n*mp.sqrt(mp.factorial(2*(n+d))/(2**(n+d)*mp.factorial(n+d)))/mp.factorial(n)
    bound=term(R)/(1-ratio(R))
    partial=mp.fsum(term(n) for n in range(R,R+150))
    rest=term(R+150)/(1-ratio(R+150))
    assert partial<=bound and partial+rest<=bound
    tails.append(dict(T=str(T),insertion_degree=d,R=R,ratio_at_R=mp.nstr(ratio(R),35),
        bound=mp.nstr(bound,35),sum_150=mp.nstr(partial,35),remaining_majorant=mp.nstr(rest,35)))
if OUT.exists():raise FileExistsError(OUT)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
inputs=[Path(__file__),ROOT/'results/susskind-slater-joint-words/site-checks.json']
result=dict(status='passed',created_utc=datetime.now(timezone.utc).isoformat(),
    covariance_scalar_checks=count,covariance_rows=rows,tail_rows=tails,
    parameters=dict(fixtures=[(4,1),(6,1),(6,2),(8,2)],T=['0.5','1','3'],degrees=[0,2,4],decimal_precision=70),
    environment=dict(python=platform.python_version(),mpmath=mp.__version__),
    input_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in inputs},
    scope='Exact finite occupation counts for row/column covariance; high-precision scalar tail diagnostics, not interval certificates. All-degree diagram contraction and time convergence are analytic proofs.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',covariance_checks=count,tail_cases=len(tails),sha256=sha(OUT))))
