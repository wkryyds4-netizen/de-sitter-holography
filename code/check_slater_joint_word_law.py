"""Exact consequences of the all-forward joint-word rule; no old main runs."""
from collections import Counter
from functools import lru_cache
from itertools import product
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,platform
import sympy as S
from slater_joint_site_polynomial import pairings,penalty

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-slater-joint-words/word-checks.json'
g=S.Symbol('gamma',positive=True)
@lru_cache(None)
def poly(w):
    if not w:return {0:1}
    if len(w)%2:return {}
    counts=Counter()
    for pi in pairings(tuple(range(len(w)))):
        kinds=[];mixed=0
        for i,j in pi:
            a,b=w[i],w[j]
            if a==b=='H':kinds.append('G')
            elif a=='+' or b=='-' or (a==b):break
            else:kinds.append('F');mixed+=('H' in (a,b))
        else:counts[mixed+4*penalty(pi,kinds)]+=1
    return dict(counts)
def expr(w):return sum(c*g**d for d,c in poly(w).items())
def adj(w):return ''.join({'+':'-','-':'+','H':'H'}[x] for x in reversed(w))
def gram_audit(value,degree):
    words=['']+[''.join(w) for n in range(1,degree+1) for w in product('H+-',repeat=n)]
    def ev(w):return sum(S.Integer(c)*value**d for d,c in poly(w).items())
    matrix=[[ev(adj(a)+b) for b in words] for a in words]
    remaining=list(range(len(words)));pivots=[]
    while remaining:
        pivot=next((i for i in remaining if matrix[i][i]),None)
        if pivot is None:
            assert all(matrix[i][j]==0 for i in remaining for j in remaining)
            break
        v=matrix[pivot][pivot];assert v>0
        pivots.append(dict(word=words[pivot],squared_norm=str(v)))
        remaining.remove(pivot)
        for i in remaining:
            for j in remaining:
                matrix[i][j]-=matrix[i][pivot]*matrix[pivot][j]/v
    return dict(gamma=str(value),degree=degree,raw_words=len(words),rank=len(pivots),pivots=pivots)

if OUT.exists():raise FileExistsError(OUT)
assert S.expand(expr('HHHHHH')-(5+6*g**4+3*g**8+g**12))==0
assert expr('---+++')==6*g**12
g_moments={}
for n in (2,4,6):
    moment=S.expand(sum((-g)**sum(x!='H' for x in w)*expr(''.join(w)) for w in product('H+-',repeat=n)))
    g_moments[str(n)]=dict(raw=str(S.factor(moment)),normalized=str(S.factor(moment/(1-g*g)**(n//2))))
assert S.factor(S.sympify(g_moments['4']['raw'],locals={'gamma':g})-(1-g*g)**2*(2+g*g+2*g**4))==0
audits=[gram_audit(v,d) for v in (S.Rational(1,2),S.Rational(2,3)) for d in (2,3)]
files=[Path(__file__),ROOT/'code/slater_joint_site_polynomial.py',ROOT/'results/susskind-slater-joint-words/site-checks.json']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
result=dict(status='passed',created_utc=datetime.now(timezone.utc).isoformat(),
 parameters=dict(gamma='q^(1/4)',rational_points=['1/2','2/3'],max_degree=6),
 sixth_words=[dict(word=''.join(w),gamma_polynomial=poly(''.join(w))) for w in product('H+-',repeat=6)],
 complement_moments=g_moments,gram_audits=audits,
 environment=dict(python=platform.python_version(),sympy=S.__version__),
 input_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in files},
 scope='All 729 sixth words as exact polynomials from the analytic rule; rational Gram positivity/rank at two points through degree3. Not an independent proof of all-word convergence or generic degree3 rank.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',g_moments=g_moments,gram=[{k:v for k,v in x.items() if k!='pivots'} for x in audits],sha256=sha(OUT))))
