"""Exact graded canonical-undressing checks on the open transfer core."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from datetime import datetime,timezone
from collections import Counter
import hashlib,json,platform
from slater_open_transfer import step,add,scale,inner,qnum,qfac
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-slater-transfer/canonical-checks.json'
def beta(v,g,dagger=False):
    out={}
    for w,c in v.items():
        if dagger:part=scale(g**(-2*len(w)),step({w:c},'P',g))
        else:part=scale(g**(-2*(len(w)-1)),step({w:c},'m',g))
        out=add(out,part)
    return out
if OUT.exists():raise FileExistsError(OUT)
counts=Counter();rows=[]
for g in [F(1,2),F(2,3)]:
    q=g**4
    for n in range(5):
        for letters in product('HP',repeat=n):
            w=''.join(letters);v={w:F(1)}
            Q=add(step(step(v,'P',g),'m',g),scale(-q,step(step(v,'m',g),'P',g)))
            assert Q==scale(q**n,v);counts['grading_from_M']+=1
            ccr=add(beta(beta(v,g,True),g),scale(-1,beta(beta(v,g),g,True)))
            assert ccr==v;counts['canonical_CCR']+=1
            cross=add(step(beta(v,g,True),'a',g),scale(-g*g,beta(step(v,'a',g),g,True)))
            assert cross==scale(g**(2*n+1),v);counts['same_J_cross']+=1
    for n in range(9):
        v={'P'*n:F(1)}
        actual=inner(step(v,'P',g),step(v,'P',g),g)/inner(v,v,g)
        assert actual==(n+1)*q**n;counts['selected_shift_norm']+=1
        rows.append(dict(gamma=str(g),level=n,selected_creation_norm_squared=str(actual)))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
files=[Path(__file__),ROOT/'code/slater_open_transfer.py',ROOT/'results/susskind-slater-transfer/exact-checks.json']
j=dict(status='passed',created_utc=datetime.now(timezone.utc).isoformat(),counts=dict(counts),rows=rows,
 parameters=dict(gamma=['1/2','2/3'],word_levels=list(range(5)),norm_levels=list(range(9))),
 environment=dict(python=platform.python_version(),dependencies='standard library'),
 input_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in files},
 scope='Exact finite-core algebra checks; Type I closure, general norms and spectral decomposition are proved analytically in companion text. No microscopic occupation identification beyond the proved joint representation.')
with OUT.open('x',encoding='utf-8') as f:json.dump(j,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',counts=dict(counts),sha256=sha(OUT))))
