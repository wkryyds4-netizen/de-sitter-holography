"""Exact transfer vs saved moments, algebra relations and primitive towers."""
from pathlib import Path
from itertools import product,combinations
from fractions import Fraction as F
from datetime import datetime,timezone
from collections import Counter
import hashlib,json,platform,time
from slater_open_transfer import step,add,scale,moment,basis_inner,inner,qnum,qfac,primitive
from slater_joint_site_polynomial import pairings,penalty
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-slater-transfer/exact-checks.json'
SOURCE=ROOT/'results/susskind-slater-joint-words/word-checks.json'
def independent(w,g):
    total=F(0);q=g**4
    for pi in pairings(tuple(range(len(w)))):
        kinds=[];b=0
        for i,j in pi:
            if w[i]==w[j]=='H':kinds.append('G')
            elif w[i]=='+' or w[j]=='-' or w[i]==w[j]:break
            else:kinds.append('F');b+=('H' in (w[i],w[j]))
        else:total+=g**b*q**penalty(pi,kinds)
    return total
def assert_same(a,b,context):assert a==b,(context,a,b)
if OUT.exists():raise FileExistsError(OUT)
start=time.perf_counter();counts=Counter();summaries=[]
saved=json.loads(SOURCE.read_text())
for g in (F(1,2),F(2,3)):
    q=g**4;s=g*g
    for row in saved['sixth_words']:
        expect=sum(c*g**int(d) for d,c in row['gamma_polynomial'].items())
        assert_same(moment(row['word'],g),expect,row['word']);counts['saved_sixth_words']+=1
    words8={'H'*8,'----++++','--HHHH++','-H-H+H+H'}
    for positions in combinations(range(8),2):
        for signs in product('+-',repeat=2):
            w=list('H'*8)
            for i,c in zip(positions,signs):w[i]=c
            words8.add(''.join(w))
    for w in sorted(words8):
        assert_same(moment(w,g),independent(w,g),w);counts['independent_eighth_words']+=1
    for n in range(5):
        for letters in product('HP',repeat=n):
            w=''.join(letters);v={w:F(1)}
            for low,up,rhs in [('a','A',1),('a','P',g*q**n),('m','A',g*q**n),('m','P',q**n)]:
                left=add(step(step(v,up,g),low,g),scale(-q,step(step(v,low,g),up,g)))
                assert_same(left,scale(rhs,v),(w,low,up));counts['q_relations']+=1
            pr=primitive(v,g)
            assert_same(step(pr,'a',g),{},('primitive',w))
            assert_same(primitive(pr,g),pr,('idempotent',w));counts['primitive_projection']+=2
            expansion={};lower=v
            for r in range(n+1):
                component=scale(F(1)/qfac(r,q),primitive(lower,g))
                for _ in range(r):component=step(component,'A',g)
                expansion=add(expansion,component);lower=step(lower,'a',g)
            assert_same(expansion,v,('tower decomposition',w));counts['tower_decomposition']+=1
    for n in range(4):
        left=[''.join(w) for w in product('HP',repeat=n)]
        right=[''.join(w) for w in product('HP',repeat=n+1)]
        for u in right:
            for v in left:
                for low,up in [('a','A'),('m','P')]:
                    assert_same(inner({u:1},step({v:1},up,g),g),inner(step({u:1},low,g),{v:1},g),(u,v,up))
                    counts['adjoint_entries']+=1
        for u in left:
            for v in left:
                assert_same(basis_inner(u,v,g),basis_inner(v,u,g),(u,v));counts['metric_symmetry']+=1
    for n in range(9):
        z={'P'+'H'*n:F(1),'H'*(n+1):-g}
        expect=add({'P'+'H'*(n+1):F(1),'H'*(n+2):-g},
                   scale(qnum(n,q),{'P'+'H'*(n-1):F(1),'H'*n:-g}) if n else {})
        assert_same(step(z,'H',g),expect,('residual Jacobi',n));counts['residual_Jacobi']+=1
        numerator=basis_inner('H'*(n+1),'H'*n+'P',g)
        assert_same(numerator,g*(n+1)*q**n*qfac(n,q),('compressed selected',n));counts['compressed_selected']+=1
    spectral=[]
    for r in range(9):
        actual=moment('-'+'H'*r+'+',g)
        expected=(1-s)*moment('H'*r,g)+s*moment('H'*(r+2),g)
        assert_same(actual,expected,('selected spectral',r));counts['selected_spectral_moments']+=1
        spectral.append(dict(power=r,moment=str(actual)))
    summaries.append(dict(gamma=str(g),q=str(q),selected_spectral_moments=spectral))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
files=[Path(__file__),ROOT/'code/slater_open_transfer.py',ROOT/'code/slater_joint_site_polynomial.py',SOURCE,
 ROOT/'results/susskind-slater-joint-words/local-evidence.json']
result=dict(status='passed',created_utc=datetime.now(timezone.utc).isoformat(),counts=dict(counts),summaries=summaries,
 parameters=dict(gamma=['1/2','2/3'],relations_max_level=4,adjoints_max_level=3,spectral_max_power=8),
 environment=dict(python=platform.python_version(),dependencies='standard library'),
 input_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in files},elapsed_seconds=time.perf_counter()-start,
 scope='Exact transfer vs sealed sixth moments and independently enumerated selected eighth moments; coefficient-level Wick relations, projection/tower decomposition, low-level metric adjunction and selected spectral identities. Analytic general-level proofs are separate.')
OUT.parent.mkdir(parents=True,exist_ok=True)
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',counts=dict(counts),elapsed=result['elapsed_seconds'],sha256=sha(OUT))))
