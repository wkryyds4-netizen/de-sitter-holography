"""Exact diagram/transfer and new-sector checks, without old driver execution."""
from pathlib import Path
from itertools import product
from fractions import Fraction as F
from datetime import datetime,timezone
import hashlib,json,platform
from slater_covariant_transfer import step,moment,paired_moment,basis_inner,inner,add,scale
from slater_open_transfer import moment as previous_moment
ROOT=Path(__file__).resolve().parent.parent
DIR=ROOT/'results/susskind-slater-covariant';OUT=DIR/'word-checks.json'
if OUT.exists():raise FileExistsError(OUT)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def op(v,name,g,d):return step(v,name,g,d)
def apply(v,sequence,g,d):
    for name in reversed(sequence):v=op(v,name,g,d)
    return v
def positive_pivots(words,g,d):
    mat=[[basis_inner(a,b,g,d) for b in words] for a in words]
    assert all(mat[i][j]==mat[j][i] for i in range(len(words)) for j in range(len(words)))
    pivots=[]
    for k in range(len(words)):
        pivot=mat[k][k];assert pivot>0,(g,d,k,pivot)
        pivots.append(str(pivot))
        for i in range(k+1,len(words)):
            for j in range(i,len(words)):
                mat[i][j]-=mat[i][k]*mat[k][j]/pivot;mat[j][i]=mat[i][j]
    return pivots
counts=dict(sixth_words=0,previous_restriction=0,algebra_relations=0,projection_orthogonality=0,Gram_entries=0)
rows=[]
for g,d in ((F(1,2),F(1,2)),(F(2,3),F(1,2))):
    q=g**4;s=g*g;c=-d*g
    for letters in product('HU+-',repeat=6):
        word=''.join(letters)
        assert moment(word,g,d)==paired_moment(word,g,d),word
        counts['sixth_words']+=1
    for letters in product('H+-',repeat=6):
        word=''.join(letters)
        assert moment(word,g,d)==previous_moment(word,g)
        counts['previous_restriction']+=1
    gram=[]
    for n in (1,2,3):
        words=[''.join(w) for w in product('HUP',repeat=n)]
        gram.append(dict(level=n,dimension_checked=len(words),positive_pivots=positive_pivots(words,g,d)))
        counts['Gram_entries']+=len(words)**2
    for n in range(4):
        for letters in product('HUP',repeat=n):
            w=''.join(letters);v={w:F(1)}
            for ai,Ai in (('a','A'),('v','V')):
                for aj,Aj in (('a','A'),('v','V')):
                    difference=add(apply(v,[ai,Aj],g,d),scale(-q,apply(v,[Aj,ai],g,d)))
                    assert difference==(v if ai==aj else {})
                    counts['algebra_relations']+=1
                cc=g if ai=='a' else -d*g
                assert add(apply(v,[ai,'P'],g,d),scale(-q,apply(v,['P',ai],g,d)))==scale(cc*q**n,v)
                assert add(apply(v,['m',Ai],g,d),scale(-q,apply(v,[Ai,'m'],g,d)))==scale(cc*q**n,v)
                counts['algebra_relations']+=2
            assert add(apply(v,['m','P'],g,d),scale(-q,apply(v,['P','m'],g,d)))==scale(q**n,v)
            counts['algebra_relations']+=1
    residual={'U':F(1),'H':-d*s/(1-s),'P':d*g/(1-s)}
    residual_norm=1-d*d*s/(1-s)
    assert inner(residual,residual,g,d)==residual_norm>0
    assert step(residual,'a',g,d)=={}
    for n in range(5):
        for letters in product('HP',repeat=n):
            assert inner(residual,{''.join(letters):F(1)},g,d)==0
            counts['projection_orthogonality']+=1
    formulas={'-UU+':1+c*c*(1+q),'U+-U':c*c,'-U+U':q*(1+c*c),
              '-+UU':F(1),'UHUH':q,'UUHH':F(1)}
    for word,value in formulas.items():assert moment(word,g,d)==value,word
    comm=moment('UHHU',g,d)-moment('UHUH',g,d)-moment('HUHU',g,d)+moment('HUUH',g,d)
    assert comm==2*(1-q)
    rows.append(dict(gamma=str(g),sqrt_lambda=str(d),q=str(q),first_level_residual_norm=str(residual_norm),
                     mixed_fourth={w:str(v) for w,v in formulas.items()},commutator_norm=str(comm),gram=gram))
inputs=[Path(__file__),ROOT/'code/slater_covariant_transfer.py',ROOT/'code/slater_open_transfer.py',
        DIR/'site-checks.json',ROOT/'results/susskind-slater-transfer/local-evidence.json']
result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),counts=counts,rows=rows,
            parameters=dict(rational_pairs=[['1/2','1/2'],['2/3','1/2']],physical_relation_imposed=False,
                            meaning='Exact independent gamma and d test points for algebraic identities; microscopic results impose gamma=exp(-lambda), d=sqrt(lambda).'),
            environment=dict(python=platform.python_version(),arithmetic='fractions.Fraction'),
            input_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in inputs},
            scope='Exact transfer vs independent pairing enumeration and low-level relations/positive Grams. Positivity at all physical parameters follows from microscopic moment convergence, not finite test points.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',counts=counts,sha256=sha(OUT))))
