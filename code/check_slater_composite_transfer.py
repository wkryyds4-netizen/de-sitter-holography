"""Fresh exact checks; refuses to overwrite the scientific output."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from datetime import datetime,timezone
from math import comb
import hashlib,json,platform,random,time
import slater_composite_sources as source
import slater_composite_transfer as transfer
import slater_covariant_transfer as old
import bilocal_multisource_chords as central

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-slater-composite-transfer/exact-checks.json'
if OUT.exists():raise FileExistsError(OUT)
OUT.parent.mkdir(parents=True,exist_ok=True)
START=time.perf_counter()
PARAMS=[(F(1,2),F(1,3)),(F(2,3),F(1,4))]
checks=0;source_checks=0;max_order=0
words=[]
for n in range(1,5):
    for tup in product('HWJ+-',repeat=n):
        word=''.join(tup)
        if (n+word.count('J'))%2==0:words.append(word)
rng=random.Random(20261108)
extra=set()
while len(extra)<72:
    word=''.join(rng.choice('HWJ+-') for _ in range(6))
    if word.count('J')==2:extra.add(word)
words+=sorted(extra)+['J+HJ+H','-JHHJ+','JWWJHH','JJJJ','WJWJHH','-JJJH+']
words=list(dict.fromkeys(words))
for index,word in enumerate(words):
    slots,blocks,types=source.layout(word)
    for pi in source.pairings(tuple(range(len(slots)))):
        if any(a<=i<j<b for a,b,_ in blocks for i,j in pi):continue
        kinds=[]
        for i,j in pi:
            left,right=slots[i],slots[j]
            if left==right=='H':kinds.append('G')
            elif left=='+' or right=='-':break
            else:kinds.append('F')
        else:
            kinds=tuple(kinds)
            actual=source.source_exponent(pi,kinds,blocks)
            predicted=transfer.geometric_exponent(pi,kinds,blocks)
            assert actual==predicted,(word,pi,kinds,actual,predicted)
            source_checks+=1
            max_order=max(max_order,max((m.bit_count() for m in actual),default=0))
    for gamma,lam in PARAMS:
        a=transfer.moment(word,gamma,lam);b=source.evaluate(word,gamma,lam)
        assert a==b,(word,gamma,lam,a,b);checks+=1
    # Bound cache memory without rerunning any sealed driver.
    source.local.cache_clear();source.source_exponent.cache_clear()
    if index%100==0:print('source/transfer',index,len(words),flush=True)

old_checks=0
for tup in product('HW+-',repeat=6):
    word=''.join(tup);gamma=F(1,2);d=F(1,3);lam=d*d
    a=transfer.moment(word,gamma,lam)
    b=old.moment(word.replace('W','U'),gamma,d)*d**word.count('W')
    assert a==b,(word,a,b);old_checks+=1
print('old transfer',old_checks,flush=True)

long=[]
for word in ['H'*24,'J'*8,'JHHJWWJHHJ','-JHHHHJ+','--JHHJ++']:
    stats={};gamma,lam=PARAMS[0];a=transfer.moment(word,gamma,lam,stats)
    if set(word)=={'H'}:
        n=len(word)//2;q=gamma**4
        b=sum((-1)**k*(comb(2*n,n-k)-(comb(2*n,n-k-1) if n>k else 0))*q**(k*(k+1)//2)
              for k in range(n+1))/(1-q)**n
    elif set(word)<=set('HWJ'):
        polynomial=central.compressed_moment(word)
        b=central.evaluate(polynomial,gamma**4,lam,word.count('J'))
    else:b=source.evaluate(word,gamma,lam)
    assert a==b,(word,a,b)
    long.append(dict(word=word,value=str(a),**stats))
    print('long',word,stats,flush=True)

preps=['','H','W','+','J','HH','HW','H+','WH','WW','W+','+H','+W','++']
gram_cases=[]
for gamma,lam in PARAMS:
    G=transfer.gram(preps,gamma,lam);n=len(G)
    assert all(G[i][j]==G[j][i] for i in range(n) for j in range(n))
    work=[row[:] for row in G];pivots=[]
    for i in range(n):
        pivot=work[i][i];assert pivot>0,(gamma,lam,i,pivot);pivots.append(str(pivot))
        for j in range(i+1,n):
            for k in range(i+1,n):work[j][k]-=work[j][i]*work[i][k]/pivot
    gram_cases.append(dict(gamma=str(gamma),lam=str(lam),LDL_pivots=pivots))

# State equality is stronger than vacuum equality. Pruning must preserve
# moments, but the unpruned reachable module is used to define its Gram.
unpruned=0
for word in ['JHHJ','-JJ+','J+HJ+H','WJWJ','JJJJ']:
    gamma,lam=PARAMS[0]
    assert transfer.vector(word,gamma,lam).get(transfer.EMPTY,0)==transfer.moment(word,gamma,lam)
    unpruned+=1

paths=['code/slater_composite_transfer.py','code/check_slater_composite_transfer.py',
       'code/slater_composite_sources.py','code/slater_covariant_transfer.py',
       'code/slater_joint_site_polynomial.py','code/bilocal_multisource_chords.py',
       'code/memory_body_source_chords.py','code/singlet_coefficient.py']
result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),
    environment=dict(python=platform.python_version(),arithmetic='fractions.Fraction; no external numerical dependency'),
    parameters=dict(formal_points=[[str(a),str(b)] for a,b in PARAMS],
        physical_relation='gamma=exp(-lambda); independent rational points are algebraic diagnostics, not physical DS samples',
        seed=20261108,random_words=sorted(extra)),
    source_exponents=source_checks,maximum_nonzero_source_order=max_order,
    source_moment_comparisons=checks,source_words=words,old_six_letter_words=old_checks,
    long_word_checks=long,Gram_preparations=preps,Gram_checks=gram_cases,
    unpruned_checks=unpruned,
    input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
    runtime_seconds=time.perf_counter()-START,
    scope='Exact finite checks of geometry and recurrence. Arbitrary fixed-word equivalence is proved analytically in the manuscript; no finite-N convergence rate or Hilbert metric on arbitrary automaton states is inferred from these tests.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({k:result[k] for k in ['status','source_exponents','source_moment_comparisons','old_six_letter_words','runtime_seconds']}),flush=True)
