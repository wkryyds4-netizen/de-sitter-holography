"""Exact common-source versus three-species pairing identities.

Pure J/H moments are checked only as a connection to the sealed old
multimark result. New content is W and its joint H/J moments. C=1.
"""
from collections import Counter
from datetime import datetime,timezone
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import hashlib,json,platform,time

from bilocal_multisource_chords import raw_moment,compressed_moment,evaluate,encode
from singlet_multimark import raw_moment as old_J_moment

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-bilocal-four-point/joint-moment-checks.json'
PARAMETERS=dict(C=1,counts_J_W_H=[[0,2,0],[0,4,0],[0,2,2],[0,4,2],[0,6,0],
    [2,2,0],[2,2,2],[1,2,0],[1,2,2],[3,2,0],[1,1,1],[2,1,1],
    [0,3,1],[1,1,3],[1,3,1],[2,0,0],[2,0,2],[4,0,0]],
    cyclic_representatives=True,formal_points=[['1/3','1/4'],['1/2','1/2'],['3/4','1/8']],
    source='Exact local D/L ratios with source-dependent singleton denominators; all nilpotent orders retained',
    expected_rule='Pair equal species only, lengths H=1,W=1,J=2, crossing q^(length product)',
    scope='Formal rational lambda,q are independent. Physical use imposes q=exp(-4lambda).')

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def digest(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def words(j,w,h):
    n=j+w+h
    for posj in combinations(range(n),j):
        remaining=[i for i in range(n) if i not in posj]
        for posw in combinations(remaining,w):
            yield ''.join('J' if i in posj else 'W' if i in posw else 'H' for i in range(n))

def canonical(word):return min(word[i:]+word[:i] for i in range(len(word)))

def main():
    if OUT.exists():raise FileExistsError(OUT)
    started=time.perf_counter();rows=[];statistics=Counter();allwords=set();oldchecks=0
    for counts in PARAMETERS['counts_J_W_H']:
        distinct=set(words(*counts));reps=sorted(set(map(canonical,distinct)))
        allwords.update(distinct)
        for word in reps:
            actual,stats=raw_moment(word);expected=compressed_moment(word)
            assert actual==expected,(word,encode(actual),encode(expected))
            if 'W' not in word:
                prior,_=old_J_moment(word.replace('J','X'))
                assert actual==prior,(word,'old-multimark')
                oldchecks+=1
            statistics.update({k:v for k,v in stats.items() if k!='maximum_source_order'})
            rows.append(dict(word=word,cyclic_orbit_size=len({word[i:]+word[:i] for i in range(len(word))}),
                actual_numerator=encode(actual),denominator_power=word.count('J'),
                expected_numerator=encode(expected),zero_moment=not actual,statistics=stats,exact=True))
        print(json.dumps(dict(counts=counts,cyclic_representatives=len(reps),elapsed_seconds=time.perf_counter()-started)),flush=True)
    diagnostics=[]
    for qs,ls in PARAMETERS['formal_points']:
        q,lam=F(qs),F(ls);w=4*lam+16*lam*lam*q/(1+q)
        values={word:evaluate(raw_moment(word)[0],q,lam,word.count('J'))
                for word in ['HH','WW','JJ','JW','WH','JJWW','JWJW','WWWW','JJJJ','WHWH','WWHH']}
        assert values['WW']==lam and values['JJ']==w and values['JW']==values['WH']==0
        assert values['WWWW']==(2+q)*lam*lam and values['JJJJ']==(2+q**4)*w*w
        assert values['JJWW']==w*lam and values['JWJW']==q*q*w*lam
        comm=2*(values['JJWW']-values['JWJW']); gram=values['JJWW']**2-values['JWJW']**2
        assert comm==2*w*lam*(1-q*q)>0 and gram>0
        assert 2*(values['WWHH']-values['WHWH'])==2*lam*(1-q)>0
        diagnostics.append(dict(q=qs,lambda_value=ls,moments={k:str(v) for k,v in values.items()},
            commutator_squared_JW=str(comm),product_Gram_determinant_JW=str(gram),
            normalized_commutator_squared_JW=str(comm/(w*lam))))
    names=['code/bilocal_multisource_chords.py','code/check_bilocal_multisource_moments.py',
           'code/memory_body_source_chords.py','code/memory_size_marked_chords.py',
           'code/singlet_coefficient.py','code/singlet_multimark.py',
           'results/susskind-singlet-multimark/local-evidence.json',
           'results/susskind-memory-joint-regression/local-evidence.json']
    result=dict(schema_version=1,status='passed',created_utc=datetime.now(timezone.utc).isoformat(),
        parameters=PARAMETERS,parameter_sha256=digest(PARAMETERS),records=rows,diagnostics=diagnostics,
        counts=dict(cyclic_representatives=len(rows),words_covered=sum(r['cyclic_orbit_size'] for r in rows),
            zero_representatives=sum(r['zero_moment'] for r in rows),old_JH_interface_checks=oldchecks,
            maximum_source_order_observed=max(r['statistics']['maximum_source_order'] for r in rows),
            formal_Gram_diagnostics=len(diagnostics),**statistics),
        input_sha256={p:sha(ROOT/p) for p in names},
        environment=dict(python=platform.python_version(),platform=platform.platform(),arithmetic='Fraction and exact integer polynomials'),
        elapsed_seconds=time.perf_counter()-started,
        scope='Exact fixed-word model diagram identities, retaining all source orders. '
            'The arbitrary-word theorem and actual pinching transfer need the proofs; '
            'full finite-relative-time bilocal four-point convergence is not inferred.')
    assert result['counts']['words_covered']==len(allwords)
    OUT.parent.mkdir(parents=True,exist_ok=True)
    with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(status='passed',counts=result['counts'],elapsed_seconds=result['elapsed_seconds'],sha256=sha(OUT))))

if __name__=='__main__':main()
