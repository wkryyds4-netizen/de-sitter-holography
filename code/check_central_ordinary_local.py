"""Independent supplemental local coverage after detecting a vacuous layout loop.

The first checker/output is retained unchanged; this driver requires nonzero
two-, three- and four-source coverage and writes a separate exclusive file.
"""
from pathlib import Path
from itertools import product,combinations
from datetime import datetime,timezone
from collections import Counter
import hashlib,json,platform,time
from central_ordinary_sources import central_exponent,exponent
from memory_joint_regression import pairings,region,cr

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-central-ordinary/local-coverage.json'

def pairs_of_intervals(n):
    for a in range(n):
        for b in range(a+1,n+1):
            for c in range(b,n):
                for d in range(c+1,n+1):yield ((a,b),(c,d))

def expected(pi,marks):
    out={0:-4*cr(pi)};sets=[];types=[]
    for k,(a,b,t) in enumerate(marks):
        inside,opened=region(pi,a,b);types.append(t);sets.append(opened)
        linear=-4*cr(opened)-2*len(inside) if t=='D' else (b-a)**2-4*cr(opened)
        out[1<<k]=linear
    for a,b in combinations(range(len(marks)),2):
        common=sets[a]&sets[b];d=len(common)
        coefficient=2*d*(d-1) if types[a]==types[b]=='D' else d*d if types[a]==types[b]=='L' else d*(d-1)
        out[(1<<a)|(1<<b)]=coefficient-4*cr(common)
    return {k:v for k,v in out.items() if v}

def main():
    if OUT.exists():raise FileExistsError(OUT)
    start=time.perf_counter();counts=Counter();samples=[]
    for n in (2,4,6,8):
        regions=list(pairs_of_intervals(n)) if n<=6 else [((0,2),(2,4),(4,8)),((0,2),(2,4),(4,6),(6,8))]
        for rs in regions:
            for types in product(('D','L'),repeat=len(rs)):
                marks=tuple((a,b,t) for (a,b),t in zip(rs,types))
                for pi in pairings(tuple(range(n))):
                    if any(t=='L' and a<=x<y<b for a,b,t in marks for x,y in pi):continue
                    actual=central_exponent(n,pi,marks,True);predicted=expected(pi,marks)
                    assert actual==predicted,(n,marks,pi,actual,predicted)
                    counts[f'{len(marks)}_source_cases']+=1
                    counts['with_internal_pairs']+=any(a<=x<y<b for a,b,t in marks for x,y in pi)
                    counts['negative_nonconstant_coefficients']+=sum(v<0 for k,v in actual.items() if k)
                    if len(samples)<12 and any(a<=x<y<b for a,b,t in marks for x,y in pi):
                        samples.append(dict(length=n,marks=marks,pairing=pi,central={str(k):str(v) for k,v in actual.items()}))
                exponent.cache_clear()
        print('Covered slots',n,dict(counts),flush=True)
    assert all(counts[f'{r}_source_cases']>0 for r in (2,3,4))
    assert counts['with_internal_pairs']>0 and counts['negative_nonconstant_coefficients']>0
    paths=['code/check_central_ordinary_local.py','code/central_ordinary_sources.py',
           'code/memory_body_source_chords.py','code/memory_joint_regression.py',
           'code/memory_size_marked_chords.py']
    result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),counts=dict(counts),
        parameters={'all_disjoint_two_intervals_slots':[2,4,6],
                    'eight_slot_regions':[[[0,2],[2,4],[4,8]],[[0,2],[2,4],[4,6],[6,8]]],
                    'types':'all D/L assignments with L-internal pairs excluded'},
        samples=samples,environment={'python':platform.python_version(),'arithmetic':'fractions.Fraction'},
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        runtime_seconds=time.perf_counter()-start,
        scope='Nonvacuous exact unrestricted-site vs central geometric source comparisons, including all masks at three/four sources. Complements, does not rewrite, the first output.')
    with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({'status':'passed','counts':dict(counts),'seconds':result['runtime_seconds']}),flush=True)

if __name__=='__main__':main()
