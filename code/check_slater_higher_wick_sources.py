"""Fresh finite local certificate and higher-Wick density checks; exclusive output."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from collections import Counter
from datetime import datetime,timezone
import hashlib,json,platform,time
import slater_higher_wick_sources as h
from slater_composite_sources import source_exponent,local,evaluate,layout

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-slater-higher-wick/exact-checks.json'

def layouts(n,pi,forbidden,pos=0):
    yield ()
    for a in range(pos,n):
        for b in range(a+1,n+1):
            if set(range(a,b))&forbidden or any(a<=u<v<b for u,v in pi):continue
            for kind in ('D','L'):
                for tail in layouts(n,pi,forbidden,b):yield ((a,b,kind),)+tail

def clean(x):return {m:v for m,v in x.items() if v}

def main():
    if OUT.exists():raise FileExistsError(OUT)
    start=time.perf_counter();counts=Counter();table=[]
    for n,pis in [(2,[((0,1),)]),(4,[((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2))])]:
        for pi in pis:
            for le,re in product((False,True),repeat=2):
                if not(le or re):continue
                exterior=({0} if le else set())|({n-1} if re else set())
                kinds=tuple('F' if set(edge)&exterior else 'G' for edge in pi)
                for marks in layouts(n,pi,exterior):
                    counts['local_layouts']+=1
                    # All nonconstant nilpotent coefficients, not only linear marks.
                    if n==4:
                        a=h.collision_part(pi,kinds,marks,(0,1))
                        b=h.collision_part(pi,('G','G'),marks,(0,1))
                        assert a==b,(pi,exterior,marks,a,b)
                        counts['cross_source_equalities']+=1
                    for e,edge in enumerate(pi):
                        a=h.collision_part(pi,kinds,marks,(e,e))
                        b=Counter(h.collision_part(pi,('G',)*len(pi),marks,(e,e)))
                        if kinds[e]=='F':
                            for k,(u,v,t) in enumerate(marks):
                                if t=='L':b[1<<k]-=sum(u<=x<v for x in edge)
                        assert a==clean(b),(pi,exterior,marks,e,a,b)
                        counts['own_source_equalities']+=1
                    # Independent all-G boundary-W insertion, including both boundaries.
                    extended=tuple(sorted(marks+tuple((x,x+1,'L') for x in exterior)))
                    boundary={i for i,(a,b,t) in enumerate(extended) if a in exterior}
                    boundary_mask=sum(1<<i for i in boundary)
                    if n==4:
                        cross=h.collision_part(pi,('G','G'),extended,(0,1))
                        assert not any(m&boundary_mask for m in cross),(pi,extended,cross)
                        counts['boundary_W_cross_zero']+=1
                    for e,edge in enumerate(pi):
                        own=h.collision_part(pi,('G',)*len(pi),extended,(e,e))
                        actual={m:v for m,v in own.items() if m&boundary_mask}
                        expected=Counter()
                        for k in boundary:
                            x=extended[k][0]
                            if x not in edge:continue
                            expected[1<<k]=1;y=edge[0] if x==edge[1] else edge[1]
                            for j,(a,b,t) in enumerate(extended):
                                if t=='L' and a<=y<b:expected[(1<<k)|(1<<j)]=1
                        assert actual==dict(expected),(pi,extended,actual,expected)
                        counts['boundary_W_own_equalities']+=1
                    table.append({'pairing':pi,'exterior':sorted(exterior),'marks':marks})
                source_exponent.cache_clear();local.cache_clear()
    print('Universal local certificate',dict(counts),flush=True)

    cases=[(0,(),(),()),(2,((0,2,'D'),),((0,2),),()),
           (4,((0,3,'D'),),((0,3),),()),
           (6,((0,3,'D'),(3,6,'D')),((0,3),(3,6)),()),
           (6,((0,4,'D'),),((0,4),),()),
           (8,((0,4,'D'),(4,8,'D')),((0,4),(4,8)),()),
           (6,((0,2,'L'),),((0,3),),()),
           (6,((1,3,'L'),),((0,3),),()),
           (6,((0,2,'L'),(3,5,'L')),((0,3),(3,6)),()),
           (6,((0,1,'L'),(2,4,'D'),(4,6,'L')),((0,1),(2,4),(4,6)),()),
           (8,((1,3,'L'),(5,7,'D')),((1,4),(4,7)),()),
           (6,((2,4,'L'),),((2,5),),((0,5,1),)),
           (6,((2,4,'L'),),((2,5),),((0,5,3),)),
           (6,((2,4,'L'),),((2,5),),((0,5,5),)),
           (8,((1,3,'L'),(5,7,'L')),((1,4),(5,8)),((0,4,2),(4,8,2)))]
    rows=[]
    for gamma,lam in [(F(1,2),F(1,3)),(F(2,3),F(1,4))]:
        for length,marks,wr,cr in cases:
            shifts=tuple(F(i+1,7) for i in range(len(marks)))
            values={mode:h.moment(length,marks,gamma,lam,wr,cr,shifts,mode)
                    for mode in ('trace','selected','chi')}
            expected=(1-gamma**2*(1+lam))*values['trace']+gamma**2*values['chi']
            assert values['selected']==expected,(length,marks,wr,cr,values,expected)
            rows.append(dict(length=length,marks=marks,wick_regions=wr,chaos_regions=cr,
                gamma=str(gamma),lam=str(lam),shifts=list(map(str,shifts)),
                moments={k:str(v) for k,v in values.items()}))
            source_exponent.cache_clear();local.cache_clear()
        print('Higher density cases',str(gamma),len(rows),flush=True)
    legacy=0
    for word in ('','H','W','J','WW','JJ','JHW','WJHW','JHHJ'):
        slots,marks,types=layout(word);g=F(1,2);lam=F(1,3)
        shifts=tuple(lam if t=='W' else 2*lam*(1-g**4)/(1+g**4) for t in types)
        for mode,prefix,suffix in [('trace','',''),('selected','-','+')]:
            assert h.moment(len(slots),marks,g,lam,shifts=shifts,mode=mode)==evaluate(prefix+word+suffix,g,lam)
            legacy+=1
        source_exponent.cache_clear();local.cache_clear()
    paths=['code/slater_higher_wick_sources.py','code/check_slater_higher_wick_sources.py',
           'code/slater_composite_sources.py','code/slater_joint_site_polynomial.py',
           'code/memory_body_source_chords.py','code/memory_size_marked_chords.py',
           'code/bilocal_multisource_chords.py']
    result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),
        environment={'python':platform.python_version(),'arithmetic':'fractions.Fraction'},
        parameters={'local_endpoint_counts':[2,4],'local_source_masks':'all nonconstant masks',
                    'moments':'two independent rational (gamma,lambda) points; not physical positivity samples'},
        counts=dict(counts),local_cases=table,density_rows=rows,legacy_comparisons=legacy,
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        runtime_seconds=time.perf_counter()-start,
        scope='Finite exhaustive local certificate plus selected higher-degree examples; proof of reduction is in companion manuscript. No ordinary internal-pair marked block or uniform-in-degree actual tail is certified.')
    OUT.parent.mkdir(parents=True,exist_ok=True)
    with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({'status':'passed','density':len(rows),'legacy':legacy,'runtime_seconds':result['runtime_seconds']}),flush=True)

if __name__=='__main__':main()
