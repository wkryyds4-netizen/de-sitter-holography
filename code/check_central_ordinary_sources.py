"""Fresh exact central-source, Gram and finite-N calibration evidence."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement,product
from math import prod
from collections import Counter
from datetime import datetime,timezone
import hashlib,json,platform,time
import central_ordinary_sources as c
from memory_joint_regression import pairings,region,cr
from memory_body_source_chords import mul
from slater_higher_wick_sources import exponential

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-central-ordinary/exact-checks.json'

def gaussian_transform(internals,lam):
    out={}
    for mask in range(1<<len(internals)):
        chosen=[i for a,i in enumerate(internals) if mask>>a&1]
        value=c.limiting_singleton_coefficient(lam,chosen)
        if value:out[mask]=value
    return out

def layouts(n,pos=0):
    yield ()
    for a in range(pos,n):
        for b in range(a+1,n+1):
            for typ in ('D','L'):
                for tail in layouts(n,b):
                    if len(tail)<1:yield ((a,b,typ),)+tail

def main():
    if OUT.exists():raise FileExistsError(OUT)
    start=time.perf_counter();counts=Counter();gram=[];finite=[]
    # Independent Gaussian moments of the Fourier variable vs nilpotent exponential.
    for r in range(7):
        for ints in product((0,1,2),repeat=r):
            for lam in (F(1,4),F(2,3)):
                K={ (1<<a)|(1<<b):-4*ints[a]*ints[b]
                    for a,b in combinations(range(r),2) if ints[a]*ints[b]}
                assert exponential(K,r,lam)==gaussian_transform(ints,lam)
                counts['fourier_nilpotent_transforms']+=1
    print('Fourier coefficients',counts['fourier_nilpotent_transforms'],flush=True)
    # Exact two-source geometric formula vs the unrestricted site-matrix engine.
    for n in (2,4,6):
        for marks in layouts(n):
            if len(marks)!=2:continue
            for pi in pairings(tuple(range(n))):
                if any(t=='L' and a<=x<y<b for a,b,t in marks for x,y in pi):continue
                K=c.central_exponent(n,pi,marks,True)
                expected={0:-4*cr(pi)};sets=[];types=[]
                for k,(a,b,typ) in enumerate(marks):
                    inside,opened=region(pi,a,b);types.append(typ)
                    if typ=='D':linear=-4*cr(opened)-2*len(inside);incident=opened
                    else:linear=(b-a)**2-4*cr(opened);incident=opened
                    if linear:expected[1<<k]=linear
                    sets.append(incident)
                common=sets[0]&sets[1];v=len(common)
                connected=(2*v*(v-1) if types==['D','D'] else v*v if types==['L','L'] else v*(v-1))-4*cr(common)
                if connected:expected[3]=connected
                expected={m:v for m,v in expected.items() if v}
                assert K==expected,(n,marks,pi,K,expected)
                counts['site_vs_geometric_central_sources']+=1
        c.exponent.cache_clear()
    print('Source comparisons',counts['site_vs_geometric_central_sources'],flush=True)
    for q,lam in ((F(1,4),F(1,3)),(F(2,3),F(1,4))):
        m=c.CentralModel(q,lam)
        for n in range(2,6):
            for degree in range(5):
                if n+degree>8 or (n+degree)%2:continue
                assert m.moment(m.H(degree),m.G(n))==m.moment(m.H(degree),m.g(n))
                counts['energy_regressions']+=1
        for n,k in combinations_with_replacement(range(2,6),2):
            for a,b in ((0,0),(1,1),(2,0)):
                if n+k+a+b>10 or (n+k+a+b)%2:continue
                R,S=m.residual(n),m.resummed(n);T,U=m.residual(k),m.resummed(k)
                values=[m.moment(R,T,a,b),m.moment(R,U,a,b),m.moment(S,T,a,b),m.moment(S,U,a,b)]
                assert len(set(values))==1,(q,lam,n,k,a,b,values)
                gram.append(dict(q=str(q),lam=str(lam),n=n,k=k,a=a,b=b,RR_RS_SR_SS=list(map(str,values))))
                counts['gram_equalities']+=3
        # The omitted connected charge term really changes raw ordinary moments.
        full=c.JointModel(q,lam)
        gap=full.moment(full.G(2),full.G(2))-m.moment(m.G(2),m.G(2))
        assert gap==4*lam
        counts['nonzero_full_vs_central_G2_gaps']+=1
        print('Central Gram point',str(q),len(gram),flush=True)
    for p in (16,32,64,128,256):
        N=4*p*p;lam=F(1,4)
        for r in range(1,7):
            ints=tuple(1+(i%2) for i in range(r));k=sum(ints)
            exact=c.central_singleton_coefficient(N,p,k,ints)
            limit=c.limiting_singleton_coefficient(lam,ints)
            finite.append(dict(N=N,p=p,chords=k,internals=ints,exact=str(exact),limit=str(limit),error=str(exact-limit)))
        z=c.central_G2_zero_chaos(N,p)
        ratio=F((N//2-p)*(N//2+p+1),(N//2)*(N//2+1))
        assert z==(N//2)*(ratio-1)
        assert z+2*lam==F(1-2*p,4*p*p+2)
        counts['finite_N_G2_scalar_bridges']+=1
    paths=['code/central_ordinary_sources.py','code/check_central_ordinary_sources.py',
           'code/memory_body_source_chords.py','code/memory_joint_regression.py',
           'code/slater_higher_wick_sources.py','code/slater_composite_sources.py',
           'code/slater_joint_site_polynomial.py','code/memory_size_marked_chords.py',
           'code/bilocal_multisource_chords.py']
    result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),
        environment={'python':platform.python_version(),'arithmetic':'fractions.Fraction'},
        parameters={'nilpotent_sources_max':6,'internal_counts':[0,1,2],
                    'local_slots':[2,4,6],'Gram_max_raw_slots':10,
                    'formal_points':[['1/4','1/3'],['2/3','1/4']],
                    'finite_N_p':[[4*p*p,p] for p in (16,32,64,128,256)]},
        counts=dict(counts),grams=gram,finite_singleton_rows=finite,
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        runtime_seconds=time.perf_counter()-start,
        scope='Exact source and Gaussian Fourier identities, same-copy Gram tests, and finite-N conditional singleton/scalar calibrations. General Fourier domination and fixed-degree convergence require the companion proof; no degree-uniform time-tail certification.')
    OUT.parent.mkdir(parents=True,exist_ok=True)
    with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({'status':'passed','counts':dict(counts),'seconds':result['runtime_seconds']}),flush=True)

if __name__=='__main__':main()
