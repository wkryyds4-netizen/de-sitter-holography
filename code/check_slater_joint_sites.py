"""Independent finite CAR vs exact single-site polynomial and six-edge tables."""
from pathlib import Path
from itertools import combinations,product
from collections import Counter
from datetime import datetime,timezone
import hashlib,json,platform,time
from slater_joint_site_polynomial import pairings,polynomial,local_coefficient,coefficient,penalty

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-slater-joint-words/site-checks.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def act(state,P,Q):
    sign=1
    for i,creation in [(i,False) for i in Q]+[(i,True) for i in reversed(P)]:
        bit=1<<i
        if bool(state&bit)==creation:return 0
        if (state&(bit-1)).bit_count()%2:sign=-sign
        state=state|bit if creation else state&~bit
    return sign*(state+1)

def fixture(n,p):
    vac=(1<<(n//2))-1
    subsets=list(combinations(range(n),p));labels=list(product(subsets,repeat=2))
    index={x:i for i,x in enumerate(labels)};rev=[index[Q,P] for P,Q in labels]
    acts=[[act(s,P,Q) for s in range(1<<n)] for P,Q in labels]
    memberships=[[int(i in P)+2*int(i in Q) for i in range(n)] for P,Q in labels]
    allowed=[('G','F') if min(P)>=n//2 and max(Q)<n//2 else
             ('G','B') if max(P)<n//2 and min(Q)>=n//2 else ('G',) for P,Q in labels]
    rows=[];checks=Counter();signs=Counter();failures=[]
    for pi in pairings(tuple(range(6))):
        polye=[local_coefficient(pi,'GGG',0,x) for x in range(64)]
        polyo=[local_coefficient(pi,'GGG',1,x) for x in range(64)]
        totals=Counter()
        for ls in product(range(len(labels)),repeat=3):
            slots=[0]*6
            for e,(i,j) in enumerate(pi):slots[i]=rev[ls[e]];slots[j]=ls[e]
            s=vac;v=1
            for lab in reversed(slots):
                z=acts[lab][s]
                if not z:v=0;break
                v*=1 if z>0 else -1;s=abs(z)-1
            if s!=vac:v=0
            local=1
            for site in range(n):
                mask=sum(memberships[ls[e]][site]<<(2*e) for e in range(3))
                local*=(polyo if site<n//2 else polye)[mask]
                if not local:break
            checks['CAR_vs_site_assignments']+=1;signs[str(v)]+=1
            if v!=local and len(failures)<10:failures.append(dict(pi=pi,labels=ls,CAR=v,site=local))
            if v:
                for kinds in product(*(allowed[l] for l in ls)):totals[''.join(kinds)]+=v
        # Every G/F/B edge specification, all fifteen three-edge pairings.
        for letters in product('GFB',repeat=3):
            kinds=''.join(letters);expected=coefficient(n,p,pi,kinds);actual=totals[kinds]
            checks['polynomial_coefficients']+=1
            if expected!=actual:failures.append(dict(pi=pi,kinds=kinds,CAR=actual,coefficient=expected))
            rows.append(dict(pairing=pi,kinds=kinds,raw=actual))
    return dict(N=n,p=p,checks=dict(checks),signs=dict(signs),rows=rows,failures=failures)

def quadratic_tables():
    rows=[]
    for k in range(1,5):
        for pi in pairings(tuple(range(2*k))):
            for letters in product('GF',repeat=k):
                kinds=''.join(letters)
                singleton=[sum(v for mask,v in polynomial(pi,kinds,o).items() if mask.bit_count()==1) for o in (0,1)]
                doubles=[sum(v for mask,v in polynomial(pi,kinds,o).items() if mask.bit_count()==2) for o in (0,1)]
                assert singleton==[k,k]
                expected=kinds.count('G')+k*(k-1)-2*penalty(pi,kinds)
                assert sum(doubles)==expected,(pi,kinds,doubles,expected)
                rows.append(dict(pairing=pi,kinds=kinds,double_sum=sum(doubles),penalty=penalty(pi,kinds)))
    return rows

def main():
    if OUT.exists():raise FileExistsError(OUT)
    t=time.perf_counter();fixtures=[]
    for n,p in [(4,1),(6,1),(4,2)]:
        x=fixture(n,p);fixtures.append(x)
        print(json.dumps(dict(N=n,p=p,checks=x['checks'],failures=len(x['failures']),elapsed=time.perf_counter()-t)),flush=True)
    rows=quadratic_tables()
    files=[Path(__file__),ROOT/'code/slater_joint_site_polynomial.py',ROOT/'results/susskind-slater-joint-fourth/local-evidence.json']
    failures=[f for x in fixtures for f in x['failures']]
    out=dict(status='passed' if not failures else 'failed',created_utc=datetime.now(timezone.utc).isoformat(),
        parameters=dict(fixtures=[(4,1),(6,1),(4,2)],edges=3,all_GFB_types=True,quadratic_edges=[1,2,3,4]),
        fixtures=fixtures,quadratic_rows=rows,quadratic_checks=len(rows),failures=failures,
        environment=dict(python=platform.python_version(),dependencies='standard library'),
        input_sha256={x.relative_to(ROOT).as_posix():sha(x) for x in files},elapsed_seconds=time.perf_counter()-t,
        scope='Exact integer comparison of microscopic CAR with site polynomial for all three-edge label assignments in three fixtures; exact coefficient extraction for every G/F/B specification. Quadratic cluster table checked through four edges. Does not alone prove a double-scaled limit.')
    OUT.parent.mkdir(parents=True,exist_ok=True)
    with OUT.open('x',encoding='utf-8') as f:json.dump(out,f,indent=2);f.write('\n')
    print(json.dumps(dict(status=out['status'],quadratic_checks=len(rows),sha256=sha(OUT))))
    assert not failures
if __name__=='__main__':main()
