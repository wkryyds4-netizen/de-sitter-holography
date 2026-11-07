"""Independent finite CAR matrices vs exact marked-site polynomials."""
from pathlib import Path
from itertools import combinations,product
from collections import defaultdict
from datetime import datetime,timezone
import hashlib,json,platform,time
import numpy as np
from slater_joint_site_polynomial import pairings,polynomial
from slater_covariant_transfer import marked_site,poly_mul
ROOT=Path(__file__).resolve().parent.parent
DIR=ROOT/'results/susskind-slater-covariant'
OUT=DIR/'site-checks.json'
if OUT.exists():raise FileExistsError(OUT)
t0=time.perf_counter()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def annihilator(n,i):
    a=np.zeros((1<<n,1<<n),dtype=np.int64)
    for x in range(1<<n):
        if x>>i&1:a[x^(1<<i),x]=(-1)**((x&((1<<i)-1)).bit_count())
    return a
fixtures=[]
for n,p in ((4,1),(4,2),(6,1)):
    cs=[annihilator(n,i) for i in range(n)];ident=np.eye(1<<n,dtype=np.int64)
    sets=list(combinations(range(n),p));labels=list(product(sets,repeat=2))
    mats=[];marked=[]
    for A,B in labels:
        X=ident.copy()
        for i in A:X=X@cs[i].T
        for i in reversed(B):X=X@cs[i]
        # Full-Fock Majorana/Clifford loss, independent of local source formula.
        DX=(2*p-n)*X+sum((c.T@X@c+c@X@c.T for c in cs),np.zeros_like(X))
        mats.append(X);marked.append(DX)
    rev=[labels.index((B,A)) for A,B in labels]
    vac=(1<<(n//2))-1;checks=negative=0
    for pi in pairings(tuple(range(4))):
        for labs in product(range(len(labels)),repeat=2):
            slots=[0]*4
            for e,(i,j) in enumerate(pi):slots[i]=rev[labs[e]];slots[j]=labs[e]
            local={0:1}
            for site in range(n):
                mask=sum(((int(site in labels[labs[e]][0])) | (int(site in labels[labs[e]][1])<<1))<<(2*e) for e in range(2))
                local=poly_mul(local,marked_site(pi,int(site<n//2),mask))
            for sources in range(16):
                v=ident[:,vac]
                for slot in reversed(range(4)):
                    v=(marked if sources>>slot&1 else mats)[slots[slot]]@v
                value=int(v[vac]);expected=local.get(sources,0)
                assert value==expected,(n,p,pi,labs,sources,value,expected)
                checks+=1;negative+=value<0
    fixtures.append(dict(N=n,p=p,comparisons=checks,negative_results=negative))
    print(json.dumps(fixtures[-1]),flush=True)

# All source coefficients of every two-label cluster through three edges.
cluster_checks=0;nonzero_source_coefficients=0
for k in (1,2,3):
    for pi in pairings(tuple(range(2*k))):
        for kinds in product('GF',repeat=k):
            totals=defaultdict(int)
            for occupied in (0,1):
                for inds in combinations(range(2*k),2):
                    mask=sum(1<<i for i in inds)
                    forbidden=any(kind=='F' and ((occupied and mask>>(2*e)&1) or (not occupied and mask>>(2*e+1)&1)) for e,kind in enumerate(kinds))
                    if forbidden:continue
                    for z,v in marked_site(pi,occupied,mask).items():totals[z]+=v
            expected={0:sum(sum(v for mask,v in polynomial(pi,kinds,o).items() if mask.bit_count()==2) for o in (0,1))}
            for e,(i,j) in enumerate(pi):
                if kinds[e]=='G':expected[1<<i]=1;expected[1<<j]=1;expected[(1<<i)|(1<<j)]=2
            for z in range(1<<(2*k)):
                assert totals.get(z,0)==expected.get(z,0),(pi,kinds,z,totals.get(z),expected.get(z))
                cluster_checks+=1;nonzero_source_coefficients+=z!=0 and totals.get(z,0)!=0

DIR.mkdir(parents=True,exist_ok=True)
inputs=[Path(__file__),ROOT/'code/slater_covariant_transfer.py',ROOT/'code/slater_joint_site_polynomial.py',ROOT/'results/susskind-slater-stationary/local-evidence.json']
result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),fixtures=fixtures,
            cluster_checks=cluster_checks,nonzero_source_coefficients=nonzero_source_coefficients,
            parameters=dict(fixtures=[[4,1],[4,2],[6,1]],finite_edges=2,all_endpoint_masks=True,cluster_edges=[1,2,3],source='2(p-half-Clifford-size)'),
            environment=dict(python=platform.python_version(),numpy=np.__version__,arithmetic='int64 with small matrices; exact integer comparisons'),
            input_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in inputs},elapsed_seconds=time.perf_counter()-t0,
            scope='Finite exact CAR/site comparisons and local source coefficient table. Arbitrary-order double scaling requires the companion normalized-defect proof.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',CAR_comparisons=sum(x['comparisons'] for x in fixtures),cluster_checks=cluster_checks,sha256=sha(OUT))))
