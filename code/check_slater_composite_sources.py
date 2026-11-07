"""Independent finite CAR/source and same-J composite-word checks."""
from pathlib import Path
from itertools import combinations,product
from fractions import Fraction as F
from functools import lru_cache
import hashlib,json,platform,time
import numpy as np
import sympy as s
from slater_composite_sources import local,layout,source_exponent,raw_moment,symbolic,evaluate
from slater_joint_site_polynomial import pairings
from slater_covariant_transfer import moment as old_transfer
from memory_body_source_chords import mul,exponent as old_exponent
from bilocal_multisource_chords import raw_moment as old_raw,evaluate as old_eval

ROOT=Path(__file__).resolve().parent.parent
DIR=ROOT/'results/susskind-slater-composite';OUT=DIR/'exact-checks.json'
if OUT.exists():raise FileExistsError(OUT)
t0=time.perf_counter();finite=[]
for N,p in [(2,1),(4,1),(4,2)]:
    m=N//2;dim=1<<N;vac=(1<<m)-1;I=np.eye(dim,dtype=np.int64)
    cs=[]
    for i in range(N):
        c=np.zeros_like(I)
        for st in range(dim):
            if st>>i&1:c[st^(1<<i),st]=(-1)**((st&((1<<i)-1)).bit_count())
        cs.append(c)
    Q=sum((c.T@c for c in cs),np.zeros_like(I))
    subsets=list(combinations(range(N),p));labs=list(product(subsets,repeat=2))
    mats=[]
    for A,B in labs:
        X=I.copy()
        for i in A:X=X@cs[i].T
        for i in reversed(B):X=X@cs[i]
        mats.append(X)
    rev=[labs.index((B,A)) for A,B in labs]
    @lru_cache(None)
    def factors(seq,kind):
        X=I.copy()
        for z in seq:X=X@mats[z]
        remove=sum((c.T@X@c for c in cs),np.zeros_like(I))
        if kind=='D':marked=2*p*len(seq)*X-2*Q@X+2*remove
        else:marked=(2*p*len(seq)-N)*X+remove+sum((c@X@c.T for c in cs),np.zeros_like(I))
        return X,marked
    layouts=[((0,2,'D'),(2,4,'D')),((0,2,'D'),(2,3,'L'),(3,4,'L')),((1,3,'D'),)]
    count=negative=selected=0;maximum=0
    for pi in pairings(tuple(range(4))):
        for labels in product(range(len(labs)),repeat=2):
            slots=[0]*4
            for e,(a,b) in enumerate(pi):slots[a]=rev[labels[e]];slots[b]=labels[e]
            masks=[sum((int(site in labs[z][0])|(int(site in labs[z][1])<<1))<<(2*e)
                       for e,z in enumerate(labels)) for site in range(N)]
            for blocks in layouts:
                starts={a:(i,b,kind) for i,(a,b,kind) in enumerate(blocks)}
                prepared=[];pos=0
                while pos<4:
                    mark,end,kind=starts.get(pos,(None,pos+1,None))
                    if mark is None:prepared.append((None,mats[slots[pos]],None))
                    else:
                        X,Y=factors(tuple(slots[pos:end]),kind);prepared.append((mark,X,Y))
                        maximum=max(maximum,int(np.max(np.abs(Y))))
                    pos=end
                exact={}
                for mask in range(1<<len(blocks)):
                    v=I[:,vac]
                    for mark,X,Y in reversed(prepared):v=(Y if mark is not None and mask>>mark&1 else X)@v
                    exact[mask]=int(v[vac]);maximum=max(maximum,int(np.max(np.abs(v))))
                for kinds in product('GFB',repeat=2):
                    allowed=all(kind=='G' or (all(i>=m for i in labs[z][0]) and all(i<m for i in labs[z][1]) if kind=='F'
                                else all(i<m for i in labs[z][0]) and all(i>=m for i in labs[z][1])) for kind,z in zip(kinds,labels))
                    value={0:F(1)}
                    for site,mask in enumerate(masks):value=mul(value,local(pi,kinds,blocks,int(site<m),mask))
                    for mark,(a,b,kind) in enumerate(blocks):value=mul(value,{0:F(1),1<<mark:F(p*(b-a))})
                    for mask in range(1<<len(blocks)):
                        observed=exact[mask] if allowed else 0
                        predicted=value.get(mask,0)*2**mask.bit_count()
                        assert observed==predicted,(N,p,pi,labels,kinds,blocks,mask,observed,predicted)
                        count+=1;negative+=observed<0;selected+=observed!=0 and kinds!=('G','G')
    assert maximum<2**40
    finite.append(dict(N=N,p=p,comparisons=count,negative=negative,nonzero_selected=selected,max_abs_integer=maximum))
    print(json.dumps(finite[-1]),flush=True)

allG=allGwords=oldwords=0
for size in range(1,5):
    for letters in product('HWJ',repeat=size):
        word=''.join(letters);slots,blocks,types=layout(word)
        if len(slots)%2 or len(slots)>6:continue
        for pi in pairings(tuple(range(len(slots)))):
            if any(a<=i<j<b for a,b,_ in blocks for i,j in pi):continue
            assert source_exponent(pi,('G',)*len(pi),blocks)==old_exponent(len(slots),pi,blocks)
            allG+=1
        assert evaluate(word,F(1,2),F(1,4))==old_eval(old_raw(word)[0],F(1,16),F(1,4),word.count('J'))
        allGwords+=1
for size in range(1,5):
    for letters in product('HW+-',repeat=size):
        word=''.join(letters)
        assert evaluate(word,F(1,2),F(1,4))==old_transfer(word.replace('W','U'),F(1,2),F(1,2))*F(1,2)**word.count('W')
        oldwords+=1

g,l=s.symbols('gamma lambda',positive=True);q=g**4;c=2*l*(1-q)/(1+q);v=4*l+16*l*l*q/(1+q)
expected={'JJ':v,'JHH':0,'J++':-2*c*g*g*q,'JH+':0,'J+H':-2*c*g*q,
          'JW+':0,'J+W':0,'-J+':0,'--J':-2*c*g*g*q,'-JJ+':v}
symbols={}
for word,want in expected.items():
    actual=s.factor(symbolic(word,g,l));assert s.factor(actual-want)==0
    symbols[word]=str(actual)

preps=['','H','W','+','HH','HW','WH','WW','+H','H+','+W','W+','++','J']
def adj(word):return ''.join({'H':'H','W':'W','J':'J','+':'-','-':'+'}[a] for a in reversed(word))
grams=[]
for gamma,lam in [(F(1,2),F(1,4)),(F(2,3),F(1,9))]:
    A=s.Matrix([[evaluate(adj(a)+b,gamma,lam) for b in preps] for a in preps])
    assert A==A.T
    _,diag=A.LDLdecomposition(hermitian=False);pivots=list(diag.diagonal());assert all(x>0 for x in pivots)
    grams.append(dict(gamma=str(gamma),lambda_value=str(lam),rank=len(preps),pivots=[str(x) for x in pivots],
                      scope='Exact formal rational parameters, not points on gamma=exp(-lambda).'))

inputs=['code/check_slater_composite_sources.py','code/slater_composite_sources.py','code/slater_joint_site_polynomial.py',
        'code/slater_covariant_transfer.py','code/memory_body_source_chords.py','code/memory_size_marked_chords.py',
        'code/bilocal_multisource_chords.py','code/singlet_multimark.py','code/singlet_coefficient.py',
        'code/complex_syk_dictionary.py','results/susskind-slater-polynomial-density/local-evidence.json']
result=dict(status='passed',scope='Finite CAR identities, source reduction and selected mixed-word diagnostics. General limits and pinching require the companion proofs.',
    finite_CAR=finite,central_source_comparisons=allG,central_words=allGwords,old_selected_words=oldwords,
    symbolic_selected_moments=symbols,Gram_preparations=preps,Gram_checks=grams,
    parameters=dict(CAR_systems=[[2,1],[4,1],[4,2]],local_sources='D two-slot, L one-slot',
                    Gaussian_edges=2,all_G_max_slots=6,old_selected_max_letters=4),
    environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__,arithmetic='bounded int64, Fraction and symbolic rational functions'),
    input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},elapsed_seconds=time.perf_counter()-t0)
DIR.mkdir(exist_ok=True)
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',CAR=sum(r['comparisons'] for r in finite),allG=allG,oldwords=oldwords,
                     sha256=hashlib.sha256(OUT.read_bytes()).hexdigest())))
