"""Fresh exact one-selected-state checks; exclusive scientific output."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from datetime import datetime,timezone
from collections import Counter
import hashlib,json,platform,time
import slater_composite_sources as source
import slater_composite_transfer as transfer
import slater_selected_covariant_state as state

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-slater-selected-state/exact-checks.json'
if OUT.exists():raise FileExistsError(OUT)
OUT.parent.mkdir(parents=True,exist_ok=True)
started=time.perf_counter();points=[(F(1,2),F(1,3)),(F(2,3),F(1,4))]
words=[''.join(t) for n in range(7) for t in product('HWJ',repeat=n)]
count=0
for gamma,lam in points:
    for word in words:
        actual=transfer.moment('-'+word+'+',gamma,lam)
        predicted=state.selected(word,gamma,lam)
        assert actual==predicted,(word,gamma,lam,actual,predicted);count+=1
    print('all words',str(gamma),count,flush=True)

# The proof's stronger diagram-level statement: an outer selected pair
# changes no quadratic source and adds no inner-F nesting.
diagrams=direct=attached=0
for n in range(5):
    for letters in product('HWJ',repeat=n):
        word='-'+''.join(letters)+'+';slots,blocks,types=source.layout(word)
        if len(slots)%2:continue
        last=len(slots)-1
        for pi in source.pairings(tuple(range(len(slots)))):
            if any(a<=i<j<b for a,b,_ in blocks for i,j in pi):continue
            kinds=tuple('F' if i==0 or j==last else 'G' for i,j in pi)
            selected=source.source_exponent(pi,kinds,blocks)
            ordinary=source.source_exponent(pi,tuple('G' for _ in pi),blocks)
            predicted=Counter(ordinary)
            for b,(a,_,kind) in enumerate(blocks):
                if kind=='L' and any(a in edge and k=='F' for edge,k in zip(pi,kinds)):
                    predicted[1<<b]-=1
            predicted={k:v for k,v in predicted.items() if v}
            assert selected==predicted,(word,pi,selected,predicted)
            diagrams+=1
            if (0,last) in pi:direct+=1
            else:attached+=1
        source.local.cache_clear();source.source_exponent.cache_clear()
print('boundary source diagrams',diagrams,flush=True)

reduced=0;curve_checks=[];nonstationary=[]
for gamma,lam in points:
    for n in range(9):
        for tup in product('HJ',repeat=n):
            word=''.join(tup)
            assert state.selected(word,gamma,lam)==state.reduced_hj(word,gamma,lam),word
            reduced+=1
    for n in range(13):
        w='J'*n
        assert state.selected(w,gamma,lam)==state.trace(w,gamma,lam),n
    q=gamma**4;v=4*lam+16*lam*lam*q/(1+q)
    comm=(transfer.moment('-JHHJ+',gamma,lam)-transfer.moment('-JHJH+',gamma,lam)
          -transfer.moment('-HJHJ+',gamma,lam)+transfer.moment('-HJJH+',gamma,lam))/v
    assert comm==state.curvature_commutator(gamma)
    curve_checks.append(dict(gamma=str(gamma),lam=str(lam),commutator_norm=str(comm)))
    rho={'HH':gamma**2,'HW':-gamma**2,'WH':-gamma**2,'WW':gamma**2}
    poly=Counter()
    for w,c in rho.items():poly['H'+w]+=c;poly[w+'H']-=c
    norm=sum(x*y*state.trace(a[::-1]+b,gamma,lam) for a,x in poly.items() for b,y in poly.items())
    assert norm==state.density_nonstationarity(gamma,lam)>0
    nonstationary.append(dict(gamma=str(gamma),lam=str(lam),density_commutator_norm=str(norm)))

paths=['code/slater_selected_covariant_state.py','code/check_slater_selected_covariant_state.py',
       'code/slater_composite_transfer.py','code/slater_composite_sources.py',
       'code/slater_joint_site_polynomial.py','code/bilocal_multisource_chords.py',
       'code/memory_body_source_chords.py','code/singlet_coefficient.py']
result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),
    environment=dict(python=platform.python_version(),arithmetic='fractions.Fraction'),
    parameters=dict(points=[[str(g),str(l)] for g,l in points],all_words_max_length=6,
        boundary_source_internal_max_length=4,reduced_hj_max_length=8,pure_j_max_power=12,
        note='Independent rational points check identities; positivity theorem uses physical gamma=exp(-lambda).'),
    density_words=len(words),density_comparisons=count,boundary_source_diagrams=diagrams,
    direct_boundary_pairs=direct,attached_boundary_pairs=attached,
    reduced_hj_comparisons=reduced,pure_curvature_comparisons=26,
    curvature_commutators=curve_checks,nonstationarity=nonstationary,
    input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
    runtime_seconds=time.perf_counter()-started,
    scope='Finite exact tests of the all-word boundary-source proof, reduced state and dynamic coefficients; not a numerical proof of arbitrary degree or a finite-N pinching convergence rate.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({k:result[k] for k in ['status','density_comparisons','boundary_source_diagrams','reduced_hj_comparisons','runtime_seconds']}),flush=True)
