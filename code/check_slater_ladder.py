"""Exact CAR/Wick tests for a state-adapted complex-SYK creation block.

Independent of the previous scientific drivers. Finite backward Wick pairs
are retained, not set to zero. Creates a new output exactly once.
"""
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import comb, factorial, exp
from pathlib import Path
import hashlib, json, platform, time

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT/'results/susskind-slater-ladder/exact-checks.json'
FIXTURES = [(4,1,4),(6,1,3),(6,2,3),(8,2,3)]
PARAMETERS = dict(fixtures_N_p_maximum_pairs=FIXTURES, filling='m=N/2',
    normalization='unit Gaussian variance, sigma_squared=L=binom(m,p)binom(N-m,p)',
    CAR='{c_i,c_j^dagger}=delta_ij', creators='ascending operator order',
    word_order='left to right product; rightmost acts first; + creates a particle-hole block',
    expected_maximum_label_assignments=1600000,
    asymptotic_parameters=dict(lambda_N='1/2',p_values=[2,4,8,16,32,64],maximum_r=8),
    scope='Finite exact CAR and Wick identities; large-p points are diagnostics, not proof of convergence')

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def digest(x): return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def choose(n,k): return comb(n,k) if 0<=k<=n else 0

def act(state, sites, creation):
    sign=1
    for i in sites:
        bit=1<<i
        if bool(state&bit)==creation: return None
        if (state&(bit-1)).bit_count()%2: sign=-sign
        state=state|bit if creation else state&~bit
    return state,sign

def monomial(state,A,B):
    first=act(state,B,False)
    if first is None:return None
    last=act(first[0],tuple(reversed(A)),True)
    return None if last is None else (last[0],first[1]*last[1])

def dyck(word):
    height=0;births=[];multiplicity=1
    for char in reversed(word):
        if char=='+':births.append(height);height+=1
        else:
            if height==0:return None
            multiplicity*=height;height-=1
    return (births,multiplicity) if height==0 else None

def fixture(n,p,maxk):
    m=n//2;vac=(1<<m)-1
    labels=[(a,b) for a in combinations(range(m,n),p) for b in combinations(range(m),p)]
    masks=[sum(1<<i for i in a+b) for a,b in labels];L=len(labels)
    # Table every possible occupation input; this is raw CAR, not a Wick rule.
    actions={}
    for index,(a,b) in enumerate(labels):
        for char,A,B in [('+',a,b),('-',b,a)]:
            actions[index,char]=[monomial(s,A,B) for s in range(1<<n)]
    rows=[];counts=Counter();witness=None
    for k in range(1,maxk+1):
        for positions in combinations(range(2*k),k):
            plus=set(positions);word=''.join('+' if j in plus else '-' for j in range(2*k))
            path=dyck(word)
            if path is None:
                counts['non_Dyck_words_zero_by_right_action']+=1
                continue
            births,mult=path;minus=[i for i in range(2*k) if i not in plus]
            forward_sum=backward_sum=0;fwd_diagrams=bwd_diagrams=0
            for paired_minus in permutations(minus):
                edges=list(zip(sorted(plus),paired_minus))
                forward=all(a<b for b,a in edges)
                slots={}
                for lab,(b,a) in enumerate(edges):slots[a]=lab;slots[b]=lab
                overlaps=[(i,j) for i in range(k) for j in range(i+1,k)
                    if max(edges[i][1],edges[j][1])<min(edges[i][0],edges[j][0])]
                subtotal=0
                for assignment in product(range(L),repeat=k):
                    state=vac;sign=1
                    for pos in range(2*k-1,-1,-1):
                        out=actions[assignment[slots[pos]],word[pos]][state]
                        if out is None:sign=0;break
                        state,phase=out;sign*=phase
                    value=sign if state==vac else 0
                    subtotal+=value;counts['Gaussian_label_assignments']+=1
                    if forward:
                        predicted=int(all(not masks[assignment[i]]&masks[assignment[j]] for i,j in overlaps))
                        assert value==predicted,(n,p,word,edges,assignment,value,predicted)
                        counts['forward_CAR_interval_comparisons']+=1
                    elif value and witness is None:
                        witness=dict(word=word,pairing=edges,label_assignment=assignment,CAR_amplitude=value)
                if forward:forward_sum+=subtotal;fwd_diagrams+=1
                else:backward_sum+=subtotal;bwd_diagrams+=1
            expected=F(mult)
            for h in births:expected*=F(choose(n-m-h*p,p)*choose(m-h*p,p),L)
            assert F(forward_sum,L**k)==expected,(n,p,word)
            assert fwd_diagrams==mult
            full=F(forward_sum+backward_sum,L**k)
            assert full>=0
            counts['word_moments']+=1
            counts['nonzero_backward_word_corrections']+=int(backward_sum!=0)
            rows.append(dict(word=word,creation_heights=births,forward_pairings=fwd_diagrams,
                backward_pairings=bwd_diagrams,forward_moment=str(expected),
                backward_correction=str(F(backward_sum,L**k)),full_moment=str(full),
                limiting_shift_moment=f'{mult} * q^{sum(births)}'))
    # Resolve the complete original p-body H acting on the fixed Slater vector.
    weights=Counter();total=0
    allsets=list(combinations(range(n),p))
    for A in allsets:
        for B in allsets:
            out=monomial(vac,A,B)
            if out:
                k=(out[0]>>m).bit_count();weights[k]+=out[1]**2;total+=1
    cm=choose(m,p)*choose(n-m+p,p)
    assert total==cm
    for k in range(p+1):
        assert weights[k]==choose(m,p)*choose(p,k)*choose(n-m,k)
    assert weights[p]==L
    # E H=0 in the joint disorder/state GNS; all other creation levels have
    # different excitation numbers. Hence the cyclic-space projection is V+Ω.
    first_step=dict(C_m_at_w_1=cm,sigma_squared_at_w_1=L,
        weights_by_created_pairs={str(k):weights[k] for k in range(p+1)},
        captured_weight=str(F(L,cm)),leakage_fraction=str(F(cm-L,cm)),
        H_Vdagger_overlap_at_w_1=L,
        normalized_overlap_squared=str(F(L,cm)))
    return dict(N=n,p=p,m=m,L=L,counts=dict(counts),records=rows,
                nonzero_backward_pair_witness=witness,full_H_first_step=first_step)

def asymptotic():
    rows=[]
    for p in PARAMETERS['asymptotic_parameters']['p_values']:
        n=2*p*p;m=n//2;den=comb(n-m+p,p)
        distribution=[F(choose(p,r)*choose(n-m,p-r),den) for r in range(p+1)]
        assert sum(distribution)==1
        factorial_moments=[]
        for t in range(1,min(4,p)+1):
            lhs=sum(F(factorial(r),factorial(r-t))*distribution[r] for r in range(t,p+1))
            rhs=F((factorial(p)//factorial(p-t))**2,factorial(n-m+p)//factorial(n-m+p-t))
            assert lhs==rhs
            factorial_moments.append(str(lhs))
        rows.append(dict(N=n,p=p,lambda_N='1/2',captured_weight=str(distribution[0]),
            captured_weight_decimal=float(distribution[0]),limit_captured_weight=exp(-1),
            saved_r_probabilities=[str(v) for v in distribution[:9]],
            exact_mean_deficit=str(sum(r*v for r,v in enumerate(distribution))),
            factorial_moments=factorial_moments,
            total_variation_to_Poisson1=0.5*(sum(abs(float(v)-exp(-1)/factorial(r)) for r,v in enumerate(distribution))+
                1-sum(exp(-1)/factorial(r) for r in range(p+1)))))
    return rows

def linear_diagnostic():
    rows=[]
    # This is a new symbolic consequence of the sealed mixed-Q word law.
    for q in [F(1,4),F(1,2),F(3,4)]:
        Q=[[q,q,q*q],[q,q,q*q],[q*q,q*q,q**4]]
        for a,b,c in product([F(-1),F(0),F(1),F(2)],repeat=3):
            if not(a or b or c):continue
            coeff=[a,b,c];variance=sum(v*v for v in coeff)
            moment=F(0)
            for i,j,k,l in product(range(3),repeat=4):
                wick=int(i==j and k==l)+int(i==l and j==k)+Q[i][j]*int(i==k and j==l)
                moment+=coeff[i]*coeff[j]*coeff[k]*coeff[l]*wick
            A=a*a+b*b;B=c*c
            cumulant=-(1-q)*A*A-2*(1-q*q)*A*B-(1-q**4)*B*B
            assert moment-3*variance*variance==cumulant<0
            assert (q**4-1)*variance**2<=cumulant<=(q-1)*variance**2
            rows.append(dict(q=str(q),coefficients=list(map(str,coeff)),variance=str(variance),
                fourth_moment=str(moment),fourth_cumulant=str(cumulant)))
    return rows

def main():
    if OUT.exists():raise FileExistsError(OUT)
    estimated=sum(sum(comb(2*k,k)//(k+1)*factorial(k)*(comb(n//2,p)**2)**k
        for k in range(1,kmax+1)) for n,p,kmax in FIXTURES)
    assert estimated<PARAMETERS['expected_maximum_label_assignments'],estimated
    started=time.perf_counter();fixtures=[]
    for cfg in FIXTURES:
        row=fixture(*cfg);fixtures.append(row)
        print(json.dumps(dict(fixture=cfg,counts=row['counts'],elapsed=time.perf_counter()-started)),flush=True)
    asym=asymptotic();linear=linear_diagnostic()
    count=Counter()
    for f in fixtures:count.update(f['counts'])
    assert count['Gaussian_label_assignments']==estimated
    files=['code/check_slater_ladder.py','results/susskind-bilocal-fourth-tail/local-evidence.json',
           'results/susskind-bilocal-four-point/local-evidence.json']
    output=dict(status='passed',created_utc=datetime.now(timezone.utc).isoformat(),
        parameters=PARAMETERS,parameter_sha256=digest(PARAMETERS),
        counts=dict(count,finite_systems=len(fixtures),asymptotic_rows=len(asym),linear_cumulant_checks=len(linear)),
        fixtures=fixtures,hypergeometric_asymptotic=asym,linear_cumulants=linear,
        input_sha256={f:sha(ROOT/f) for f in files},
        environment=dict(python=platform.python_version(),platform=platform.platform(),
            arithmetic='exact integers/Fraction; asymptotic decimal summaries use Python float'),
        elapsed_seconds=time.perf_counter()-started,
        scope=PARAMETERS['scope'],failures=[])
    OUT.parent.mkdir(parents=True,exist_ok=True)
    with OUT.open('x',encoding='utf-8') as f:json.dump(output,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps(dict(status='passed',counts=output['counts'],sha256=sha(OUT))))

if __name__=='__main__':main()
