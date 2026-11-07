"""Complete exact same-J Slater mixed fourth moments, independent CAR enumeration.

H, +, - denote unnormalized original H, pure excitation R^dagger and R.
Every Gaussian pair has unit covariance. No external scientific driver runs.
"""
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, exp
from pathlib import Path
import hashlib, json, platform, time

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-slater-joint-fourth/exact-checks.json'
FIXTURES=[(4,1),(6,1),(6,2),(8,2)]
PAIRINGS=[((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2))]
PARAMETERS=dict(fixtures=FIXTURES, alphabet='H+-', order='written left to right, rightmost acts first',
    covariance='E J_AB J_CD = delta_AD delta_BC', filling='m=N/2',
    arithmetic='integer CAR and Fraction', words='all 3^4 fourth-order words and all 3^2 second-order words',
    asymptotic_p=[2,4,8,16,32,64], asymptotic_lambda='1/2')

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def choose(n,k):return comb(n,k) if 0<=k<=n else 0

def act(state,A,B):
    sign=1
    for i,creation in [(i,False) for i in B]+[(i,True) for i in reversed(A)]:
        bit=1<<i
        if bool(state&bit)==creation:return 0
        if (state&(bit-1)).bit_count()%2:sign=-sign
        state=state|bit if creation else state&~bit
    return sign*(state+1)

def pair_shape(pairing):
    return 'disjoint' if pairing==PAIRINGS[0] else 'crossed' if pairing==PAIRINGS[1] else 'nested'

def edge_type(word,edge):
    i,j=edge
    if word[i]==word[j]=='H':return 'G'
    if word[i]=='+' or word[j]=='-':
        if word[i]=='-' or word[j]=='+':return 'forbidden'
        return 'backward'
    return 'forward'

def formula(word,pairing,L,D,R,T,U,K):
    types=[edge_type(word,e) for e in pairing]
    if 'forbidden' in types:return 0
    shape=pair_shape(pairing)
    if types==['G','G']:return K if shape=='crossed' else D*D
    if 'G' not in types:
        if types==['forward','forward']:return L*L if shape=='disjoint' else L*L*R
        return L if shape=='nested' and types==['forward','backward'] else 0
    gs=types.index('G'); ss=1-gs
    if types[ss]=='backward':return L if shape=='nested' and gs==0 else 0
    if shape=='disjoint':return D*L
    if shape=='crossed':return L*T
    return U if gs==0 else D*L

def scalar_data(n,p):
    m=a=n//2
    L=comb(m,p)*comb(a,p);D=comb(m,p)*comb(a+p,p)
    R=F(choose(m-p,p)*choose(a-p,p),L)
    T=choose(m-p,p)*comb(a,p)
    W=[comb(m,p)*choose(p,k)*choose(a,k) for k in range(p+1)]
    U=sum(W[k]*choose(a-k,p)*choose(m-k,p) for k in range(p+1))
    return L,D,R,T,U,W

def fixture(n,p):
    m=n//2;vac=(1<<m)-1
    subsets=list(combinations(range(n),p))
    labels=list(product(subsets,repeat=2));index={x:i for i,x in enumerate(labels)}
    reverse=[index[B,A] for A,B in labels]
    kinds=[('H','+') if min(A)>=m and max(B)<m else
           ('H','-') if max(A)<m and min(B)>=m else ('H',) for A,B in labels]
    states=[s for s in range(1<<n) if s.bit_count()==m]
    table=[[act(s,A,B) if s.bit_count()==m else 0 for s in range(1<<n)] for A,B in labels]
    words=[''.join(w) for w in product('H+-',repeat=4)]
    sums={w:[0,0,0] for w in words};counts=Counter();signed=Counter()
    for pi,pairing in enumerate(PAIRINGS):
        for ia in range(len(labels)):
            for ib in range(len(labels)):
                slots=[0]*4
                for (i,j),lab in zip(pairing,(ia,ib)):
                    slots[i]=lab;slots[j]=reverse[lab]
                state=vac;sign=1
                for slot in reversed(slots):
                    out=table[slot][state]
                    if not out:sign=0;break
                    if out<0:sign=-sign
                    state=abs(out)-1
                value=sign if state==vac else 0
                counts['Gaussian_pair_label_assignments']+=1
                signed[str(value)]+=1
                if value:
                    for letters in product(*(kinds[i] for i in slots)):
                        sums[''.join(letters)][pi]+=value
                        counts['nonzero_word_assignments']+=1
    L,D,R,T,U,W=scalar_data(n,p)
    K=sums['HHHH'][1]
    rows=[]
    for word in words:
        predicted=[formula(word,pairing,L,D,R,T,U,K) for pairing in PAIRINGS]
        assert list(map(F,sums[word]))==predicted,(n,p,word,sums[word],predicted)
        rows.append(dict(word=word,pairing_contributions=sums[word],total=sum(sums[word]),
                         normalized_denominator=dict(D_power=word.count('H'),L_power=4-word.count('H'),outer_power='1/2')))
    second={''.join(w):0 for w in product('H+-',repeat=2)}
    for lab in range(len(labels)):
        out=table[reverse[lab]][vac]
        if not out:continue
        back=table[lab][abs(out)-1]
        value=(1 if out>0 else -1)*(1 if back>0 else -1) if abs(back)-1==vac else 0
        for w in product(kinds[lab],kinds[reverse[lab]]):second[''.join(w)]+=value
    assert second=={'HH':D,'H+':L,'H-':0,'+H':0,'++':0,'+-':0,'-H':L,'-+':L,'--':0}
    # Direct independent two-label action checks crossed covariance on one p-hop.
    alpha=index[subsets[-1],subsets[0]]
    X=table[alpha]
    covariance=[]
    for state in states:
        values=Counter()
        for lab in range(len(labels)):
            o=table[reverse[lab]][state]
            if not o:continue
            z=X[abs(o)-1]
            if not z:continue
            y=table[lab][abs(z)-1]
            if not y:continue
            values[abs(y)-1]+=(1 if o*z*y>0 else -1)
        values={s:v for s,v in values.items() if v}
        expected={} if not X[state] or not T else {abs(X[state])-1:T*(1 if X[state]>0 else -1)}
        assert values==expected,(n,p,state,values,expected)
        covariance.append(dict(state=state,nonzero=values))
    return dict(N=n,p=p,L=L,D=D,ratio_R=str(R),crossed_eigenvalue_T=T,
        nested_U=U,H_crossed_K=K,H_first_step_weights=W,counts=dict(counts),
        raw_CAR_sign_counts=dict(signed),second_moments=second,words=rows,
        covariance_columns_verified=len(states),covariance_columns=covariance)

def asymptotic():
    rows=[]
    for p in PARAMETERS['asymptotic_p']:
        n=2*p*p;L,D,R,T,U,W=scalar_data(n,p)
        gamma2=F(L,D);cross=F(T,D);nest=F(U,D*L)
        rows.append(dict(N=n,p=p,gamma_squared=str(gamma2),R=str(R),cross=str(cross),nest=str(nest),
            inverse_D=str(F(1,D)),decimal=dict(gamma_squared=float(gamma2),R=float(R),cross=float(cross),nest=float(nest)),
            limiting=dict(gamma_squared=exp(-1),R=exp(-2),cross=exp(-2),nest=exp(-2))))
    return rows

def main():
    if OUT.exists():raise FileExistsError(OUT)
    start=time.perf_counter();fixtures=[]
    for n,p in FIXTURES:
        f=fixture(n,p);fixtures.append(f)
        print(json.dumps(dict(N=n,p=p,counts=f['counts'],elapsed=time.perf_counter()-start)),flush=True)
    rows=asymptotic()
    source=Path(__file__)
    parents=['results/susskind-slater-ladder/local-evidence.json',
             'results/susskind-bilocal-four-point/local-evidence.json']
    result=dict(status='passed',created_utc=datetime.now(timezone.utc).isoformat(),
        parameters=PARAMETERS,parameter_sha256=digest(PARAMETERS),fixtures=fixtures,asymptotic=rows,
        counts=dict(finite_systems=len(fixtures),fourth_words=81*len(fixtures),pairing_formulas=243*len(fixtures),
            second_words=9*len(fixtures),label_assignments=sum(f['counts']['Gaussian_pair_label_assignments'] for f in fixtures),
            covariance_columns=sum(f['covariance_columns_verified'] for f in fixtures)),
        environment=dict(python=platform.python_version(),platform=platform.platform(),dependencies='standard library'),
        input_sha256={source.relative_to(ROOT).as_posix():sha(source),**{p:sha(ROOT/p) for p in parents}},
        scope='Exact all fourth/second words for finite fixtures. H-only crossed K is enumerated, not independently predicted here. Other pairing formulas and covariance eigenvalue are independently checked. Large-p scalar rows diagnose analytic limits, not interval certificates.',
        elapsed_seconds=time.perf_counter()-start,failures=[])
    OUT.parent.mkdir(parents=True,exist_ok=True)
    with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps(dict(status='passed',output=OUT.relative_to(ROOT).as_posix(),counts=result['counts'],sha256=sha(OUT))))

if __name__=='__main__':main()
