"""Exact checks for polynomial Slater orbit densities; no old drivers imported."""
from pathlib import Path
from itertools import combinations, product
from collections import defaultdict, Counter
from math import comb, factorial, prod
from fractions import Fraction
import hashlib, json, platform
import sympy as s
from slater_covariant_transfer import moment

ROOT=Path(__file__).resolve().parent.parent
DIR=ROOT/'results/susskind-slater-polynomial-density'
OUT=DIR/'exact-checks.json'
if OUT.exists(): raise FileExistsError(OUT)

def hop(state,A,B):
    sign=1
    for i in B:
        if not state>>i&1: return None
        sign*=(-1)**((state&((1<<i)-1)).bit_count());state^=1<<i
    for i in reversed(A):
        if state>>i&1: return None
        sign*=(-1)**((state&((1<<i)-1)).bit_count());state|=1<<i
    return state,sign

rows=[];seed_checks=partial_checks=0
for N,m,p in [(4,2,1),(6,3,1),(4,2,2),(6,3,2),(8,4,2)]:
    sets=list(combinations(range(N),p)); labels=list(product(sets,repeat=2))
    selected=[i for i,(A,B) in enumerate(labels) if min(A)>=m and max(B)<m]
    rev={i:labels.index((B,A)) for i,(A,B) in enumerate(labels)}
    O=(1<<m)-1;w=tuple(int(i<m) for i in range(N))
    weight=[tuple(int(t in B)-int(t in A) for t in range(N)) for A,B in labels]
    L=comb(m,p)*comb(N-m,p);C=comb(m,p)*comb(N-m+p,p)
    for prep in ['MdagH','Mdag2']:
        raw={};contraction=defaultdict(int)
        for j in (range(len(labels)) if prep=='MdagH' else selected):
            first=hop(O,*labels[j])
            if first is None:continue
            for i in selected:
                second=hop(first[0],*labels[i])
                if second is not None:
                    t,sign=second;raw[i,j,t]=sign*first[1]
                    if rev[i]==j:contraction[t]+=sign*first[1]
        assert all(v==0 for v in contraction.values())
        sym=defaultdict(Fraction)
        for (i,j,t),v in raw.items():
            sym[i,j,t]+=Fraction(v,2);sym[j,i,t]+=Fraction(v,2)
        sym={key:v for key,v in sym.items() if v}
        for (i,j,t),v in sym.items():
            actual=tuple(weight[i][a]+weight[j][a]+int(t>>a&1) for a in range(N))
            assert actual==w;seed_checks+=1
            for removed in [(),(i,),(j,),(i,j)]:
                nu=tuple(sum(weight[z][a] for z in removed) for a in range(N))
                z=tuple(w[a]-nu[a] for a in range(N))
                assert sum(z)==m and sum(abs(z[a]-w[a]) for a in range(N))<=2*len(removed)*p
                partial_checks+=1
        denom=L*C if prep=='MdagH' else L*L
        tensor_norm=sum(v*v for v in sym.values())
        independent=sum(v*v+v*raw.get((j,i,t),0) for (i,j,t),v in raw.items())
        assert independent==2*tensor_norm
        rows.append(dict(N=N,m=m,p=p,preparation=prep,raw_terms=len(raw),symmetric_terms=len(sym),
                         chaos0_norm_squared='0',unnormalized_symmetric_norm=str(tensor_norm),
                         normalized_chaos2_norm=str(Fraction(independent,denom))))

# Rational bounds used in the factorial lemma, including negative shifts.
factorial_checks=0
for m in range(2,201):
    for delta in range(-(m//2),m//2+1):
        if delta>=0: upper=sum((Fraction(j,m) for j in range(1,delta+1)),Fraction(0))
        else: upper=sum((Fraction(2*j,m) for j in range(-delta)),Fraction(0))
        assert upper<=Fraction(delta*delta+abs(delta),m)
        factorial_checks+=1

weight_rows=[]
for N,p,K in [(48,1,2),(96,2,2),(144,2,3)]:
    m=N//2;w=[1]*m+[0]*m
    assert p<=Fraction(N,12*K)
    # Asymmetric bad weights, repeated transfers and occupied/unoccupied swaps.
    edges=[(0,1),(0,m),(m,0),(m,m+1),(1,m+1),(m+1,1)]
    for ell in range(1,K+1):
        for pattern in product(range(len(edges)),repeat=ell):
            z=w.copy()
            for e in pattern:
                a,b=edges[e];z[a]+=p;z[b]-=p
            bad=[x for x in z if x not in (0,1)];b=len(bad);S=sum(bad)
            n0=z.count(0);n1=z.count(1);D=3*ell*p
            assert n1==m-S and n0==m-b+S
            assert b<=2*ell*p and abs(S)<=D and abs(S-b)<=D
            assert D<=m//2 and b<=m
            orbit=factorial(N)//prod(factorial(v) for v in Counter(z).values())
            assert prod(factorial(v) for v in Counter(bad).values())<=factorial(b)
            prefactor=Fraction(m**b,factorial(b));assert prefactor>=1
            exponent=Fraction(2*D*D+2*D,m)
            assert exponent<=Fraction(36*ell*ell*p*p,N)+1
            weight_rows.append(dict(N=N,p=p,ell=ell,bad=b,S=S,
                                   orbit_ratio=str(Fraction(orbit,comb(N,m))),exponent=str(exponent)))

# Conditional moments of sqrt(k) g^(j), whose off-diagonal conditional
# covariance is -delta_ij. Compare independently with univariate Hermites.
x=s.symbols('x0:3')
def conditional(labels):
    if not labels:return s.Integer(1)
    a,*tail=labels
    value=x[a]*conditional(tuple(tail))
    for j,b in enumerate(tail):
        if a==b:value-=conditional(tuple(tail[:j]+tail[j+1:]))
    return s.expand(value)
decoupling_checks=0
for k in range(1,7):
    for labels in product(range(2),repeat=k):
        want=prod(s.hermite_prob(labels.count(a),x[a]) for a in range(2))
        assert s.expand(conditional(labels)-want)==0;decoupling_checks+=1

def gaussian_expect(expr,variables):
    value=0
    for powers,c in s.Poly(s.expand(expr),*variables).terms():
        if any(n%2 for n in powers):continue
        value+=c*prod(s.factorial2(n-1) if n else 1 for n in powers)
    return value

# Finite rectangular tensors test every orientation in the iterated bound.
tensor_rows=[]
for k in (1,2,3):
    labels=list(product(range(2),repeat=k))
    coeff={i:s.Matrix(2,3,lambda a,b:(sum(i)+1)*(a-b)+(a+1)*(b+1)+sum(i)**2) for i in labels}
    g=s.symbols('g:'+str(2*k));z=s.zeros(2,3)
    coupled=s.zeros(2,3)
    for i,A in coeff.items():
        z+=A*prod(g[2*j+i[j]] for j in range(k))
        coupled+=A*prod(s.hermite_prob(i.count(a),x[a]) for a in range(2))
    flattenings=[]
    for mask in range(1<<k):
        I=[j for j in range(k) if mask>>j&1];J=[j for j in range(k) if not mask>>j&1]
        ir=list(product(range(2),repeat=len(I)));jc=list(product(range(2),repeat=len(J)))
        F=s.zeros(2*len(ir),3*len(jc))
        for ri,ii in enumerate(ir):
            for cj,jj in enumerate(jc):
                index=[0]*k
                for pos,a in zip(I,ii):index[pos]=a
                for pos,a in zip(J,jj):index[pos]=a
                F[2*ri:2*ri+2,3*cj:3*cj+3]=coeff[tuple(index)]
        flattenings.append(F)
    for r in range(1,3 if k==3 else 4):
        lhs=gaussian_expect(s.trace((z*z.T)**r),g)
        rhs=s.Rational((2*r-1)**(k*r),2**k)*sum(s.trace((F*F.T)**r) for F in flattenings)
        wick=gaussian_expect(s.trace((coupled*coupled.T)**r),x[:2])
        assert lhs<=rhs and wick<=k**(k*r)*lhs
        tensor_rows.append(dict(k=k,r=r,decoupled_moment=str(lhs),iterated_bound=str(rhs),
                                coupled_Wick_moment=str(wick),Jensen_bound=str(k**(k*r)*lhs)))

gamma,d=s.symbols('gamma d',real=True)
transfer={w:str(s.factor(moment(w,gamma,d))) for w in ['H-+H','--++','H-UU+H','--UU++']}
assert s.expand(moment('H-+H',gamma,d)-gamma**4*(1+gamma**2))==0
assert s.expand(moment('--++',gamma,d)-2*gamma**4)==0
result=dict(status='passed',scope='Exact finite diagnostics, not a numerical proof of the general asymptotic theorem.',
            CAR_rows=rows,seed_weights=seed_checks,partial_weights=partial_checks,
            factorial_rational_checks=factorial_checks,near_weight_rows=weight_rows,
            conditional_Wick_checks=decoupling_checks,rectangular_tensors=tensor_rows,
            limiting_transfer=transfer,parameters=dict(CAR_fixtures=5,max_conditional_degree=6,half_filling=True),
            environment=dict(python=platform.python_version(),sympy=s.__version__,arithmetic='exact integers, Fraction, SymPy'),
            input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in
              ['code/check_slater_polynomial_density.py','code/slater_covariant_transfer.py',
               'code/slater_joint_site_polynomial.py','results/susskind-slater-orbit-density/local-evidence.json']})
DIR.mkdir(exist_ok=True)
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',seed_weights=seed_checks,partial_weights=partial_checks,
                     factorial_checks=factorial_checks,near_weights=len(weight_rows),decoupling=decoupling_checks,
                     tensor_moments=len(tensor_rows),sha256=hashlib.sha256(OUT.read_bytes()).hexdigest())))
