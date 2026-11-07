"""Exact finite CAR weights, orbit covariance and projective density moments.

No Haar sampling, old generator execution or asymptotic fitting.
"""
from pathlib import Path
from itertools import combinations, permutations, product
from functools import lru_cache
from math import comb, prod
from datetime import datetime,timezone
import hashlib,json,platform
import sympy as s

ROOT=Path(__file__).resolve().parent.parent
DIR=ROOT/'results/susskind-slater-orbit-density';OUT=DIR/'exact-checks.json'
if OUT.exists():raise FileExistsError(OUT)

def hop(state,A,B):
    sign=1
    for i in B:
        if not (state>>i)&1:return None
        sign*=(-1)**((state&((1<<i)-1)).bit_count());state^=1<<i
    for i in reversed(A):
        if (state>>i)&1:return None
        sign*=(-1)**((state&((1<<i)-1)).bit_count());state|=1<<i
    return state,sign

def wedge(U,k):
    inds=list(combinations(range(U.rows),k))
    return s.Matrix([[U.extract(a,b).det() for b in inds] for a in inds])

seed_rows=[];weight_checks=0;frame_counts=0;stabilizer_entries=0
for N,m,p in ((4,2,1),(4,2,2),(6,3,1),(6,3,2),(8,4,2)):
    subsets=list(combinations(range(N),p));states=list(combinations(range(N),m))
    bits=[sum(1<<i for i in x) for x in states];index={x:i for i,x in enumerate(bits)}
    O=(1<<m)-1;L=comb(m,p)*comb(N-m,p);d=len(states)
    seed={}
    for ai,A in enumerate(subsets):
        for bi,B in enumerate(subsets):
            if not all(i>=m for i in A) or not all(i<m for i in B):continue
            out,sign=hop(O,A,B);t=index[out]
            weight=tuple(int(i in B)-int(i in A)+int(out>>i&1) for i in range(N))
            assert weight==tuple(int(i<m) for i in range(N))
            seed[ai,bi,t]=s.Integer(sign);weight_checks+=1
    assert len(seed)==L
    # Exact permutation twirl of the diagonal output covariance Q_p(O)/L.
    for target in bits:
        count=sum((target&~source).bit_count()==p for source in bits)
        assert count==L;frame_counts+=1
    record=dict(N=N,m=m,p=p,sector_dimension=d,selected_labels=L,weight_orbit_size=comb(N,m))
    if N==4:
        rotations=[]
        U=s.eye(N)
        for a,b in ((0,1),(m,m+1)):
            U[a,a]=U[b,b]=s.Rational(3,5);U[a,b]=-s.Rational(4,5);U[b,a]=s.Rational(4,5)
        rotations.append(U)
        rotations.append(s.diag(*[s.I**i for i in range(N)]))
        for U in rotations:
            Up=wedge(U,p);Um=wedge(U,m);character=U[:m,:m].det()
            tensors=[s.zeros(len(subsets)) for _ in states]
            for (a,b,t),v in seed.items():tensors[t][a,b]=v
            transformed=[Up.conjugate()*T*Up.T for T in tensors]
            for st in range(d):
                T=sum((Um[st,t]*transformed[t] for t in range(d)),s.zeros(len(subsets)))
                assert T==character*tensors[st]
                stabilizer_entries+=len(subsets)**2
        record['exact_stabilizer_matrices']=2
    seed_rows.append(record)

@lru_cache(None)
def projector_moment(N,pairs):
    k=len(pairs)
    numerator=sum(all(pairs[i][0]==pairs[pi[i]][1] for i in range(k))
                  for pi in permutations(range(k)))
    return s.Rational(numerator,prod(N+j for j in range(k)))

def closed_density(J):
    N=J.rows;tr=s.trace(J)
    return (N*J*J-2*tr*J+(tr**2+s.trace(J*J))/(N+1)*s.eye(N))/((N-1)*(N+2))

haar_rows=[];haar_entries=0
for N in (2,3,4):
    for fixture in (0,1):
        A=s.Matrix(N,N,lambda i,j:(i+1)*(j+2)+(i-j)**2+fixture*(i+2*j))
        J=A+A.T+s.I*(A-A.T)
        actual=s.zeros(N)
        for i,j in product(range(N),repeat=2):
            value=(J*J)[i,j]/N
            for a,b,c in product(range(N),repeat=3):
                value-=J[a,b]*J[c,j]*projector_moment(N,((i,a),(b,c)))
                value-=J[i,a]*J[b,c]*projector_moment(N,((a,b),(c,j)))
            for a,b,c,d in product(range(N),repeat=4):
                value+=J[a,b]*J[c,d]*projector_moment(N,((i,a),(b,c),(d,j)))
            actual[i,j]=s.simplify(N*value/(N-1))
        difference=(actual-closed_density(J)).applyfunc(s.simplify)
        assert difference==s.zeros(N),(N,fixture,difference)
        assert all(s.simplify(actual.extract(ids,ids).det())>=0
                   for k in range(1,N+1) for ids in combinations(range(N),k))
        haar_entries+=N*N
        haar_rows.append(dict(N=N,m=1,p=1,fixture=fixture,
                              trace_density=str(s.trace(actual)/N),exact_identity=True))

# Independent Gaussian coordinate Wick contraction of the closed matrix density.
moment_rows=[];gaussian_pairs=0
for N in (2,3,4,5):
    bs=[];variances=[]
    for i in range(N):
        B=s.zeros(N);B[i,i]=1;bs.append(B);variances.append(s.Integer(1))
    for i in range(N):
        for j in range(i):
            B=s.zeros(N);B[i,j]=B[j,i]=1;bs.append(B);variances.append(s.Rational(1,2))
            B=s.zeros(N);B[i,j]=s.I;B[j,i]=-s.I;bs.append(B);variances.append(s.Rational(1,2))
    ks=[B-s.trace(B)/N*s.eye(N) for B in bs]
    coeff=[];mean=s.zeros(N)
    for a,A in enumerate(ks):
        line=[]
        for b,B in enumerate(ks):
            C=s.Rational(N,(N-1)*(N+2))*((A*B+B*A)/2+s.trace(A*B)/(N*(N+1))*s.eye(N))
            line.append(C)
            if a==b:mean+=variances[a]*C
        coeff.append(line)
    assert mean==s.eye(N)
    second=s.trace(mean*mean)/N
    for a,b in product(range(len(bs)),repeat=2):
        second+=2*variances[a]*variances[b]*s.trace(coeff[a][b]*coeff[a][b])/N
        gaussian_pairs+=1
    expected=(2*N**4-5*N*N+3+s.Rational(2*N+3,(N+1)**2)*(N**4-1))/((N-1)*(N+2))**2
    assert s.simplify(second-expected)==0 and 1<=second<=3
    moment_rows.append(dict(N=N,real_Gaussian_coordinates=len(bs),exact_second_moment=str(second),bound=3))
chi_rows=[]
for r in range(1,9):
    exact=s.Rational(prod(3+2*j for j in range(r)),3**r)
    bound=prod(2*j+1 for j in range(r))
    assert exact<=bound
    chi_rows.append(dict(r=r,N2_exact_density_moment=str(exact),universal_bound=bound))

DIR.mkdir(parents=True,exist_ok=True)
files=[Path(__file__),ROOT/'results/susskind-slater-neutral-slope/local-evidence.json']
result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),
    seed_systems=seed_rows,projective_integrals=haar_rows,Gaussian_Wick_moments=moment_rows,N2_high_moments=chi_rows,
    counts=dict(seed_weights=weight_checks,permutation_frame_counts=frame_counts,stabilizer_entries=stabilizer_entries,
                exact_Haar_entries=haar_entries,Gaussian_coordinate_pairs=gaussian_pairs),
    parameters=dict(seed_fixtures=[[4,2,1],[4,2,2],[6,3,1],[6,3,2],[8,4,2]],
                    projective_N=[2,3,4],Gaussian_Wick_N=[2,3,4,5],N2_moment_orders=list(range(1,9)),
                    sampling='none; exact Haar monomials through third order'),
    environment=dict(python=platform.python_version(),sympy=s.__version__,arithmetic='SymPy Rational and exact complex matrices'),
    input_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
    scope='Finite diagnostics for seed torus weight/stabilizer and density normalization. General orbit-frame and moment bounds require the analytic representation-theoretic proof. The m=1 fixtures beyond N=2 are general finite-N checks, not half-filled DS examples.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',counts=result['counts'],sha256=hashlib.sha256(OUT.read_bytes()).hexdigest())))
