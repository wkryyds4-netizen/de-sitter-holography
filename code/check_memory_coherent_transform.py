"""Independent finite identities for the adjacent-charge coherent transform."""
from pathlib import Path
from math import comb
from fractions import Fraction as F
import hashlib, json, platform, time
from memory_coherent_cross import (choose, charged, exchanged, rates,
                                  coherent_table, mixed_moments)
from cubic_linear_dynamics import weights

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT/'results/susskind-memory-coherent-cross'


def main():
    OUT.mkdir(exist_ok=True)
    target = OUT/'transform-checks.json'
    if target.exists():
        raise FileExistsError(target)
    start = time.time()
    embeddings = orthogonality = moments = recurrences = 0
    rows = []
    # General charge: compare pure-rank decomposition with the independently
    # known full embedded charged covariance eigenvalues (binomial formula).
    for N in range(3,17):
        for m in range(2,N):
            R = min(m-1,N-m)
            for p in range(R+1):
                C = comb(m,p)*comb(N-m+p,p)
                omega = weights(N,m,p)
                for t in range(R+1):
                    actual = sum(omega[a]*charged(N,m,a,t) for a in range(p+1))
                    expected = F(m*choose(m-t-1,p)*choose(N-m+p-t,p),(m-p)*C)
                    assert actual == expected, (N,m,p,t,actual,expected)
                    embeddings += 1
    for N in [4,6,8,10,12,14]:
        m = N//2
        ds,es,ell2,T,mu = coherent_table(N,m)
        dp,dm = comb(N,m),comb(N,m-1)
        x = lambda a: a*(N+1-a)
        y = lambda t: t*(N-t)
        K = (m-1)*(N+1)
        for t in range(m):
            f0 = sum((-1)**r*ds[r]*ell2[r]*T[r][t] for r in range(m))
            f1 = sum((-1)**r*ds[r]*ell2[r]*x(r)*T[r][t] for r in range(m))
            assert f0 == (-1)**(m-1+t)*m*dm*(m-t)
            assert f1 == (K-y(t))*f0
        for a in range(m):
            B,D = rates(N,a)
            assert 2*(B+D) == K-x(a)
            for b in range(m):
                gram = sum(es[t]*T[a][t]*T[b][t] for t in range(m))
                assert gram == (F(m*m*dp*dm,ds[a])/ell2[a] if a==b else 0)
                orthogonality += 1
                actual = B*(T[a][b]-(T[a+1][b] if a+1<m else 0))+D*(T[a][b]-(T[a-1][b] if a else 0))
                assert actual == y(b)*T[a][b]
                recurrences += 1
                predicted = mixed_moments(N,a,b)
                direct = [[sum(((-1 if (a+b-r)%2 else 1) if signed else 1)*
                               mu[a,b][r]*(x(a+b)-x(r))**k for r in range(m))
                           for k in range(2)] for signed in [False,True]]
                assert direct == list(predicted), (N,a,b,direct,predicted)
                assert sum(mu[a,b]) == 1
                assert all(v==0 for r,v in enumerate(mu[a,b]) if r>a+b)
                xi = exchanged(N,a,b)
                contiguous = T[a][b]-(F(b,N+1-2*b)*(T[a][b-1]-T[a][b]) if b else 0)
                assert xi == contiguous
                moments += 4
        rows.append(dict(N=N,m=m,pairs=m*m,minimum_coherent_weight=str(min(v for values in mu.values() for v in values))))
        print('exact transform',N,'passed',flush=True)
    names = ['code/memory_coherent_cross.py','code/check_memory_coherent_transform.py',
             'code/cubic_linear_dynamics.py']
    data = dict(status='passed',rows=rows,embedded_covariance_identities=embeddings,
                orthogonality_identities=orthogonality,signed_unsigned_moments=moments,
                difference_equations=recurrences,elapsed_seconds=time.time()-start,
                scope='Exact finite rational identities; coherent weights are not assumed positive. These tests are finite checks of the derivation, not a proof by extrapolation.',
                environment=dict(python=platform.python_version()),
                input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print({k:data[k] for k in ['embedded_covariance_identities','signed_unsigned_moments']})


if __name__ == '__main__':
    main()
