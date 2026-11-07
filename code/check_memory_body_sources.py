"""Fresh exact finite-CAR and whole-field Poisson-source checks.

No historical main is executed; output is exclusive-create.
"""
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import hashlib
import json
import platform
import time

import numpy as np
from memory_body_source_chords import exponent, G_exponent, local, mul, sign_of_grouping
from singlet_coefficient import pairings

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT/'results/susskind-memory-oriented-transfer/body-source-checks.json'
PARAMETERS = dict(
    finite_systems=[[2, 1, 4, 'all'], [3, 1, 4, 'all'], [4, 2, 4, 36], [3, 1, 6, 24]],
    finite_blocks='two contiguous blocks with first length 1, floor(L/2), L-1',
    b0_layouts=[[2,2,2,0],[3,3,2,0],[4,2,2,0],[2,4,1,1],[2,2,2,2],
                [3,1,2,2],[4,4,0,0],[1,1,3,3],[2,2,0,0],[3,3,1,1]],
    arithmetic='Fraction local nilpotents; exact int64 CAR matrices with per-step bounds',
    finite_scope='Fixed p-set words for every listed pairing, not annealed N-limit approximation',
    source='p*d_F-D; d_F=L_F-2*(whole-field internal Gaussian pairs)',
    G_source='p*d_F-D-lambda_N*d_F*(d_F-1)',
    normalization='normalized full Fock trace; raw integer CAR words',
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def finite():
    rows = []
    maximum = 0
    for N, p, length, limit in PARAMETERS['finite_systems']:
        dimension = 1 << N
        identity = np.eye(dimension, dtype=np.int64)
        cs = []
        for i in range(N):
            c = np.zeros_like(identity)
            for state in range(dimension):
                if (state >> i) & 1:
                    c[state ^ (1 << i), state] = (-1)**((state & ((1 << i)-1)).bit_count())
            cs.append(c)
        Q = np.diag([s.bit_count() for s in range(dimension)])

        def matrix_product(a, b):
            bound = dimension*int(np.max(np.abs(a), initial=0))*int(np.max(np.abs(b), initial=0))
            assert bound < 2**55
            return a@b

        def D(x):
            return Q@x-sum(c.T@x@c for c in cs)

        subsets = list(combinations(range(N), p))
        creators = []
        for subset in subsets:
            operator = identity
            for i in subset:
                operator = matrix_product(operator, cs[i].T)
            creators.append(operator)
        if limit == 'all':
            assignments = list(product(range(len(subsets)), repeat=length))
        else:
            assignments = []
            for case in range(limit):
                assignments.append(tuple(((case+1)*(slot+1)+case//len(subsets)+slot*slot)
                                         %len(subsets) for slot in range(length)))
            assignments = sorted(set(assignments))
        layouts = [((0, j, 'D'), (j, length, 'D')) for j in sorted({1, length//2, length-1})]
        count = 0
        zero_components = negative_components = internal_cases = 0
        for pairing in pairings(tuple(range(length))):
            pi = dict(pairing) | {b:a for a,b in pairing}
            for assignment in assignments:
                sets = [subsets[index] for index in assignment]
                words = [matrix_product(creators[assignment[slot]], creators[assignment[pi[slot]]].T)
                         for slot in range(length)]
                sequence = [site for slot in range(length)
                            for site in (*sets[slot], *reversed(sets[pi[slot]]))]
                sign = sign_of_grouping(sequence)
                for blocks in layouts:
                    marked = []
                    degree = []
                    for a,b,_ in blocks:
                        x = identity
                        for word in words[a:b]:
                            x = matrix_product(x, word)
                        d = b-a-2*sum(a <= u < v < b for u,v in pairing)
                        degree.append(d)
                        marked.append([x, p*d*x-D(x)])
                    predicted = {0:F(sign)}
                    for site in range(N):
                        occupied = tuple(slot for slot, values in enumerate(sets) if site in values)
                        predicted = mul(predicted, local(length, pi, blocks, occupied))
                    for index,d in enumerate(degree):
                        predicted = mul(predicted, {0:F(1), 1 << index:F(p*d)})
                    for mask in range(4):
                        value = matrix_product(marked[0][mask & 1], marked[1][(mask >> 1) & 1])
                        maximum = max(maximum, int(np.max(np.abs(value), initial=0)))
                        actual = F(int(np.trace(value)), dimension)
                        assert actual == predicted.get(mask, F(0)), (N,p,pairing,assignment,blocks,mask,actual,predicted)
                        count += 1
                        zero_components += int(actual == 0)
                        negative_components += int(actual < 0)
                    internal_cases += int(any(d < b-a for d,(a,b,_) in zip(degree,blocks)))
        rows.append(dict(N=N,p=p,slots=length,pairings=len(list(pairings(tuple(range(length))))),
                         assignments=len(assignments),layouts=len(layouts),
                         scalar_coefficients=count,with_internal_pairs=internal_cases,
                         zero_components=zero_components,negative_components=negative_components,
                         exact=True))
    return rows, maximum


def no_bridge():
    rows = []
    for j,k,a,b in PARAMETERS['b0_layouts']:
        length=j+k+a+b
        blocks=((0,j,'D'),(j+a,j+a+k,'D'))
        accepted=[]
        for pairing in pairings(tuple(range(length))):
            if any(x<j and j+a<=y<j+a+k for x,y in pairing):
                continue
            i1=sum(x<y<j for x,y in pairing)
            i2=sum(j+a<=x<y<j+a+k for x,y in pairing)
            source=exponent(length,pairing,blocks)
            assert source.get(3,0)==4*i1*i2
            accepted.append(dict(pairing=pairing,internal=[i1,i2],
                                 source={str(mask):str(value) for mask,value in source.items()}))
        rows.append(dict(layout=[j,k,a,b],cases=len(accepted),
                         nonzero_connected=sum(row['internal'][0]*row['internal'][1]>0 for row in accepted),
                         diagrams=accepted))
    return rows


def main():
    if OUT.exists():
        raise FileExistsError(OUT)
    start=time.perf_counter()
    finite_rows, maximum=finite()
    b0=no_bridge()
    example=G_exponent(4,((0,1),(2,3)),((0,2,'D'),(2,4,'D')))
    assert example=={1:F(-2),2:F(-2),3:F(4)}
    inputs=['code/check_memory_body_sources.py','code/memory_body_source_chords.py',
            'code/memory_size_marked_chords.py','code/singlet_coefficient.py',
            'code/singlet_multimark.py','code/complex_syk_dictionary.py']
    record=dict(status='passed',schema_version=1,created_utc=datetime.now(timezone.utc).isoformat(),
                parameters=PARAMETERS,parameter_sha256=digest(PARAMETERS),finite=finite_rows,
                no_bridge=b0,maximum_absolute_integer=maximum,
                counts=dict(finite_scalar_coefficients=sum(r['scalar_coefficients'] for r in finite_rows),
                            no_bridge_diagrams=sum(r['cases'] for r in b0),
                            nonzero_no_bridge_connected=sum(r['nonzero_connected'] for r in b0)),
                static_G2=dict(exponent={str(i):str(x) for i,x in example.items()},
                               second_moment_over_C_squared='4*lambda+4*lambda^2',
                               centered_variance_over_C_squared='4*lambda'),
                environment=dict(python=platform.python_version(),numpy=np.__version__,platform=platform.platform()),
                input_sha256={name:sha(ROOT/name) for name in inputs},
                elapsed_seconds=time.perf_counter()-start,
                scope='Exact local CAR checks and fixed-diagram no-bridge identity checks. '
                      'The all-degree pair-Poisson and no-bridge proofs are separate. '
                      'No charged dynamic spectrum or exact pinching result is inferred from this file.')
    OUT.parent.mkdir(parents=True,exist_ok=True)
    with OUT.open('x',encoding='utf-8') as stream:
        json.dump(record,stream,indent=2)
        stream.write('\n')
    print(json.dumps(dict(status='passed',counts=record['counts'],maximum_absolute_integer=maximum,
                          output_sha256=sha(OUT),elapsed_seconds=record['elapsed_seconds'])))


if __name__=='__main__':
    main()
