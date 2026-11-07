"""Independent finite charge polynomials and joint-source tilt checks.

All listed finite support assignments are exhaustive. Local 2x2 nilpotent
products are compared with full-Fock CAR matrices charge by charge. No
historical scientific module or main is imported or run.
"""
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
from pathlib import Path
import hashlib
import json
import platform
import time

import mpmath as mp
import numpy as np


ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT/'results/susskind-memory-joint-regression/charge-tilt-checks.json'
PARENTS = ['results/susskind-memory-oriented-transfer/local-evidence.json',
           'results/susskind-memory-oriented-transfer/charge-background-checks.json']
LAYOUTS = [
    ('GG', [(0, 2, 'D'), (2, 4, 'D')], []),
    ('GL', [(0, 2, 'D'), (2, 4, 'L')], [(2, 4)]),
    ('LG', [(0, 2, 'L'), (2, 4, 'D')], [(0, 2)]),
    ('G_Lsub', [(0, 2, 'D'), (2, 3, 'L')], [(2, 4)]),
    ('Lsub', [(0, 1, 'L')], [(0, 2)]),
    ('G_only', [(0, 2, 'D')], []),
]
PARAMETERS = dict(
    finite_systems=[[2, 1], [3, 1], [4, 2]], slots=4,
    support_assignment_scope='All binom(N,p)^4 ordered p-set assignments in each listed system',
    layouts=LAYOUTS, finite_charge_polynomial='Every coefficient of u^Q and every source subset',
    G_source='p*d-D-lambda_N*d*(d-1), lambda_N=p^2/N',
    neutral_source='p*mark_length-L; marked subblock remains inside the recorded Wick forbidden block',
    exact_occupations=['1/3', '1/2', '3/5', '2/3'], exact_singleton_p=[1, 2, 7, 32],
    recorded_charge_moments=list(range(5)),
    scalar_p=[8, 32, 128, 512], scalar_lambda=['1/4', '1/2', '1'], scalar_chords=[1, 2, 4],
    scalar_t=[['-2','0'], ['-1','0'], ['0','0'], ['1','0'], ['2','0'], ['1/2','3/4']],
    G2_systems=[[256,8], [4096,32], [65536,128], [1024,32], [4160,32]],
    G2_t=[['-1/2','0'], ['0','0'], ['1/2','0'], ['0','1/2']],
    precisions=[70,110],
    finite_normalization='Raw standard-CAR support words; normalized full Fock trace is 2^-N',
    G2_normalization='Divide zero-chaos G2 by exact common-variance C_(N/2); every nonzero charge sector summed',
)
ZERO = (F(0),)*4
IDENTITY = (F(1), F(0), F(0), F(1))
CREATOR = (F(0), F(0), F(1), F(0))
ANNIHILATOR = (F(0), F(1), F(0), F(0))
MAX_INTEGER = 0


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix_mul(a, b):
    global MAX_INTEGER
    ma, mb = int(np.max(np.abs(a), initial=0)), int(np.max(np.abs(b), initial=0))
    assert a.shape[1]*ma*mb < 2**60
    value = a@b
    MAX_INTEGER = max(MAX_INTEGER, int(np.max(np.abs(value), initial=0)))
    return value


def mm(a, b):
    return (a[0]*b[0]+a[1]*b[2], a[0]*b[1]+a[1]*b[3],
            a[2]*b[0]+a[3]*b[2], a[2]*b[1]+a[3]*b[3])


def loss(a, kind):
    if kind == 'D':
        return (F(0), a[1]/2, a[2]/2, a[3]-a[0])
    return ((a[0]-a[3])/2, a[1]/2, a[2]/2, (a[3]-a[0])/2)


def matrix_poly_mul(a, b):
    out = {}
    for i, x in a.items():
        for j, y in b.items():
            if not i & j:
                value = mm(x, y)
                out[i | j] = tuple(u+v for u, v in zip(out.get(i | j, ZERO), value))
    return out


def scalar_mul(a, b):
    out = Counter()
    for i, x in a.items():
        for j, y in b.items():
            if not i & j:
                out[i | j] += x*y
    return {mask: value for mask, value in out.items() if value}


def scalar_add(a, b):
    out = Counter(a)
    out.update(b)
    return {mask: value for mask, value in out.items() if value}


def scalar_scale(a, value):
    return {mask: F(value)*x for mask, x in a.items() if value*x}


def inverse(a):
    constant = a[0]
    rest = scalar_scale({mask: value for mask, value in a.items() if mask}, -1/constant)
    return scalar_scale(scalar_add(scalar_add({0:F(1)}, rest), scalar_mul(rest, rest)), 1/constant)


def serial(a):
    return {str(mask):str(value) for mask,value in sorted(a.items())}


def pairings(slots):
    if not slots:
        yield ()
    else:
        for i in range(1, len(slots)):
            for rest in pairings(slots[1:i]+slots[i+1:]):
                yield ((slots[0], slots[i]),)+rest


def make_car(N, p):
    dimension = 1 << N
    cs = []
    for index in range(N):
        c = np.zeros((dimension, dimension), dtype=np.int64)
        for state in range(dimension):
            if state & (1 << index):
                c[state ^ (1 << index), state] = (-1)**((state & ((1 << index)-1)).bit_count())
        cs.append(c)
    eye = np.eye(dimension, dtype=np.int64)
    charges = np.array([state.bit_count() for state in range(dimension)], dtype=np.int64)
    Q = np.diag(charges)
    for i in range(N):
        for j in range(N):
            assert np.array_equal(cs[i]@cs[j].T+cs[j].T@cs[i], eye if i == j else 0*eye)
    subsets = list(combinations(range(N), p))
    create = []
    for subset in subsets:
        value = eye
        for index in subset:
            value = matrix_mul(value, cs[index].T)
        create.append(value)
    return cs, eye, charges, Q, subsets, create


def local_charge_polynomial(N, p, pairing, blocks, sets):
    pi = dict(pairing) | {b:a for a,b in pairing}
    starts = {a:(index,b,kind) for index,(a,b,kind) in enumerate(blocks)}
    sequence = [site for slot in range(4) for site in (*sets[slot], *reversed(sets[pi[slot]]))]
    sign = (-1)**sum(sequence[i] > sequence[j] for i in range(len(sequence)) for j in range(i+1,len(sequence)))
    total = {(0,0):F(sign)}
    for site in range(N):
        local, position = {0:IDENTITY}, 0
        while position < 4:
            index, end, kind = starts.get(position, (None,position+1,None))
            value = IDENTITY
            for slot in range(position, end):
                if site in sets[slot]:
                    value = mm(value, CREATOR)
                if site in sets[pi[slot]]:
                    value = mm(value, ANNIHILATOR)
            factor = {0:value}
            if index is not None:
                factor[1 << index] = tuple(-x for x in loss(value,kind))
            local = matrix_poly_mul(local, factor)
            position = end
        factor = {(mask,charge):value for mask,matrix in local.items()
                  for charge,value in [(0,matrix[0]),(1,matrix[3])] if value}
        new = Counter()
        for (mask,charge),value in total.items():
            for (other,occupation),coefficient in factor.items():
                if not mask & other:
                    new[mask | other,charge+occupation] += value*coefficient
        total = {key:value for key,value in new.items() if value}
    for index,(a,b,kind) in enumerate(blocks):
        degree = b-a-2*sum(a <= x < y < b for x,y in pairing)
        scalar = F(p*degree)-F(p*p,N)*degree*(degree-1) if kind == 'D' else F(p*(b-a))
        new = Counter(total)
        for (mask,charge),value in total.items():
            if not mask & (1 << index):
                new[mask | (1 << index),charge] += scalar*value
        total = {key:value for key,value in new.items() if value}
    return total


def finite_checks():
    summaries, examples = [], []
    for N,p in PARAMETERS['finite_systems']:
        cs, eye, charges, Q, subsets, creators = make_car(N,p)
        den = 2*N
        counts = Counter()
        for assignment in product(range(len(subsets)), repeat=4):
            sets = [subsets[index] for index in assignment]
            for pairing in pairings(tuple(range(4))):
                pi = dict(pairing) | {b:a for a,b in pairing}
                operators = [matrix_mul(creators[assignment[slot]], creators[assignment[pi[slot]]].T)
                             for slot in range(4)]
                for name,blocks,wicks in LAYOUTS:
                    if any(a <= x < y < b for a,b in wicks for x,y in pairing):
                        continue
                    starts = {a:(index,b) for index,(a,b,kind) in enumerate(blocks)}
                    ordinary, marked = [], []
                    for a,b,kind in blocks:
                        value = eye
                        for operator in operators[a:b]:
                            value = matrix_mul(value,operator)
                        assert np.array_equal(Q@value,value@Q)
                        degree = b-a-2*sum(a <= x < y < b for x,y in pairing)
                        if kind == 'D':
                            D = Q@value-sum(c.T@value@c for c in cs)
                            mark = den*(p*degree*value-D)-2*p*p*degree*(degree-1)*value
                        else:
                            assert degree == b-a
                            twiceL = N*value-sum(c.T@value@c+c@value@c.T for c in cs)
                            mark = N*(2*p*(b-a)*value-twiceL)
                        ordinary.append(value)
                        marked.append(mark)
                    predicted = local_charge_polynomial(N,p,pairing,blocks,sets)
                    observed = {}
                    for mask in range(1 << len(blocks)):
                        value, position = eye, 0
                        while position < 4:
                            if position in starts:
                                index,end = starts[position]
                                value = matrix_mul(value,marked[index] if mask >> index & 1 else ordinary[index])
                                position = end
                            else:
                                value = matrix_mul(value,operators[position])
                                position += 1
                        diagonal = np.diag(value)
                        for charge in range(N+1):
                            actual = F(int(diagonal[charges == charge].sum()), den**mask.bit_count())
                            assert actual == predicted.get((mask,charge),F()), (N,p,assignment,pairing,name,mask,charge)
                            observed[mask,charge] = actual
                            counts['charge_source_coefficients'] += 1
                            counts['nonzero_coefficients'] += actual != 0
                            counts['negative_coefficients'] += actual < 0
                    counts['support_pairing_layout_cases'] += 1
                    counts['cases_'+name] += 1
                    if len([row for row in examples if row['N'] == N]) < 6 and any(
                            value for (mask,charge),value in observed.items() if mask):
                        moments = []
                        for mask in range(1 << len(blocks)):
                            values = [sum(observed[mask,m]*F(2*m-N,2*p)**r for m in range(N+1))/2**N
                                      for r in PARAMETERS['recorded_charge_moments']]
                            moments.append(dict(mask=mask, X_moments=list(map(str,values))))
                        examples.append(dict(N=N,p=p,assignment=assignment,pairing=pairing,layout=name,
                            charge_coefficients={str(mask):[str(observed[mask,m]) for m in range(N+1)]
                                                 for mask in range(1 << len(blocks))}, moments=moments))
        expected_assignments = comb(N,p)**4
        assert counts['support_pairing_layout_cases'] == 14*expected_assignments
        assert counts['charge_source_coefficients'] == 46*(N+1)*expected_assignments
        summaries.append(dict(N=N,p=p,assignments=expected_assignments,counts=dict(counts),
                              all_charge_and_source_coefficients_equal=True))
        print(json.dumps(dict(finite=[N,p],counts=dict(counts),status='passed')),flush=True)
    return dict(systems=summaries,examples=examples,
                counts={key:sum(row['counts'][key] for row in summaries)
                        for key in ['support_pairing_layout_cases','charge_source_coefficients',
                                    'nonzero_coefficients','negative_coefficients']})


def exact_singletons():
    rows, collisions = [], []
    for occupation in PARAMETERS['exact_occupations']:
        rho = F(occupation)
        ratio = rho/(1-rho)
        for p in PARAMETERS['exact_singleton_p']:
            factor = scalar_mul({0:F(1),1:F(-1)},{0:F(1),1:ratio})
            value = {0:F(1)}
            for _ in range(p):
                value = scalar_mul(value,factor)
            expected = {0:F(1)}
            if p*(ratio-1):
                expected[1] = p*(ratio-1)
            assert value == expected
            # Two independent nilpotents retain their product; no square of a
            # single source is allowed.
            second = {0:F(1),2:p*(ratio-1)}
            combined = scalar_mul(value,second)
            assert combined.get(3,F()) == p*p*(ratio-1)**2
            assert (4*rho*(1-rho))**p == (4*ratio/(1+ratio)**2)**p
            rows.append(dict(rho=occupation,p=p,exp_t_over_p=str(ratio),one=serial(value),two=serial(combined)))
        def dressed(kind,bit):
            return ({},{0:F(1),bit:F(-1)}) if kind == 'n' else ({0:F(1)},{bit:F(1)})
        def tr(matrix):
            return scalar_add(scalar_scale(matrix[0],1-rho),scalar_scale(matrix[1],rho))
        all_ratios = {}
        for first in ['n','bar']:
            for second in ['n','bar']:
                a,b = dressed(first,1),dressed(second,2)
                numerator = tr(tuple(scalar_mul(x,y) for x,y in zip(a,b)))
                denominator = scalar_mul(tr(a),tr(b))
                value = scalar_mul(numerator,inverse(denominator))
                all_ratios[first+','+second] = value
        expected = {'n,n':{0:1/rho},'n,bar':{2:1/(1-rho)},'bar,n':{1:1/(1-rho)},
                    'bar,bar':{0:1/(1-rho),1:-rho/(1-rho)**2,2:-rho/(1-rho)**2,
                               3:rho/(1-rho)**3}}
        assert all_ratios == expected
        collisions.append(dict(rho=occupation,ratios={key:serial(value) for key,value in all_ratios.items()}))
    return dict(singleton_rows=rows,collision_rows=collisions,
                singleton_cases=len(rows),collision_ratios=4*len(collisions))


def mfrac(value):
    value = F(value)
    return mp.mpf(value.numerator)/value.denominator


def complex_record(value):
    return dict(real=mp.nstr(mp.re(value),60),imag=mp.nstr(mp.im(value),60))


def scalar_tilt(p,lam,ell,t,precision):
    with mp.workdps(precision):
        t = mfrac(t[0])+1j*mfrac(t[1])
        lam = mfrac(lam)
        N = p*p/lam
        z = t/(2*p)
        singleton_log = -2*p*ell*mp.log(mp.cosh(z))
        gaussian_log = t*t/(8*lam)
        charge_log = N*mp.log(mp.cosh(z))
        source = p*mp.expm1(t/p)
        source_bound = abs(t)**2*mp.exp(abs(t)/p)/(2*p)
        singleton_bound = p*ell*abs(z)**2/(1-abs(z)**2)
        charge_bound = N*abs(z)**4/(4*(1-abs(z)**2))
        tolerance = mp.mpf(10)**(-precision+8)
        assert abs(source-t) <= source_bound+tolerance
        assert abs(singleton_log) <= singleton_bound+tolerance
        assert abs(charge_log-gaussian_log) <= charge_bound+tolerance
        return dict(singleton_log=singleton_log,charge_log_error=charge_log-gaussian_log,
                    source_coefficient=source,source_error=source-t,
                    combined_MGF=mp.exp(charge_log+singleton_log),target_MGF=mp.exp(gaussian_log),
                    source_bound=source_bound,singleton_log_bound=singleton_bound,charge_log_bound=charge_bound)


def scalar_checks():
    rows = []
    with mp.workdps(110):
        for p,lam,ell,t in product(PARAMETERS['scalar_p'],PARAMETERS['scalar_lambda'],
                                   PARAMETERS['scalar_chords'],PARAMETERS['scalar_t']):
            low = scalar_tilt(p,lam,ell,t,70)
            high = scalar_tilt(p,lam,ell,t,110)
            error = max(abs(high[key]-low[key]) for key in high)
            assert error < mp.mpf('1e-58')
            rows.append(dict(p=p,lambda_N=lam,chords=ell,t=t,
                             values={key:complex_record(value) for key,value in high.items()},
                             maximum_precision_difference=mp.nstr(error,12),all_three_bounds=True))
    return dict(rows=rows,cases=len(rows),precisions=[70,110])


def G2_tilt(N,p,t,precision):
    with mp.workdps(precision):
        t = mfrac(t[0])+1j*mfrac(t[1])
        lam = mp.mpf(p*p)/N
        current = mp.mpf(1)
        for j in range(p):
            current *= mp.mpf((j+1)*(N-j))/((N//2-j)*(N//2+p-j))
        probability = mp.mpf(comb(N,p))/mp.power(2,N)*mp.exp(t*mp.mpf(2*p-N)/(2*p))
        ratio = mp.exp(t/p)
        first = second = mp.mpc(0)
        for m in range(p,N+1):
            G0 = -current*mp.mpf(p*(N-2*m+p+1))/(N-m+1)
            first += probability*G0
            second += probability*G0*G0
            if m<N:
                current *= mp.mpf((m+1)*(N-m))/((m+1-p)*(N-m+p))
                probability *= mp.mpf(N-m)/(m+1)*ratio
        M = mp.exp(t*t/(8*lam))
        return dict(mean=first,second=second,finite_charge_MGF=mp.cosh(t/(2*p))**N,
                    target_charge_MGF=M,target_mean=(-2*lam+t)*M,
                    target_second=((-2*lam+t)**2+4*lam)*M)


def G2_checks():
    rows = []
    with mp.workdps(110):
        for (N,p),t in product(PARAMETERS['G2_systems'],PARAMETERS['G2_t']):
            low = G2_tilt(N,p,t,70)
            high = G2_tilt(N,p,t,110)
            error = max(abs(high[key]-low[key]) for key in high)
            assert error < mp.mpf('1e-57')
            rows.append(dict(N=N,p=p,lambda_N=str(F(p*p,N)),t=t,
                all_charge_sectors=N+1,exactly_zero_low_sectors=p,
                values={key:complex_record(value) for key,value in high.items()},
                mean_absolute_error=mp.nstr(abs(high['mean']-high['target_mean']),40),
                second_absolute_error=mp.nstr(abs(high['second']-high['target_second']),40),
                maximum_precision_difference=mp.nstr(error,12)))
            print(json.dumps(dict(G2=[N,p],t=t,status='passed',second_error=rows[-1]['second_absolute_error'])),flush=True)
    return dict(rows=rows,cases=len(rows),precisions=[70,110],
                scope='Complete finite-binomial tilted zero-chaos checks, including complex t. '
                      'Finite errors remain; precision agreement is not an interval certificate or limit proof.')


def main():
    if OUT.exists():
        raise FileExistsError(OUT)
    started = time.perf_counter()
    parents = {name:sha(ROOT/name) for name in PARENTS}
    finite = finite_checks()
    singleton = exact_singletons()
    scalar = scalar_checks()
    background = G2_checks()
    record = dict(schema_version=1,status='passed',created_utc=datetime.now(timezone.utc).isoformat(),
        parameters=PARAMETERS,parameter_sha256=digest(PARAMETERS),
        finite=finite,singleton=singleton,scalar=scalar,G2_background=background,
        counts=dict(**finite['counts'],exact_singletons=singleton['singleton_cases'],
                    exact_collision_ratios=singleton['collision_ratios'],
                    scalar_complex_tilt_cases=scalar['cases'],G2_full_charge_tilt_cases=background['cases']),
        maximum_absolute_integer=MAX_INTEGER,integer_operation_guard=2**60,
        parents_sha256=parents,
        input_sha256={'code/check_memory_charge_tilt.py':sha(Path(__file__)),**parents},
        environment=dict(python=platform.python_version(),numpy=np.__version__,mpmath=mp.__version__,platform=platform.platform()),
        elapsed_seconds=time.perf_counter()-started,
        scope='Exact fixed-word full-charge CAR/source polynomials, rational occupations, and '
              'two-precision finite tilt factors. One/two sources and a neutral smaller marked '
              'subblock are retained. The separate complex-local-uniform and exponential-tail '
              'proof establishes joint polynomial moments/UI; no general joint energy regression, '
              'charged spectral kernel, growing-degree or infinite-time interchange is inferred here.')
    OUT.parent.mkdir(parents=True,exist_ok=True)
    with OUT.open('x',encoding='utf-8') as stream:
        stream.write(json.dumps(record,indent=2)+'\n')
    print(json.dumps(dict(status='passed',counts=record['counts'],maximum_absolute_integer=MAX_INTEGER,
                         elapsed_seconds=record['elapsed_seconds'],output_sha256=sha(OUT))),flush=True)


if __name__ == '__main__':
    main()
