"""Charge-static zero-chaos background of the whole-field oriented G2.

Exact binomial/overlap identities, full-charge high-precision sums and an
independent two-site-source nilpotent calculation. No cutoff or disorder
sampling is used; asymptotic claims require the separate analytic proof.
"""
from datetime import datetime, timezone
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path
import hashlib
import json
import platform
import time

import mpmath as mp


ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT/'results/susskind-memory-oriented-transfer/charge-background-checks.json'
PARAMETERS = dict(
    exact_even_N=list(range(2, 33, 2)), exact_p='1..min(6,N/2)',
    high_precision_p=[8, 32, 128], fixed_lambda=['1/4', '1/2', '1'],
    drifting_N='4p^2+2p at p=8,32,128', precisions=[70, 110],
    central_x=['-1', '-1/2', '0', '1/2', '1', '2'],
    Js=1, charge_weights='2^-N binom(N,m), every m=0..N',
    normalization='Exact calculations use C_(N/2)=1; numerical physical C0 uses the stated common w_N',
    physical_w='4^p (Js^2/lambda_N) (p!/N^p)^2, lambda_N=p^2/N',
    background='G2_zero_chaos=-D C_Q; (D C_Q)_m=m(C_m-C_(m-1))',
)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def c_ratio(N, p, m):
    return F(comb(m, p)*comb(N-m+p, p), comb(N//2, p)*comb(N//2+p, p))


def exact_checks():
    rows, sectors = [], 0
    for N in PARAMETERS['exact_even_N']:
        for p in range(1, min(6, N//2)+1):
            denominator = comb(N//2, p)*comb(N//2+p, p)
            mean_C = mean_D = second_D = shifted_D = F()
            maximum_ratio = F()
            for m in range(N+1):
                weight = F(comb(N, m), 2**N)
                current = c_ratio(N, p, m)
                previous = c_ratio(N, p, m-1) if m else F()
                difference = m*(current-previous)
                formula = current*F(p*(N-2*m+p+1), N-m+1)
                assert difference == formula
                maximum_ratio = max(maximum_ratio, current)
                mean_C += weight*current
                mean_D += weight*difference
                second_D += weight*difference*difference
                shifted_D += weight*(2*m-N)*current
                sectors += 1
            # Independent support overlap I intersect J; each term fixes p
            # occupied and p-r empty sites, leaving all others free.
            overlap_C = overlap_D = F()
            for r in range(max(0, 2*p-N), p+1):
                weight = F(comb(N, p)*comb(p, r)*comb(N-p, p-r), denominator*2**(2*p-r))
                overlap_C += weight
                overlap_D += r*weight
            assert mean_C == overlap_C and mean_D == overlap_D == shifted_D
            variance = second_D-mean_D**2
            assert variance > 0
            rows.append(dict(N=N, p=p, lambda_N=str(F(p*p, N)),
                             mean_C_over_C0=str(mean_C), mean_G0_over_C0=str(-mean_D),
                             second_G0_over_C0_squared=str(second_D),
                             variance_G0_over_C0_squared=str(variance),
                             maximum_C_over_C0=str(maximum_ratio),
                             exact_difference_identity=True, exact_trace_shift_identity=True,
                             independent_support_overlap_identity=True))
    return dict(rows=rows, sector_identities=sectors, systems=len(rows))


def p_add(a, b):
    out = dict(a)
    for mask, value in b.items():
        out[mask] = out.get(mask, F())+value
    return {mask: value for mask, value in out.items() if value}


def p_scale(a, value):
    return {mask: F(value)*coefficient for mask, coefficient in a.items() if value*coefficient}


def p_mul(a, b):
    out = {}
    for left, x in a.items():
        for right, y in b.items():
            if not left & right:
                key = left | right
                out[key] = out.get(key, F())+x*y
    return {mask: value for mask, value in out.items() if value}


def p_inverse(a):
    constant = a[0]
    rest = p_scale({key: value for key, value in a.items() if key}, -1/constant)
    return p_scale(p_add(p_add({0: F(1)}, rest), p_mul(rest, rest)), 1/constant)


def serialize(poly):
    return {str(mask): str(value) for mask, value in sorted(poly.items())}


def pair_source_checks():
    # Diagonal 2x2 matrices [vacuum, occupied]. Dn=n; D(1-n)=-n.
    def dressed(kind, bit):
        return ({}, {0: F(1), bit: F(-1)}) if kind == 'n' else ({0: F(1)}, {bit: F(1)})

    def trace(matrix):
        return p_scale(p_add(*matrix), F(1, 2))

    def ratio(first, second):
        numerator = trace(tuple(p_mul(a, b) for a, b in zip(first, second)))
        denominator = p_mul(trace(first), trace(second))
        return p_mul(numerator, p_inverse(denominator))

    cross = {}
    total_cross = {}
    for first in ['n', 'bar']:
        for second in ['n', 'bar']:
            value = ratio(dressed(first, 1), dressed(second, 2))
            cross[first+','+second] = serialize(value)
            total_cross = p_add(total_cross, p_add(value, {0: F(-1)}))
    assert total_cross == {3: F(4)}
    own = {}
    kappa = total_cross
    for bit in [1, 2]:
        # Coincident I/J within one internal chord gives local word n,
        # not the product n(1-n). Divide by the two separated-site traces.
        denominator = p_mul(trace(dressed('n', bit)), trace(dressed('bar', bit)))
        value = p_mul(trace(dressed('n', bit)), p_inverse(denominator))
        assert value == {0: F(2), bit: F(-2)}
        own[str(bit)] = serialize(value)
        kappa = p_add(kappa, p_add(value, {0: F(-2)}))
    assert kappa == {1: F(-2), 2: F(-2), 3: F(4)}
    rows = []
    for lam in [F(1, 4), F(1, 2), F(1), F(3, 2)]:
        linear = p_scale(kappa, lam)
        exponential = p_add(p_add({0: F(1)}, linear), p_scale(p_mul(linear, linear), F(1, 2)))
        assert exponential[1] == exponential[2] == -2*lam
        assert exponential[3] == 4*lam+4*lam*lam
        assert exponential[3]-exponential[1]*exponential[2] == 4*lam
        rows.append(dict(lambda_value=str(lam), exponential=serialize(exponential),
                         connected_coefficient=str(4*lam)))
    return dict(cross_collision_ratios=cross, own_collision_ratios=own,
                cross_connected_exponent=serialize(total_cross), kappa=serialize(kappa),
                evaluated_exponentials=rows,
                interpretation='Both H^2 fields have only internal Gaussian pairings; their site overlaps '
                               'retain charge-static covariance 4 lambda C^2. No Gaussian bridge does '
                               'not imply no charge-static correlation for the whole-field oriented mark.')


def numerical_sum(N, p, precision):
    with mp.workdps(precision):
        lam = mp.mpf(p*p)/N
        C0 = 1/lam
        cp = mp.mpf(1)
        for j in range(p):
            C0 *= (1-mp.mpf(2*j)/N)*(1+mp.mpf(2*(p-j))/N)
            cp *= mp.mpf((j+1)*(N-j))/((N//2-j)*(N//2+p-j))
        # m<p has C_m=0 exactly. All remaining charge sectors are summed.
        probability = mp.mpf(comb(N, p))/mp.power(2, N)
        current = cp
        mean_C = mean_G = second_G = shifted_G = mp.mpf(0)
        for m in range(p, N+1):
            d = current*mp.mpf(p*(N-2*m+p+1))/(N-m+1)
            mean_C += probability*current
            mean_G -= probability*d
            second_G += probability*d*d
            shifted_G -= probability*(2*m-N)*current
            if m < N:
                current *= mp.mpf((m+1)*(N-m))/((m+1-p)*(N-m+p))
                probability *= mp.mpf(N-m)/(m+1)
        assert abs(mean_G-shifted_G) < mp.mpf(10)**(-precision+10)
        variance = second_G-mean_G*mean_G
        central = []
        for x_string in PARAMETERS['central_x']:
            x = F(x_string)
            m_fraction = F(N, 2)+x*p
            assert m_fraction.denominator == 1
            m = int(m_fraction)
            ratio = mp.mpf(1)
            for j in range(p):
                ratio *= mp.mpf((m-j)*(N-m+p-j))/((N//2-j)*(N//2+p-j))
            observed = -ratio*mp.mpf(p*(N-2*m+p+1))/(N-m+1)
            target = -2*lam*(1-2*mp.mpf(x.numerator)/x.denominator)
            central.append(dict(x=x_string, m=m, observed=observed, target=target,
                                absolute_error=abs(observed-target)))
        return dict(C0_physical=C0, mean_C_over_C0=mean_C, mean_G0_over_C0=mean_G,
                    second_G0_over_C0_squared=second_G, variance_G0_over_C0_squared=variance,
                    target_mean=-2*lam, target_second=4*(lam+lam*lam), target_variance=4*lam,
                    central=central)


def high_precision_checks():
    cases = {(int(F(p*p)/F(lam)), p, 'fixed')
             for p in PARAMETERS['high_precision_p'] for lam in PARAMETERS['fixed_lambda']}
    cases.update((4*p*p+2*p, p, 'drifting') for p in PARAMETERS['high_precision_p'])
    rows = []
    with mp.workdps(max(PARAMETERS['precisions'])):
        for N, p, family in sorted(cases):
            low = numerical_sum(N, p, PARAMETERS['precisions'][0])
            high = numerical_sum(N, p, PARAMETERS['precisions'][1])
            differences = [abs(high[key]-low[key]) for key in high if key != 'central']
            differences += [abs(a[key]-b[key]) for a, b in zip(high['central'], low['central'])
                            for key in ['observed', 'target', 'absolute_error']]
            maximum = max(differences)
            assert maximum < mp.mpf('1e-58')
            record = {key: mp.nstr(value, 60) for key, value in high.items() if key != 'central'}
            record.update(N=N, p=p, family=family, lambda_N=str(F(p*p, N)),
                          all_charge_sectors=N+1, identically_zero_low_sectors=p,
                          maximum_precision_difference=mp.nstr(maximum, 12),
                          central=[{key: value if key in ['x', 'm'] else mp.nstr(value, 60)
                                    for key, value in item.items()} for item in high['central']])
            for key, target in [('mean_G0_over_C0', 'target_mean'),
                                ('second_G0_over_C0_squared', 'target_second'),
                                ('variance_G0_over_C0_squared', 'target_variance')]:
                record[key+'_relative_error'] = mp.nstr(abs(high[key]-high[target])/abs(high[target]), 30)
            rows.append(record)
            print(json.dumps(dict(N=N, p=p, lambda_N=record['lambda_N'], status='passed',
                                  variance_relative_error=record['variance_G0_over_C0_squared_relative_error'])), flush=True)
    return dict(rows=rows, systems=len(rows), precisions=PARAMETERS['precisions'],
                central_cases=sum(len(row['central']) for row in rows),
                numerical_scope='Two-precision full-binomial sums, no charge cutoff and no interval certificate. '
                                'Finite errors are retained; the limit follows from the analytic uniform-integrability argument.')


def main():
    if OUT.exists():
        raise FileExistsError(OUT)
    started = time.perf_counter()
    exact, pair_source = exact_checks(), pair_source_checks()
    high_precision = high_precision_checks()
    payload = dict(schema_version=1, status='passed', created_utc=datetime.now(timezone.utc).isoformat(),
                   parameters=PARAMETERS, parameter_sha256=digest(PARAMETERS),
                   exact=exact, pair_source=pair_source, high_precision=high_precision,
                   counts=dict(exact_systems=exact['systems'], exact_sector_identities=exact['sector_identities'],
                               cross_collision_types=4, own_collision_types=2, exponential_checks=4,
                               high_precision_systems=high_precision['systems'],
                               central_cases=high_precision['central_cases']),
                   environment=dict(python=platform.python_version(), mpmath=mp.__version__,
                                    platform=platform.platform()),
                   input_sha256={'code/check_memory_charge_background.py': sha(Path(__file__))},
                   elapsed_seconds=time.perf_counter()-started,
                   scope='Static charge background only, in beta=mu=0 normalized full-Fock trace; '
                         'C_m uses the same w_N in every charge. G2 Gaussian-chaos zero is -D C_Q. '
                         'This checks why two internal-only Gaussian fields can retain connected '
                         'charge fluctuations. It does not construct the charged AA/AB/BA kernel '
                         'or establish a general energy-plus-charge regression/pinching theorem.')
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(payload, indent=2)+'\n')
    print(json.dumps(dict(status='passed', counts=payload['counts'], output_sha256=sha(OUT),
                         elapsed_seconds=payload['elapsed_seconds'])), flush=True)


if __name__ == '__main__':
    main()
