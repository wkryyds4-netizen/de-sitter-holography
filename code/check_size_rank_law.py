"""Fresh exact neutral-size algebra and candidate defect-law checks.

The finite neutral-superoperator identities and the analytic candidate law are
tested independently. No identification with an actual finite-N rank MGF is
assumed. Probability tails have explicit exact-rational upper bounds.
"""
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction as F
from functools import reduce
from itertools import permutations
from math import comb, factorial, gcd, lcm
from pathlib import Path
import hashlib
import json
import platform
import sys
import time

import mpmath as mp
import numpy as np

from complex_syk_dictionary import fermions


ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'results/susskind-memory-size-bridge/rank-law-checks.json'
MAX_INTEGER_PRODUCT_BOUND = 0
POWER_CHECK = 256
TAIL_TARGET = F(1, 10**42)


def integer_matrix(value):
    value = np.asarray(value)
    if np.iscomplexobj(value):
        assert not np.any(value.imag)
        value = value.real
    assert np.array_equal(value, np.rint(value))
    assert np.max(abs(value)) < 2**48
    return value.astype(np.int64)


def imul(left, right):
    global MAX_INTEGER_PRODUCT_BOUND
    bound = left.shape[1]*int(np.max(abs(left)))*int(np.max(abs(right)))
    assert bound < 2**60
    MAX_INTEGER_PRODUCT_BOUND = max(MAX_INTEGER_PRODUCT_BOUND, bound)
    return left @ right


def projections(operator, roots):
    identity = np.eye(len(operator), dtype=np.int64)
    result = []
    for root in roots:
        numerator = identity.copy()
        denominator = 1
        for other in roots:
            if other != root:
                numerator = imul(numerator, operator-other*identity)
                denominator *= root-other
        if denominator < 0:
            numerator, denominator = -numerator, -denominator
        common = reduce(gcd, [denominator]+[abs(int(value)) for value in numerator.flat])
        numerator //= common
        denominator //= common
        assert np.array_equal(imul(numerator, numerator), denominator*numerator)
        assert np.array_equal(imul(operator, numerator), root*numerator)
        assert np.array_equal(numerator, numerator.T)
        result.append((numerator, denominator))
    common = lcm(*[denominator for _, denominator in result])
    total = sum(numerator*(common//denominator) for numerator, denominator in result)
    assert np.array_equal(total, common*identity)
    for i, (left, _) in enumerate(result):
        for right, _ in result[:i]:
            assert not np.any(imul(left, right))
    return result


def neutral_superoperators():
    rows = []
    for N in [2, 3, 4]:
        cs = [integer_matrix(c) for c in fermions(N)]
        dimension = 1 << N
        identity = np.eye(dimension, dtype=np.int64)
        charges = [state.bit_count() for state in range(dimension)]
        Q = np.diag(charges)
        pairs = [(a, b) for a in range(dimension) for b in range(dimension)
                 if charges[a] == charges[b]]
        indices = [a*dimension+b for a, b in pairs]
        restrict = lambda matrix: matrix[np.ix_(indices, indices)]
        D = restrict(np.kron(Q, identity)-sum(np.kron(c.T, c.T) for c in cs))
        Db = restrict(np.kron(N*identity-Q, identity)-sum(np.kron(c, c) for c in cs))
        L2, T2 = D+Db, D-Db
        assert np.array_equal(L2, L2.T)
        Qleft = restrict(np.kron(Q, identity))
        assert np.array_equal(T2.T, 4*Qleft-2*N*np.eye(len(pairs), dtype=np.int64)-T2)
        generators = {(a, b): imul(cs[a].T, cs[b]) for a in range(N) for b in range(N)}
        adjoint = {(a, b): restrict(np.kron(value, identity)-np.kron(identity, value.T))
                   for (a, b), value in generators.items()}
        C2 = sum(imul(adjoint[a, b], adjoint[b, a]) for a in range(N) for b in range(N))
        assert np.array_equal(C2, C2.T)
        assert np.array_equal(2*C2, 2*(N+1)*L2-imul(L2, L2)-imul(T2.T, T2))
        assert np.array_equal(imul(L2, C2), imul(C2, L2))
        assert np.array_equal(imul(T2, L2)-imul(L2, T2), 2*T2)
        rank_projectors = projections(C2, [2*r*(N+1-r) for r in range(N//2+1)])
        size_projectors = projections(L2, [2*ell for ell in range(N+1)])
        joint = []
        for r, (Pr, dr) in enumerate(rank_projectors):
            irreducible_dimension = comb(N, r)**2 - (comb(N, r-1)**2 if r else 0)
            assert F(int(np.trace(Pr)), dr) == (N-2*r+1)*irreducible_dimension
            for ell, (Pl, dl) in enumerate(size_projectors):
                product = imul(Pr, Pl)
                assert np.array_equal(product, imul(Pl, Pr))
                multiplicity = F(int(np.trace(product)), dr*dl)
                expected = irreducible_dimension if r <= ell <= N-r else 0
                assert multiplicity == expected
                if expected == 0:
                    assert not np.any(product)
                joint.append(dict(rank=r, half_Clifford_size=ell,
                                  auxiliary_spin=str(F(N, 2)-r),
                                  multiplicity=int(multiplicity)))
        for ell, (Pl, denominator) in enumerate(size_projectors):
            assert F(int(np.trace(Pl)), denominator) == comb(N, ell)**2
        rows.append(dict(N=N, complete_neutral_matrix_units=len(pairs),
                         exact_superoperator_identities=5, joint_multiplicities=joint,
                         rank_less_equal_half_size=True))
        print(f'complete neutral superoperators N={N} passed', flush=True)
    return dict(rows=rows,
                conventions='Full normalized Fock Hilbert-Schmidt inner product, restricted '
                            'to all charge-neutral matrix units. L=(D+Db)/2, T=(D-Db)/2, '
                            'Cx has rank eigenvalue r(N+1-r). D/Db retain neighboring-charge '
                            'blocks; the adjoint is the full-neutral HS adjoint.',
                identities=['Cx=L(N+1-L)-T_dagger*T', '[L,Cx]=0', '[T,L]=T',
                            'L_dagger=L', 'T_dagger=(2Q-N)-T'],
                spin_scope='On rank r the auxiliary spin is J=N/2-r, while L has eigenvalues '
                           'r,r+1,...,N-r. Thus R=N/2-J and R<=L. J is a spin label, '
                           'not the operator L itself.')


def mahonian(j):
    coefficients = [1]
    for length in range(1, j+1):
        following = [0]*(len(coefficients)+length-1)
        for degree, value in enumerate(coefficients):
            for increment in range(length):
                following[degree+increment] += value
        coefficients = following
    return coefficients


def integer_weight_checks():
    records, powers = [], {}
    for j in range(1, 7):
        counts = mahonian(j)
        inversions = [sum(perm[a] > perm[b] for a in range(j) for b in range(a+1, j))
                      for perm in permutations(range(j))]
        direct = Counter(inversions)
        assert counts == [direct[index] for index in range(len(counts))]
        assert counts == counts[::-1] and sum(counts) == factorial(j)
        top = j*(j-1)//2
        for index in range(top+1):
            assert (j*j-4*index)+(j*j-4*(top-index)) == 2*j
        values = []
        for degree in range(POWER_CHECK+1):
            value = sum(count*(j*j-4*index)**degree for index, count in enumerate(counts))
            independent = sum((j*j-4*index)**degree for index in inversions)
            assert value == independent
            assert value >= 0
            paired = 0
            for index in range((top+1)//2):
                a, b = j*j-4*index, j*j-4*(top-index)
                pair = a**degree+b**degree
                assert pair >= 0
                paired += counts[index]*pair
            if top % 2 == 0:
                paired += counts[top//2]*j**degree
            assert paired == value
            values.append(value)
        powers[j] = values
        records.append(dict(j=j, permutation_count=factorial(j), Mahonian_counts=counts,
                            tested_power_range=[0, POWER_CHECK], all_power_sums_nonnegative=True,
                            first_eight_power_sums=values[:8],
                            largest_tested_power_sum=str(values[-1])))
    return dict(rows=records,
                analytic_positivity='Mahonian reflection pairs a=j^2-4*inv with 2j-a. '
                                    'Writing the pair j+y,j-y gives '
                                    '(j+y)^d+(j-y)^d=2*sum_even_l binom(d,l)*j^(d-l)*y^l>=0. '
                                    'The unpaired central term, if present, is j^d. '
                                    'This is an all-d analytic argument; enumeration checks only d<=256.'), powers


def exp_tail_rational(rate, cutoff):
    """Exact bound for sum_{d>cutoff} rate^d/d!, rate>=0 rational."""
    assert rate >= 0 and cutoff >= 0 and rate < cutoff+2
    return rate**(cutoff+1)/factorial(cutoff+1)/(1-rate/F(cutoff+2))


def cutoff_for(j, rate):
    cutoff = max(8, int(rate)+1)
    while factorial(j)*exp_tail_rational(rate, cutoff) > TAIL_TARGET:
        cutoff += 1
    assert cutoff <= POWER_CHECK
    return cutoff


def mp_fraction(value):
    return mp.mpf(value.numerator)/value.denominator


def qfactorial(j, q):
    return mp.fprod(mp.fsum(q**index for index in range(length))
                    for length in range(1, j+1))


def evaluate_case(j, lam, powers, precision):
    with mp.workdps(precision):
        ll = mp_fraction(lam)
        q = mp.exp(-4*ll)
        Z = qfactorial(j, q)
        counts = mahonian(j)
        partition = mp.fsum(count*q**index for index, count in enumerate(counts))
        assert abs(Z-partition) <= mp.mpf('1e-'+str(precision-10))*max(1, abs(Z))
        inv_mean = mp.fsum(count*index*q**index for index, count in enumerate(counts))/Z
        inv_second = mp.fsum(count*index**2*q**index for index, count in enumerate(counts))/Z
        inv_var = inv_second-inv_mean**2
        closed_inv_mean = j*q/(1-q)-mp.fsum(length*q**length/(1-q**length)
                                           for length in range(1, j+1))
        closed_inv_var = j*q/(1-q)**2-mp.fsum(length**2*q**length/(1-q**length)**2
                                             for length in range(1, j+1))
        tolerance = mp.mpf('1e-'+str(precision-10))
        assert abs(inv_mean-closed_inv_mean) < tolerance
        assert abs(inv_var-closed_inv_var) < tolerance
        mean = ll*(j*j-4*inv_mean)
        variance = mean+16*ll**2*inv_var
        second = variance+mean**2
        direct_a2 = mp.fsum(count*(j*j-4*index)**2*q**index
                            for index, count in enumerate(counts))/Z
        assert abs(second-(mean+ll**2*direct_a2)) < tolerance*max(1, abs(second))
        mgf = []
        normalization = None
        for z in [F(1, 2), F(1), F(3, 2)]:
            zz = mp_fraction(z)
            rate = lam*j*j*z
            cutoff = cutoff_for(j, rate)
            raw_weights = [mp.exp(-ll*j*j)*ll**degree*values/factorial(degree)/Z
                           for degree, values in enumerate(powers[:cutoff+1])]
            assert min(raw_weights) >= 0
            partial = mp.fsum(weight*zz**degree for degree, weight in enumerate(raw_weights))
            exact_pgf = mp.exp(ll*j*j*(zz-1))*qfactorial(j, mp.exp(-4*ll*zz))/Z
            independent_pgf = mp.exp(-ll*j*j)*mp.fsum(
                count*mp.exp(ll*(j*j-4*index)*zz) for index, count in enumerate(counts))/Z
            assert abs(exact_pgf-independent_pgf) < tolerance*max(1, abs(exact_pgf))
            rational_bound = factorial(j)*exp_tail_rational(rate, cutoff)
            missing = exact_pgf-partial
            assert missing >= -tolerance*max(1, abs(exact_pgf))
            assert missing <= mp_fraction(rational_bound)+tolerance*max(1, abs(exact_pgf))
            row = dict(z=str(z), t='log('+str(z)+')', cutoff=cutoff,
                       normalized_MGF=mp.nstr(exact_pgf, precision-2),
                       positive_partial_sum=mp.nstr(partial, precision-2),
                       numerical_missing_tail=mp.nstr(missing, precision-2),
                       exact_rational_tail_upper=str(rational_bound),
                       tail_upper_decimal=mp.nstr(mp_fraction(rational_bound), 45))
            mgf.append(row)
            if z == 1:
                partial_mean = mp.fsum(degree*weight for degree, weight in enumerate(raw_weights))
                partial_second = mp.fsum(degree**2*weight for degree, weight in enumerate(raw_weights))
                mean_tail = factorial(j)*rate*exp_tail_rational(rate, cutoff-1)
                second_tail = factorial(j)*(rate**2*exp_tail_rational(rate, cutoff-2)
                                             +rate*exp_tail_rational(rate, cutoff-1))
                assert -tolerance <= mean-partial_mean <= mp_fraction(mean_tail)+tolerance
                assert -tolerance <= second-partial_second <= mp_fraction(second_tail)+tolerance
                normalization = dict(cutoff=cutoff,
                                     pmf=[mp.nstr(weight, 65) for weight in raw_weights],
                                     positive_partial_mean=mp.nstr(partial_mean, precision-2),
                                     positive_partial_second_moment=mp.nstr(partial_second, precision-2),
                                     missing_first_moment=mp.nstr(mean-partial_mean, precision-2),
                                     missing_second_moment=mp.nstr(second-partial_second, precision-2),
                                     exact_rational_first_moment_tail_upper=str(mean_tail),
                                     exact_rational_second_moment_tail_upper=str(second_tail))
        return dict(j=j, lambda_value=str(lam), precision=precision,
                    unnormalized_total=mp.nstr(Z, precision-2),
                    inversion_mean=mp.nstr(inv_mean, precision-2),
                    inversion_variance=mp.nstr(inv_var, precision-2),
                    mean=mp.nstr(mean, precision-2), variance=mp.nstr(variance, precision-2),
                    second_moment=mp.nstr(second, precision-2),
                    MGF=mgf, normalization=normalization)


def probability_checks(power_table):
    rows, changes = [], []
    for j in range(1, 7):
        for lam in [F(1, 10), F(1, 4), F(1)]:
            lo = evaluate_case(j, lam, power_table[j], 80)
            hi = evaluate_case(j, lam, power_table[j], 120)
            with mp.workdps(120):
                fields = ['unnormalized_total', 'mean', 'variance', 'second_moment']
                pairs = [(lo[name], hi[name]) for name in fields]
                pairs += [(a[name], b[name]) for a, b in zip(lo['MGF'], hi['MGF'])
                          for name in ['normalized_MGF', 'positive_partial_sum']]
                change = max(abs(mp.mpf(a)-mp.mpf(b))/(1+abs(mp.mpf(b))) for a, b in pairs)
                assert change < mp.mpf('1e-70')
                changes.append(change)
                hi['scaled_80_120_precision_change'] = mp.nstr(change, 35)
            rows.append(hi)
        print(f'positive candidate law j={j} passed', flush=True)
    with mp.workdps(120):
        maximum = mp.nstr(max(changes), 35)
    return dict(rows=rows, precisions=[80, 120], max_scaled_precision_change=maximum,
                formulas=dict(weight='w_(j,d)=exp(-lambda*j^2)*lambda^d/d!*'
                                     'sum_pi(j^2-4*inv(pi))^d',
                              total='Z_j=[j]_(exp(-4*lambda))!',
                              normalized_MGF='M_j(t)=exp(lambda*j^2*(exp(t)-1))*'
                                             '[j]_(exp(-4*lambda*exp(t)))!/Z_j',
                              mean='lambda*(j^2-4*E_q[inv])',
                              variance='mean+16*lambda^2*Var_q(inv)'),
                rigorous_tail='For z=exp(t)>0 rational, |j^2-4 inv|<=j^2 and '
                              'exp(-lambda*j^2)/Z_j<=1 give the exact rational bound '
                              'j!*(lambda*j^2*z)^(D+1)/(D+1)!/'
                              '(1-lambda*j^2*z/(D+2)). First/second moment tails use '
                              'd/d!=1/(d-1)! and d^2/d!=1/(d-2)!+1/(d-1)!.',
                arithmetic_scope='The rational tail bounds are exact. Exponentials, MGF '
                                 'differences and displayed probabilities are 80/120-digit '
                                 'refinement checks, not interval-rounded evaluations. '
                                 'The normalized MGF includes division by Z_j.')


def source_hashes():
    paths = {Path(__file__).resolve()}
    for module in tuple(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename:
            path = Path(filename).resolve()
            if path.parent == ROOT/'code' and path.suffix == '.py':
                paths.add(path)
    return {path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(paths)}


def main():
    if OUT.exists():
        raise FileExistsError(OUT)
    started = time.time()
    algebra = neutral_superoperators()
    weights, table = integer_weight_checks()
    probability = probability_checks(table)
    counts = dict(neutral_matrix_units=sum(row['complete_neutral_matrix_units'] for row in algebra['rows']),
                  superoperator_systems=len(algebra['rows']),
                  exact_integer_weight_powers=6*(POWER_CHECK+1),
                  probability_cases=len(probability['rows']),
                  MGF_cases=sum(len(row['MGF']) for row in probability['rows']))
    data = dict(status='passed', created_at_utc=datetime.now(timezone.utc).isoformat(),
                neutral_algebra=algebra, integer_weights=weights, probabilities=probability,
                counts=counts, largest_integer_product_bound=MAX_INTEGER_PRODUCT_BOUND,
                environment=dict(python=platform.python_version(), numpy=np.__version__, mpmath=mp.__version__),
                input_sha256=source_hashes(), elapsed_seconds=time.time()-started,
                scope='Exact finite full-neutral superoperator identities and an independently '
                      'normalized analytic candidate defect law. These checks do not establish '
                      'convergence of the actual finite-N rank distribution or its MGF, '
                      'nor a finite-time mixed memory kernel.',
                reproducibility='Deterministic enumeration and arithmetic, no random seed. '
                                'Only this fresh driver main is executed; output is exclusive-create.')
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open('x', encoding='utf-8') as handle:
        json.dump(data, handle, indent=2)
        handle.write('\n')
    print(json.dumps(dict(status='passed', counts=counts,
                          output_sha256=hashlib.sha256(OUT.read_bytes()).hexdigest(),
                          elapsed_seconds=data['elapsed_seconds'])), flush=True)


if __name__ == '__main__':
    main()
