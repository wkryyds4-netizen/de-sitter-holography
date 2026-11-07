"""Fresh fixed-degree rank-tail checks for the actual neutral-series interface.

Only this module's main is run. Historical distribution/Racah functions are
imported without executing their drivers. Every finite check uses Fraction;
large-system deficit entries are rigorous rational upper bounds obtained from
Casimir moments, not evaluations of a large full rank distribution.
"""
from datetime import datetime, timezone
from fractions import Fraction as F
from math import comb
from pathlib import Path
import hashlib
import json
import platform
import sys
import time

import mpmath as mp
import numpy as np
import sympy as sp

from cubic_linear_dynamics import crossing, weights
from memory_wick_rank import distribution, mu
from memory_sixth_moment import finite_deficit_moments
from rank_product_moments import product_distribution


ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'results/susskind-memory-neutral-marked-series/rank-checks.json'
ORDER = 6


def xrank(N, r):
    return r * (N + 1 - r)


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def raw_casimir_moments(N, m, character, order=ORDER):
    """Invert characters directly in monomials, unlike the old factorial helper."""
    assert order <= min(m, N-m)
    newton = [[F(1)]]
    for j in range(1, order+1):
        old = newton[-1]
        row = [F(0)] * (j+1)
        for ell, value in enumerate(old):
            row[ell] -= xrank(N, j-1) * value
            row[ell+1] += value
        newton.append(row)
    moments = [F(1)]
    for n in range(1, order+1):
        coefficient = F(1)
        character_polynomial = [F(0)] * (n+1)
        for j in range(n+1):
            if j:
                ell = j-1
                coefficient *= F((-n+ell)*(n-N-1+ell),
                                 (-m+ell)*(m-N+ell)*(-N-1+ell)*j)
            for ell, value in enumerate(newton[j]):
                character_polynomial[ell] += coefficient * value
        leading = character_polynomial[n]
        assert leading != 0
        moments.append((character(n) - sum(character_polynomial[j]*moments[j]
                                          for j in range(n))) / leading)
    return moments


def shifted_moments(raw, top):
    return [sum(F(comb(k, j)*(-1)**j*top**(k-j))*raw[j]
                for j in range(k+1)) for k in range(len(raw))]


def transition_table(N, p, m):
    R = min(m, N-m)
    dims = [choose(N, r)**2 - choose(N, r-1)**2 for r in range(R+1)]
    dimension = choose(N, m)
    characters = [[crossing(N, m, r, s) for s in range(R+1)]
                  for r in range(R+1)]
    eigenvalues = [mu(N, p, m, s) for s in range(R+1)]
    result = []
    for r in range(R+1):
        row = [F(dims[t], dimension**2) *
               sum(dims[s]*characters[r][s]*eigenvalues[s]*characters[t][s]
                   for s in range(R+1)) for t in range(R+1)]
        assert sum(row) == 1 and min(row) >= 0
        assert all(value == 0 for t, value in enumerate(row) if abs(t-r) > p)
        result.append(row)
    input_weights = weights(N, m, p)
    assert result[0] == input_weights + [F(0)]*(R+1-len(input_weights))
    return result


def finite_checks():
    systems = [(14, 1, 6), (16, 2, 7), (20, 2, 10), (24, 2, 11)]
    chain_rows, step_rows, mixture_rows = [], [], []
    exact_raw_equalities = 0
    for N, p, m in systems:
        transition = transition_table(N, p, m)
        R = min(m, N-m)
        for r, row in enumerate(transition):
            raw = raw_casimir_moments(N, m,
                                      lambda n: crossing(N, m, r, n)*mu(N, p, m, n))
            direct_raw = [sum(value*xrank(N, t)**j for t, value in enumerate(row))
                          for j in range(ORDER+1)]
            assert raw == direct_raw
            exact_raw_equalities += ORDER+1
            actual = [sum(value*(p+r-t)**j for t, value in enumerate(row))
                      for j in range(ORDER+1)]
            assert all(value >= 0 for value in actual)
            delta_x = shifted_moments(raw, xrank(N, r+p))
            denominator = N+1-2*(r+p)
            upper = None
            if denominator > 0:
                upper = {str(j): delta_x[j]/denominator**j for j in [4, 6]}
                assert all(actual[j] <= upper[str(j)] for j in [4, 6])
            step_rows.append(dict(N=N, p=p, m=m, input_rank=r,
                                  delta_fourth=str(actual[4]), delta_sixth=str(actual[6]),
                                  casimir_comparison_denominator=denominator,
                                  casimir_deficit_upper_bounds=None if upper is None else
                                  {j: str(v) for j, v in upper.items()},
                                  probabilities=list(map(str, row))))
        previous = [F(1)] + [F(0)]*R
        for k in range(1, 7):
            current = [sum(previous[r]*transition[r][t] for r in range(R+1))
                       for t in range(R+1)]
            if k in [3, 4, 6]:
                direct = distribution(N, p, m, k)
                assert direct == current
                raw = raw_casimir_moments(N, m, lambda n: mu(N, p, m, n)**k)
                assert raw == [sum(value*xrank(N, r)**j for r, value in enumerate(direct))
                               for j in range(ORDER+1)]
                exact_raw_equalities += ORDER+1
                moments = [sum(value*(k*p-r)**j for r, value in enumerate(direct))
                           for j in range(ORDER+1)]
                assert all(value >= 0 for value in moments)
                # Jensen bound for the actual path, keeping each earlier law.
                path = [F(1)] + [F(0)]*R
                total = F(0)
                conditional_sixth = [sum(value*(p+r-t)**6 for t, value in enumerate(line))
                                     for r, line in enumerate(transition)]
                for _ in range(k):
                    total += sum(path[r]*conditional_sixth[r] for r in range(R+1))
                    path = [sum(path[r]*transition[r][t] for r in range(R+1))
                            for t in range(R+1)]
                assert moments[6] <= k**5 * total
                chain_rows.append(dict(N=N, p=p, m=m, k=k,
                                       actual_deficit_moments=list(map(str, moments)),
                                       sixth_path_jensen_upper_bound=str(k**5*total),
                                       probabilities=list(map(str, direct))))
            previous = current
        # A separate pure-rank product mixture checks the conditional kernel.
        if N in [14, 20]:
            omega = weights(N, m, p)
            for r in [0, 1, 3]:
                mixture = [F(0)]*(R+1)
                for a, weight in enumerate(omega):
                    pure = product_distribution(N, m, r, a)
                    assert len(pure) == len(mixture)
                    mixture = [old+weight*value for old, value in zip(mixture, pure)]
                    raw = raw_casimir_moments(N, m,
                                             lambda n: crossing(N, m, r, n)*crossing(N, m, a, n))
                    scaled = [value/F(N**j) for j, value in
                              enumerate(shifted_moments(raw, xrank(N, r+a)))]
                    assert scaled == finite_deficit_moments(N, m, r, a, ORDER)
                assert mixture == transition[r]
                mixture_rows.append(dict(N=N, p=p, m=m, input_rank=r,
                                         complete_input_ranks=list(range(p+1))))
        print(f'finite rank system N={N}, p={p}, m={m} passed', flush=True)
    return dict(chain_rows=chain_rows, conditional_rows=step_rows,
                independent_product_mixtures=mixture_rows,
                raw_casimir_exact_equalities=exact_raw_equalities)


def large_checks():
    chain_rows, step_rows = [], []
    for lam in [F(1, 4), F(1, 2)]:
        for p in [16, 32, 64]:
            N = int(p*p/lam)
            for x in [-1, 0, 1]:
                m = N//2 + x*p
                for k in [3, 4, 6]:
                    raw = raw_casimir_moments(N, m, lambda n: mu(N, p, m, n)**k)
                    delta_x = shifted_moments(raw, xrank(N, k*p))
                    denominator = N+1-2*k*p
                    assert denominator > 0
                    upper = {j: delta_x[j]/denominator**j for j in [4, 6]}
                    assert all(value >= 0 for value in upper.values())
                    chain_rows.append(dict(N=N, p=p, m=m, x=x, k=k, lambda_value=str(lam),
                                           delta_casimir_fourth_over_N4=str(delta_x[4]/N**4),
                                           delta_casimir_sixth_over_N6=str(delta_x[6]/N**6),
                                           deficit_fourth_upper_bound=str(upper[4]),
                                           deficit_sixth_upper_bound=str(upper[6]),
                                           sixth_upper_decimal=float(upper[6])))
                for r in [0, p, 2*p, 4*p]:
                    raw = raw_casimir_moments(N, m,
                                             lambda n: crossing(N, m, r, n)*mu(N, p, m, n))
                    delta_x = shifted_moments(raw, xrank(N, r+p))
                    denominator = N+1-2*(r+p)
                    assert denominator > 0
                    upper = {j: delta_x[j]/denominator**j for j in [4, 6]}
                    assert all(value >= 0 for value in upper.values())
                    step_rows.append(dict(N=N, p=p, m=m, x=x, input_rank=r,
                                          lambda_value=str(lam),
                                          delta_fourth_upper_bound=str(upper[4]),
                                          delta_sixth_upper_bound=str(upper[6]),
                                          sixth_upper_decimal=float(upper[6])))
    return dict(chain_proxy_rows=chain_rows, conditional_proxy_rows=step_rows,
                scope='Exact rational Casimir moments and positive-deficit comparison bounds. '
                      'These are not direct large-system rank-deficit moments or a Poisson-law test.')


def mu_product(N, p, m, r):
    value = F(1)
    for j in range(p):
        value *= max(F(0), 1-F(r, m-j))
        value *= max(F(0), 1-F(r, N-m+p-j))
    return value


def slope_checks():
    rows, changes = [], []
    for lam in [F(1, 4), F(1, 2)]:
        for p in [16, 32, 64]:
            N = int(p*p/lam)
            for x in [-1, 0, 1]:
                m = N//2 + x*p
                assert p <= N//8 and N//4 <= m <= 3*N//4
                for k in [1, 2, 3, 4, 6]:
                    base = mu(N, p, m, k*p)
                    assert base == mu_product(N, p, m, k*p)
                    for deficit in [0, 1, 2, 3, 6]:
                        value = mu(N, p, m, k*p-deficit)
                        assert value == mu_product(N, p, m, k*p-deficit)
                        scaled = p*(value-base)
                        upper = 12*lam*deficit
                        assert 0 <= scaled <= upper
                        summaries = []
                        for precision in [80, 120]:
                            with mp.workdps(precision):
                                ll = mp.mpf(lam.numerator)/lam.denominator
                                exact = mp.mpf(scaled.numerator)/scaled.denominator
                                target = 4*ll*mp.exp(-4*ll*k)*deficit
                                error = exact/target-1 if deficit else exact-target
                                summaries.append([mp.nstr(v, precision-2)
                                                  for v in [exact, target, error]])
                        with mp.workdps(120):
                            change = max(abs(mp.mpf(a)-mp.mpf(b))/(1+abs(mp.mpf(b)))
                                         for a, b in zip(*summaries))
                            assert change < mp.mpf('1e-70')
                            changes.append(change)
                            values = [mp.nstr(mp.mpf(v), 65) for v in summaries[-1]]
                        rows.append(dict(N=N, p=p, m=m, x=x, k=k, deficit=deficit,
                                         lambda_value=str(lam), exact_scaled_difference=str(scaled),
                                         exact_domination_upper=str(upper), finite_value=values[0],
                                         limiting_slope=values[1],
                                         relative_limit_deviation=values[2] if deficit else None,
                                         zero_deficit_absolute_error=values[2] if not deficit else None))
    with mp.workdps(120):
        largest = mp.nstr(max(changes), 35)
    return dict(rows=rows, precisions=[80, 120], max_scaled_precision_change=largest,
                scope='Domination and positive-part/binomial agreement are exact rational checks. '
                      'Exponential targets use precision refinement, not interval certification; '
                      'finite-size deviations are retained without a convergence tolerance.')


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
    finite = finite_checks()
    large = large_checks()
    slopes = slope_checks()
    result = dict(status='passed', created_at_utc=datetime.now(timezone.utc).isoformat(),
                  finite=finite, large=large, slopes=slopes,
                  counts=dict(finite_chain_sixth_cases=len(finite['chain_rows']),
                              finite_complete_conditional_distributions=len(finite['conditional_rows']),
                              independent_product_mixtures=len(finite['independent_product_mixtures']),
                              large_chain_proxy_bounds=len(large['chain_proxy_rows']),
                              large_conditional_proxy_bounds=len(large['conditional_proxy_rows']),
                              exact_slope_domination_cases=len(slopes['rows'])),
                  environment=dict(python=platform.python_version(), numpy=np.__version__,
                                   sympy=sp.__version__, mpmath=mp.__version__),
                  input_sha256=source_hashes(), elapsed_seconds=time.time()-started,
                  scope='Finite independent-copy positive rank laws and one-step transitions. '
                        'The fixed-k high-moment theorem uses the existing polynomial-ring '
                        'Racah cancellations on bounded charge/rank parameter sets; finite cases '
                        'do not prove its uniformity. This driver does not prove an all-k sixth '
                        'moment rate, a Poisson law, actual-exponential convergence, or a mixed kernel.',
                  reproducibility='Deterministic exact arithmetic; no random seed or rank cutoff. '
                                  'No historical main function is executed. Output is exclusive-create.')
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps(dict(status=result['status'], counts=result['counts'],
                          output_sha256=hashlib.sha256(OUT.read_bytes()).hexdigest(),
                          elapsed_seconds=result['elapsed_seconds'])), flush=True)


if __name__ == '__main__':
    main()
