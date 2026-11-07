"""J2-W1 mixed OTOC from an existing unequal-probe kernel, with new checks.

Imports sealed Kernel/cutoff functions only; never executes their mains.
Independent colored-chord polynomials retain distinct H and W colors even
though both have length one. Outputs are exclusive-create.
"""
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
from io import BytesIO
from itertools import combinations
from math import exp, factorial, sqrt
from pathlib import Path
import hashlib
import json
import platform
import time
import traceback

import numpy as np
import sympy as sym
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from mixed_cubic_otoc import Kernel, cutoff
from singlet_otoc import taylor_tail

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT/'results/susskind-bilocal-four-point'
TARGET = OUT/'mixed-otoc-checks.json'
FIGURE = OUT/'normalized-commutator.png'
PARAMETERS = dict(Js=1., physical_lambda=[.25, .1, .05], probe_lengths=[2, 1],
    swapped_probe_lengths=[1, 2], times=sorted({i/20 for i in range(21)} | {.01}),
    tmax=1., requested_truncation_tolerance=1e-11, refinement_increment=16,
    Taylor_orders=[0, 2, 4, 6], Taylor_comparison_times=[0., .01],
    floating_comparison_tolerance=2e-11, Taylor_comparison_tolerance=3e-12,
    figure_dpi=180,
    convention='X=J2/sqrt(w), W=W1/sqrt(C*lambda), w=C^2[4lambda+16lambda^2 q/(1+q)], '
               'C=Js^2/lambda, q=exp(-4lambda), beta=mu=0, hbar=1/string time, DS first. '
               'F=trace[X(t) W X(t) W], D=2(1-F). Distinct species, one crossed channel; mutual q^2 retained.')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def matchings(items):
    """Independent recursive enumeration, not the imported chord engine."""
    if not items:
        yield ()
        return
    first = items[0]
    for index in range(1, len(items)):
        second = items[index]
        for rest in matchings(items[1:index]+items[index+1:]):
            yield ((first, second),)+rest


def intersects(c, d):
    a, b = c; u, v = d
    return a < u < b < v or u < a < v < b


def arc_moment(arcs, d1, d2):
    word = ''.join(color+'H'*n for color, n in zip('XWXW', arcs))
    hs = tuple(i for i, color in enumerate(word) if color == 'H')
    xx = tuple(i for i, color in enumerate(word) if color == 'X')
    ww = tuple(i for i, color in enumerate(word) if color == 'W')
    powers = Counter(); diagrams = 0
    for matching in matchings(hs):
        chords = [(pair, 1) for pair in matching]+[(xx, d1), (ww, d2)]
        degree = sum(a*b for (c, a), (d, b) in combinations(chords, 2) if intersects(c, d))
        powers[degree] += 1; diagrams += 1
    return powers, diagrams


def coefficient(order, d1, d2):
    assert order % 2 == 0
    polynomial = Counter(); terms = diagrams = 0
    for a in range(order+1):
        for b in range(order-a+1):
            for c in range(order-a-b+1):
                d = order-a-b-c
                moment, count = arc_moment((a, b, c, d), d1, d2)
                multiplier = Fraction((-1)**(a+c+order//2),
                                      factorial(a)*factorial(b)*factorial(c)*factorial(d))
                for degree, value in moment.items():
                    polynomial[degree] += multiplier*value
                terms += 1; diagrams += count
    return {k: v for k, v in sorted(polynomial.items()) if v}, terms, diagrams


def exact_coefficients():
    q = sym.Symbol('q'); records = []
    for order in PARAMETERS['Taylor_orders']:
        poly, arcs, diagrams = coefficient(order, 2, 1)
        swapped, swapped_arcs, swapped_diagrams = coefficient(order, 1, 2)
        assert poly == swapped and (arcs, diagrams) == (swapped_arcs, swapped_diagrams)
        if order == 0:
            assert poly == {2: Fraction(1)}  # Mutual factor, not a fitted normalization.
        expression = sum(sym.Rational(v.numerator, v.denominator)*q**k for k, v in poly.items())
        records.append(dict(order=order, coefficients={str(k): str(v) for k, v in poly.items()},
            factored_polynomial=str(sym.factor(expression)), variable='sigma*t, sigma=Js/sqrt(lambda)',
            arc_layouts=arcs, colored_pairings=diagrams, species_swap_exact=True))
    return records


def polynomial_value(record, q):
    return sum(float(Fraction(value))*q**int(power) for power, value in record['coefficients'].items())


def bound_at(lam, t, L, d1, d2, inner):
    q = exp(-4*lam); a = q**d1; b = q**d2; mutual = q**(d1*d2)
    Cb = 1/np.prod([1-q**j for j in range(1, 2*d2)])
    K = float(Cb/(1-a)); RH = 2/sqrt(lam*(1-q))
    vacuum = 2*taylor_tail(RH*abs(t), L)
    propagation = 2*taylor_tail(RH*abs(t), L-inner)
    tail = b**inner
    absolute = 2*mutual*K*(vacuum+propagation+2*tail)
    return dict(L=L, inner=inner, lengths=[d1, d2], Cb=float(Cb), K=K, RH=RH,
        vacuum_bound=vacuum, propagator_bound=propagation, source_tail=tail,
        otoc_absolute_truncation_bound=absolute, commutator_absolute_truncation_bound=2*absolute)


def compute(coefficients):
    records = []
    for lam in PARAMETERS['physical_lambda']:
        requested = [cutoff(lam, PARAMETERS['tmax'], *ds,
                            tol=PARAMETERS['requested_truncation_tolerance']) for ds in [(2, 1), (1, 2)]]
        L = max(row['L'] for row in requested); refined = L+PARAMETERS['refinement_increment']
        base = Kernel(lam, L, 2, 1); reference = Kernel(lam, refined, 2, 1)
        swap = Kernel(lam, L, 1, 2)
        bounds = [bound_at(lam, 1., L, *ds, request['inner'])
                  for ds, request in zip([(2, 1), (1, 2)], requested)]
        bounds.append(bound_at(lam, 1., refined, 2, 1, requested[0]['inner']))
        assert max(x['otoc_absolute_truncation_bound'] for x in bounds) < PARAMETERS['requested_truncation_tolerance']
        points = []; q = exp(-4*lam)
        for t in PARAMETERS['times']:
            coarse, fine, reversed_species = base.at(t), reference.at(t), swap.at(t)
            F, D = fine['otoc'], fine['commutator_squared']
            cutoff_error = abs(coarse['otoc']-F)
            swap_error = abs(coarse['otoc']-reversed_species['otoc'])
            empirical = max(cutoff_error, swap_error)
            assert empirical < PARAMETERS['floating_comparison_tolerance']
            assert -1e-12 <= F <= 1+1e-12 and -1e-12 <= D <= 2+1e-12
            assert D == 2*(1-F)
            point_bound = bound_at(lam, t, refined, 2, 1, requested[0]['inner'])['otoc_absolute_truncation_bound']
            point = dict(time=t, otoc=F, commutator_squared=D,
                coarse_otoc=coarse['otoc'], swapped_coarse_otoc=reversed_species['otoc'],
                cutoff_absolute_difference=cutoff_error, species_swap_absolute_difference=swap_error,
                cutoff_relative_otoc_difference=cutoff_error/max(abs(F), 1e-300),
                species_swap_relative_otoc_difference=swap_error/max(abs(F), 1e-300),
                empirical_absolute_commutator_difference=2*empirical,
                empirical_relative_commutator_difference=2*empirical/max(abs(D), 1e-300),
                otoc_absolute_analytic_truncation_bound=point_bound,
                otoc_relative_analytic_truncation_bound=point_bound/max(abs(F), 1e-300),
                commutator_absolute_analytic_truncation_bound=2*point_bound,
                commutator_relative_analytic_truncation_bound=2*point_bound/max(abs(D), 1e-300),
                commutator_subtraction_condition=(1+abs(F))/max(abs(1-F), 1e-300),
                vacuum_norm_residual=max(coarse['vacuum_residual'], fine['vacuum_residual'], reversed_species['vacuum_residual']))
            if t in PARAMETERS['Taylor_comparison_times']:
                terms = [polynomial_value(c, q)*(t/sqrt(lam))**c['order'] for c in coefficients]
                series = sum(terms); error = abs(F-series)
                assert error < PARAMETERS['Taylor_comparison_tolerance']
                point.update(sixth_order_chord_series=series, chord_series_absolute_difference=error,
                    chord_series_relative_difference=error/max(abs(F), 1e-300),
                    chord_series_sum_condition=sum(abs(v) for v in terms)/max(abs(series), 1e-300),
                    chord_series_terms=terms)
            if t == 0.:
                assert abs(F-q*q) < PARAMETERS['floating_comparison_tolerance']
                point['equal_time_mutual_factor_error'] = abs(F-q*q)
            points.append(point)
        records.append(dict(lambda_value=lam, q=q, sigma=1/sqrt(lam), cutoff=L,
            refined_cutoff=refined, source_cutoff_requests=requested, interval_bounds=bounds, points=points))
        print('lambda', lam, 'L', L, 'F(1)', points[-1]['otoc'], 'D(1)', points[-1]['commutator_squared'], flush=True)
    return records


def figure_bytes(records):
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
    colors = ['#b34a3c', '#30738c', '#5d5a9d']
    for r, color in zip(records, colors):
        points = [p for p in r['points'] if p['time'] != .01]
        ax.plot([p['time'] for p in points], [p['commutator_squared'] for p in points],
                color=color, lw=2, marker='o', markersize=3, label=fr'$\lambda={r["lambda_value"]}$')
    ax.set(xlabel=r'$J_s T$', ylabel=r'$\tau([X(T),W]^\dagger[X(T),W])$',
           title=r'Mixed $J_2$-$W_1$ commutator (unit-normalized fields)', xlim=(0, 1), ylim=(0, 2.05))
    ax.grid(alpha=.2); ax.legend(frameon=False, loc='lower right')
    fig.text(.5, -.015, r'$X=J_2/\sqrt{w}$, $W=W_1/\sqrt{C\lambda}$; fixed-$\lambda$ DS limit, $\beta=\mu=0$',
             ha='center', fontsize=9)
    output = BytesIO(); fig.savefig(output, format='png', dpi=PARAMETERS['figure_dpi'], bbox_inches='tight')
    plt.close(fig)
    return output.getvalue()


def main():
    if TARGET.exists() or FIGURE.exists():
        raise FileExistsError('Scientific JSON or figure already exists')
    started = time.perf_counter(); coefficients = exact_coefficients()
    records = compute(coefficients); png = figure_bytes(records)
    inputs = ['code/check_bilocal_mixed_otoc.py', 'code/mixed_cubic_otoc.py', 'code/singlet_otoc.py',
              'code/singlet_multimark.py', 'code/singlet_coefficient.py',
              'manuscript/susskind-mixed-cubic/main.tex',
              'results/susskind-bilocal-four-point/finite-four-source-checks.json',
              'results/susskind-bilocal-four-point/joint-moment-checks.json',
              'literature/susskind-bilocal-four-point-sources.md']
    all_points = [p for r in records for p in r['points']]
    output = dict(schema_version=1, status='passed', created_utc=datetime.now(timezone.utc).isoformat(),
        parameters=PARAMETERS, parameter_sha256=digest(PARAMETERS), coefficients=coefficients, records=records,
        counts=dict(Taylor_coefficients=len(coefficients), exact_species_swap_coefficients=len(coefficients),
            arc_layouts_per_species_order=sum(c['arc_layouts'] for c in coefficients),
            colored_pairings_per_species_order=sum(c['colored_pairings'] for c in coefficients),
            physical_lambda_cases=len(records), time_points=len(all_points), kernel_evaluations=3*len(all_points),
            independent_small_time_comparisons=sum('chord_series_absolute_difference' in p for p in all_points)),
        maximum_cutoff_absolute_difference=max(p['cutoff_absolute_difference'] for p in all_points),
        maximum_species_swap_absolute_difference=max(p['species_swap_absolute_difference'] for p in all_points),
        maximum_commutator_subtraction_condition=max(p['commutator_subtraction_condition'] for p in all_points),
        figure=dict(path=FIGURE.relative_to(ROOT).as_posix(), sha256=hashlib.sha256(png).hexdigest(),
                    dpi=PARAMETERS['figure_dpi'], plotted_time_points_per_curve=21),
        input_sha256={name: sha(ROOT/name) for name in inputs},
        environment=dict(python=platform.python_version(), numpy=np.__version__, sympy=sym.__version__,
                         matplotlib=matplotlib.__version__, platform=platform.platform()),
        elapsed_seconds=time.perf_counter()-started,
        scope='Existing BINT unequal-probe positive series with d1=2,d2=1 and mutual q^2, '
              'applied after the stage three-species identification. Independent colored-chord Taylor '
              'enumeration checks the time expansion. Analytic bounds cover infinite-Jacobi truncation '
              'in exact arithmetic; refinement, species exchange and roundoff-sensitive relative errors '
              'are numerical checks, not interval certification. Only actual low-coefficient A0/B1 '
              'fields are interpreted here, not the complete finite-relative-time bilocal four-point '
              'function, a vacuum CCR algebra or a bulk scattering amplitude.')
    OUT.mkdir(parents=True, exist_ok=True)
    with FIGURE.open('xb') as stream:
        stream.write(png)
    with TARGET.open('x', encoding='utf-8') as stream:
        json.dump(output, stream, indent=2, allow_nan=False); stream.write('\n')
    print('passed', output['counts'], 'elapsed', output['elapsed_seconds'], flush=True)


if __name__ == '__main__':
    try:
        main()
    except Exception:
        diagnostic = OUT/'diagnostics'; diagnostic.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        with (diagnostic/f'mixed-otoc-failure-{stamp}.json').open('x', encoding='utf-8') as stream:
            json.dump(dict(status='failed', parameters=PARAMETERS, script_sha256=sha(Path(__file__)),
                           traceback=traceback.format_exc()), stream, indent=2)
        raise
