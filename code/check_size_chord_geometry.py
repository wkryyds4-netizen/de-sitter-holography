"""Independent chord geometry against exact local 2x2 size-mark algebra.

Only exponent() is reused from memory_size_marked_chords; its main and old
outputs are not rerun.  Geometry, Mahonian recursion, q-Hermite energy moments
and differentiated q-factorials are independently constructed here.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path
import datetime
import hashlib
import json
import math
import platform
import time

import numpy as np
import sympy as sp

from memory_size_marked_chords import exponent
from singlet_coefficient import pairings


ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results/susskind-memory-size-bridge"
q = sp.Symbol("q")


def crosses(first, second):
    # Independent endpoint-separation test; no imported crossing predicate.
    a, b = sorted(first)
    c, d = sorted(second)
    return ((a < c < b) != (a < d < b)) and not ({a, b} & {c, d})


def crossing_count(chords):
    return sum(crosses(a, b) for a, b in combinations(chords, 2))


def allowed(pairing, blocks):
    return not any(start <= a < end and start <= b < end
                   for start, end in blocks for a, b in pairing)


def geometric_exponent(pairing, blocks):
    incident = [tuple(chord for chord in pairing
                      if any(start <= endpoint < end for endpoint in chord))
                for start, end in blocks]
    result = {0: -4 * crossing_count(pairing)}
    for a, (start, end) in enumerate(blocks):
        assert len(incident[a]) == end - start
        result[1 << a] = (end - start) ** 2 - 4 * crossing_count(incident[a])
    for a, b in combinations(range(len(blocks)), 2):
        bridges = tuple(chord for chord in incident[a] if chord in incident[b])
        result[(1 << a) | (1 << b)] = len(bridges) ** 2 - 4 * crossing_count(bridges)
    return result


def two_mark_layouts(length):
    half = length // 2
    layouts = {((0, 1), (1, 2)), ((0, 1), (half, half + 1)),
               ((0, half), (half, length))}
    if length >= 4:
        layouts.update({((0, 2), (half, half + 2)),
                        ((0, 1), (length - 3, length))})
    if length >= 6:
        layouts.update({((0, 2), (2, 5)), ((1, 3), (length - 2, length))})
    return sorted(layouts)


def all_layouts():
    result = [(length, blocks) for length in (2, 4, 6, 8, 10)
              for blocks in two_mark_layouts(length)]
    result += [
        (6, ((0, 2), (2, 4), (4, 6))),
        (6, ((0, 1), (2, 3), (4, 5))),
        (6, ((0, 1), (1, 3), (4, 6))),
        (8, ((0, 2), (3, 5), (6, 8))),
        (8, ((0, 1), (3, 4), (7, 8))),
        (8, ((0, 1), (2, 4), (5, 8))),
        (8, ((0, 2), (2, 4), (4, 6), (6, 8))),
        (8, ((0, 1), (2, 3), (4, 5), (6, 7))),
    ]
    return sorted(result)


def geometry_checks(layouts):
    rows = []
    for length, blocks in layouts:
        counts = Counter()
        minimum = {"single": None, "double": None}
        maximum = {"single": None, "double": None}
        for pairing in pairings(tuple(range(length))):
            counts["all_pairings"] += 1
            if not allowed(pairing, blocks):
                continue
            counts["allowed_pairings"] += 1
            actual = exponent(length, pairing, blocks)
            expected = geometric_exponent(pairing, blocks)
            for mask in range(1 << len(blocks)):
                observed = actual.get(mask, Fraction(0))
                prediction = expected.get(mask, 0)
                assert observed == prediction, (length, blocks, pairing, mask, observed, prediction)
                counts["exact_coefficient_comparisons"] += 1
                degree = mask.bit_count()
                if degree >= 3:
                    counts["connected_degree_at_least_three_zero"] += 1
                if degree in (1, 2):
                    kind = "single" if degree == 1 else "double"
                    counts[kind + ("_negative" if prediction < 0 else
                                   "_positive" if prediction > 0 else "_zero")] += 1
                    minimum[kind] = prediction if minimum[kind] is None else min(minimum[kind], prediction)
                    maximum[kind] = prediction if maximum[kind] is None else max(maximum[kind], prediction)
        assert counts["all_pairings"] == math.prod(range(1, length, 2))
        rows.append(dict(length=length, blocks=blocks, number_marks=len(blocks),
                         counts=dict(counts), minimum_signed_connected_coefficient=minimum,
                         maximum_signed_connected_coefficient=maximum))
        print(json.dumps(dict(stage="geometry", length=length, blocks=blocks,
                              allowed_pairings=counts["allowed_pairings"], status="passed")), flush=True)
    return rows


def polynomial(counter):
    return sp.expand(sum(sp.Rational(value.numerator, value.denominator) * q ** degree
                         if isinstance(value, Fraction) else value * q ** degree
                         for degree, value in counter.items()))


def qinteger(k):
    return sum(q ** index for index in range(k))


def mahonian(j):
    # Insert the new largest permutation entry: its inversion contribution is
    # each integer 0,...,j-1 exactly once.  No chord diagram enters this recursion.
    values = [1]
    for size in range(1, j + 1):
        new = [0] * (len(values) + size - 1)
        for old_degree, count in enumerate(values):
            for extra in range(size):
                new[old_degree + extra] += count
        values = new
    return values


def qfactorial(j):
    return sp.expand(sp.prod(qinteger(k) for k in range(1, j + 1)))


def derivative_check(expression):
    assert sp.cancel(expression) == 0, expression


def bridge_moment_checks():
    rows = []
    for j in range(1, 6):
        length = 2 * j
        blocks = ((0, j), (j, length))
        zeroth, first, second_linear, second_quadratic = (Counter() for _ in range(4))
        count = 0
        for pairing in pairings(tuple(range(length))):
            if not allowed(pairing, blocks):
                continue
            count += 1
            kappa = exponent(length, pairing, blocks)
            cr = crossing_count(pairing)
            zeroth[cr] += 1
            first[cr] += kappa.get(1, 0)
            # exp(lambda*(ka ta + kb tb + kab ta tb)) at coefficient ta tb.
            second_linear[cr] += kappa.get(3, 0)
            second_quadratic[cr] += kappa.get(1, 0) * kappa.get(2, 0)
        counts = mahonian(j)
        assert count == math.factorial(j)
        assert dict(zeroth) == {index: value for index, value in enumerate(counts)}
        factorial = qfactorial(j)
        Dj = lambda f: sp.expand(j * j * f - 4 * q * sp.diff(f, q))
        d1, d2 = Dj(factorial), Dj(Dj(factorial))
        assert polynomial(zeroth) == factorial
        assert polynomial(first) == d1
        assert polynomial(second_linear) == d1
        assert polynomial(second_quadratic) == d2
        # Closed forms via the first two logarithmic derivatives of [j]_q!.
        A = sum(q / (1 - q) - k * q ** k / (1 - q ** k) for k in range(1, j + 1))
        B = sum(q / (1 - q) ** 2 - k * k * q ** k / (1 - q ** k) ** 2
                for k in range(1, j + 1))
        derivative_check(d1 - factorial * (j * j - 4 * A))
        derivative_check(d2 - factorial * ((j * j - 4 * A) ** 2 + 16 * B))
        rows.append(dict(j=j, all_bridge_pairings=count, Mahonian_coefficients=counts,
                         q_factorial=str(factorial), first_lambda_coefficient=str(d1),
                         second_lambda_coefficient=str(d1),
                         second_lambda_squared_coefficient=str(d2),
                         closed_log_derivative_form_verified=True))
    return rows


def energy_coefficients(power):
    # Monic variance-one spectral q-Hermites: x W_k = W_{k+1}+[k]_q W_{k-1}.
    result = [sp.Integer(1)]
    for _ in range(power):
        new = [sp.Integer(0)] * (len(result) + 1)
        for k, value in enumerate(result):
            new[k + 1] += value
            if k:
                new[k - 1] += value * qinteger(k)
        result = list(map(sp.expand, new))
    return result


def u1_spectral_checks():
    rows = []
    for total_energy in (0, 2, 4, 6):
        length = total_energy + 2
        for a in range(total_energy + 1):
            b = total_energy - a
            position = a + 1
            blocks = ((0, 1), (position, position + 1))
            first, second_linear, second_quadratic, zeroth = (Counter() for _ in range(4))
            count = 0
            for pairing in pairings(tuple(range(length))):
                kappa = exponent(length, pairing, blocks)
                cr = crossing_count(pairing)
                count += 1
                zeroth[cr] += 1
                first[cr] += kappa.get(1, 0)
                second_linear[cr] += kappa.get(3, 0)
                second_quadratic[cr] += kappa.get(1, 0) * kappa.get(2, 0)
            # Independent fixed-probe chord enumeration: force the two marked
            # positions to pair, then pair only the remaining Hamiltonian slots.
            bridge = Counter()
            remaining = tuple(index for index in range(length) if index not in (0, position))
            bridge_count = 0
            for rest in pairings(remaining):
                pairing = ((0, position),) + rest
                bridge[crossing_count(pairing)] += 1
                bridge_count += 1
            vacuum_moment = energy_coefficients(length)[0]
            ca, cb = energy_coefficients(a), energy_coefficients(b)
            mehler = sp.expand(sum(ca[k] * cb[k] * q ** k * qfactorial(k)
                                   for k in range(min(len(ca), len(cb)))))
            assert polynomial(zeroth) == vacuum_moment
            assert polynomial(first) == vacuum_moment
            assert polynomial(second_quadratic) == vacuum_moment
            assert polynomial(second_linear) == polynomial(bridge) == mehler
            rows.append(dict(left_energy=a, right_energy=b, length=length,
                             all_pairings=count, fixed_mark_bridge_pairings=bridge_count,
                             single_mark_lambda_coefficient=str(vacuum_moment),
                             double_mark_lambda_coefficient=str(mehler),
                             double_mark_lambda_squared_coefficient=str(vacuum_moment),
                             Mehler_parameter="r=q", qHermite_variance=1,
                             independent_vacuum_and_Mehler_recurrences_verified=True))
        print(json.dumps(dict(stage="U1_spectrum", total_energy=total_energy,
                              cases=total_energy + 1, status="passed")), flush=True)
    return rows


def main():
    destination = OUT / "geometry-checks.json"
    if destination.exists():
        raise FileExistsError(destination)
    started = time.time()
    layouts = all_layouts()
    geometry = geometry_checks(layouts)
    moments = bridge_moment_checks()
    u1 = u1_spectral_checks()
    names = ["code/check_size_chord_geometry.py", "code/memory_size_marked_chords.py",
             "code/singlet_coefficient.py", "code/singlet_multimark.py",
             "code/complex_syk_dictionary.py"]
    predecessor = "results/susskind-memory-size-bridge/marked-chord-checks.json"
    parameters = dict(layouts=layouts, full_bridge_j=list(range(1, 6)),
                      U1_total_energy=[0, 2, 4, 6],
                      arithmetic="exact integers, fractions.Fraction and symbolic rational functions",
                      normalization="C=1; q=exp(-4 lambda); lambda=p^2/N")
    encoded = json.dumps(parameters, sort_keys=True, separators=(",", ":")).encode()
    totals = Counter()
    for row in geometry:
        totals.update(row["counts"])
    data = dict(
        status="passed", created_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        parameters=parameters, parameter_sha256=hashlib.sha256(encoded).hexdigest(),
        geometry_layouts=geometry, geometry_totals=dict(totals),
        all_bridge_moments=moments, U1_two_energy_moments=u1,
        formulas=dict(single="kappa_a=j_a^2-4 cr(chords incident to block a)",
                      double="kappa_ab=b_ab^2-4 cr(chords bridging blocks a and b)",
                      connected_order_at_least_three="0",
                      full_bridge_first="lambda D_j [j]_q!",
                      full_bridge_second="lambda D_j [j]_q! + lambda^2 D_j^2 [j]_q!",
                      D_j="j^2-4 q d/dq",
                      A_j="sum_{k=1}^j [q/(1-q)-k q^k/(1-q^k)]",
                      B_j="sum_{k=1}^j [q/(1-q)^2-k^2 q^k/(1-q^k)^2]",
                      D_j_closed="[j]_q! (j^2-4 A_j)",
                      D_j_squared_closed="[j]_q! [(j^2-4 A_j)^2+16 B_j]",
                      U1_single="lambda m_(a+b+2)(q)",
                      U1_double="lambda P_ab(q;r=q) + lambda^2 m_(a+b+2)(q)"),
        scope=("Exact finite enumerations compare an independent crossing geometry with "
               "the existing local 2x2 nilpotent exponent; they do not prove arbitrary-word "
               "geometry, finite-N rank-to-size convergence, or any full-time kernel. "
               "The exponent is the pair-Poisson limiting chord rule, not a finite-N "
               "CAR equality. The connected 3/4-mark coefficients vanish in the tested "
               "layouts; the exponentiated disconnected higher-mark moments need not vanish."),
        environment=dict(python=platform.python_version(), numpy=np.__version__, sympy=sp.__version__),
        input_sha256={name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in names},
        predecessor_output_sha256=hashlib.sha256((ROOT / predecessor).read_bytes()).hexdigest(),
        exponent_cache_entries=exponent.cache_info().currsize,
        elapsed_seconds=time.time() - started)
    OUT.mkdir(parents=True, exist_ok=True)
    with destination.open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(data, indent=2) + "\n")
    print(json.dumps(dict(status="passed", layouts=len(geometry),
                          geometry_totals=dict(totals), bridge_moment_cases=len(moments),
                          U1_spectral_cases=len(u1), elapsed_seconds=data["elapsed_seconds"],
                          output_sha256=hashlib.sha256(destination.read_bytes()).hexdigest())), flush=True)


if __name__ == "__main__":
    main()
