"""Fresh fixed-degree four-source checks, with exclusive-create output.

Explicit integer full-Fock CAR matrices are independent of the local source
engine. Independent Fraction 2x2 sources and partition coefficients audit the
sealed exponent helper; no historical scientific main is executed.
"""
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction as Fr
from itertools import combinations, product
from math import comb, factorial
from pathlib import Path
import hashlib
import json
import platform
import time

import numpy as np

from memory_body_source_chords import exponent as sealed_exponent
from memory_body_source_chords import local as sealed_local


ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results/susskind-bilocal-four-point/finite-four-source-checks.json"
PARAMETERS = {
    "systems": [[2, 1], [3, 1]],
    "double_slot_fields": [[0, 2], [2, 4], [4, 6], [6, 8]],
    "double_slot_kinds": ["DDDD", "LLLL", "DLDL", "LDLD"],
    "double_slot_assignments": {
        "2": [[0]*8, [0, 1]*4, [0, 0, 1, 1]*2, [0, 1, 1, 0, 1, 0, 0, 1]],
        "3": [[0]*8, [0, 1, 2, 0, 1, 2, 0, 1],
              [0, 0, 1, 1, 2, 2, 0, 1], [0, 1, 1, 2, 2, 0, 1, 0]],
    },
    "single_slot_fields": [[0, 1], [1, 2], [2, 3], [3, 4]],
    "mixed_slot_fields": [[0, 1], [1, 3], [3, 4], [4, 6]],
    "single_slot_assignments": "all N^4, LLLL, both raw and W1 subtraction",
    "W1_counterterm": "(p-L-lambda_N)H, lambda_N=p^2/N",
    "geometry": "all allowed double-slot x four kinds; single-slot x all 16 kinds; mixed 1/2/1/2 LDLD",
    "all_source_masks": list(range(16)),
    "formal_lambda_for_partition_check": "1/4",
    "arithmetic": "bounded int64 full CAR; Fraction local algebra and counterterms",
    "scope": "fixed p-set diagram coefficients, not finite-time four-point convergence",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def pairings(slots):
    if not slots:
        yield ()
        return
    for i in range(1, len(slots)):
        for rest in pairings(slots[1:i] + slots[i+1:]):
            yield ((slots[0], slots[i]),) + rest


def legal(pairing, fields):
    return not any(a <= x < y < b for a, b in fields for x, y in pairing)


def pmul(first, second):
    out = Counter()
    for a, x in first.items():
        for b, y in second.items():
            if not a & b:
                out[a | b] += x*y
    return {a: x for a, x in out.items() if x}


def pinverse(poly):
    value = poly[0]
    rest = {a: -x/value for a, x in poly.items() if a}
    out, power = Counter({0: Fr(1)}), {0: Fr(1)}
    for _ in range(4):
        power = pmul(power, rest)
        out.update(power)
    return {a: x/value for a, x in out.items() if x}


def mm(a, b):
    return tuple(sum(a[2*i+k]*b[2*k+j] for k in range(2))
                 for i in range(2) for j in range(2))


IDENTITY = tuple(map(Fr, (1, 0, 0, 1)))
ANNIHILATOR = tuple(map(Fr, (0, 1, 0, 0)))
CREATOR = tuple(map(Fr, (0, 0, 1, 0)))
NUMBER = mm(CREATOR, ANNIHILATOR)
ZERO = (Fr(0),)*4


def local_loss(matrix, kind):
    # Definitions from anticommutator and creation/annihilation, not a
    # component-wise copy of the sealed helper.
    forward = mm(CREATOR, mm(matrix, ANNIHILATOR))
    if kind == "D":
        left, right = mm(NUMBER, matrix), mm(matrix, NUMBER)
        return tuple((a+b)/2-c for a, b, c in zip(left, right, forward))
    reverse = mm(ANNIHILATOR, mm(matrix, CREATOR))
    return tuple((a-b-c)/2 for a, b, c in zip(matrix, forward, reverse))


def local_source(length, pi, blocks, occupied):
    begins = {a: (index, b, kind) for index, (a, b, kind) in enumerate(blocks)}
    result, pos = {0: IDENTITY}, 0
    while pos < length:
        mark, end, kind = begins.get(pos, (None, pos+1, None))
        value = IDENTITY
        for slot in range(pos, end):
            if slot in occupied:
                value = mm(value, CREATOR)
            if pi[slot] in occupied:
                value = mm(value, ANNIHILATOR)
        factor = {0: value}
        if mark is not None:
            factor[1 << mark] = tuple(-x for x in local_loss(value, kind))
        updated = {}
        for a, x in result.items():
            for b, y in factor.items():
                if not a & b:
                    term = mm(x, y)
                    updated[a | b] = tuple(u+v for u, v in zip(updated.get(a | b, ZERO), term))
        result, pos = updated, end
    return {a: (x[0]+x[3])/2 for a, x in result.items() if x[0]+x[3]}


def sign(sequence):
    return -1 if sum(a > b for i, a in enumerate(sequence) for b in sequence[i+1:]) % 2 else 1


def car_system(n):
    dim = 1 << n
    eye = np.eye(dim, dtype=np.int64)
    annihilators = []
    for i in range(n):
        c = np.zeros_like(eye)
        for state in range(dim):
            if state & (1 << i):
                c[state ^ (1 << i), state] = (-1)**((state & ((1 << i)-1)).bit_count())
        annihilators.append(c)
    number = np.diag([state.bit_count() for state in range(dim)]).astype(np.int64)
    return eye, annihilators, number


def imul(a, b):
    bound = len(a)*int(np.max(np.abs(a), initial=0))*int(np.max(np.abs(b), initial=0))
    assert bound < 2**55
    return a @ b


def transform(poly, counterterm):
    out = dict(poly)
    for mark in range(4):
        out = pmul(out, {0: Fr(1), 1 << mark: -counterterm})
    return out


def finite_case(n, pairing, fields, kinds, assignment, operators):
    eye, cs, number = operators
    length, dimension = len(assignment), len(eye)
    pi = dict(pairing) | {b: a for a, b in pairing}
    blocks = tuple((a, b, k) for (a, b), k in zip(fields, kinds))
    words = [imul(cs[assignment[slot]].T, cs[assignment[pi[slot]]]) for slot in range(length)]
    marked = []
    for (a, b), kind in zip(fields, kinds):
        value = eye
        for word in words[a:b]:
            value = imul(value, word)
        forward = sum(imul(c.T, imul(value, c)) for c in cs)
        if kind == "D":
            twice_loss = 2*(imul(number, value)-forward)
        else:
            reverse = sum(imul(c, imul(value, c.T)) for c in cs)
            twice_loss = n*value-forward-reverse
        marked.append((value, 2*(b-a)*value-twice_loss))
    direct = {}
    maximum = 0
    for mask in range(16):
        value = eye
        for index, variants in enumerate(marked):
            value = imul(value, variants[(mask >> index) & 1])
        maximum = max(maximum, int(np.max(np.abs(value), initial=0)))
        coefficient = Fr(int(np.trace(value)), dimension*(2**mask.bit_count()))
        if coefficient:
            direct[mask] = coefficient
    sequence = [assignment[x] for slot in range(length) for x in (slot, pi[slot])]
    predicted = {0: Fr(sign(sequence))}
    for site in range(n):
        occupied = tuple(slot for slot, index in enumerate(assignment) if index == site)
        factor = local_source(length, pi, blocks, occupied)
        assert factor == sealed_local(length, pi, blocks, occupied)
        predicted = pmul(predicted, factor)
    for index, (a, b) in enumerate(fields):
        predicted = pmul(predicted, {0: Fr(1), 1 << index: Fr(b-a)})
    assert direct == predicted, (n, pairing, kinds, assignment, direct, predicted)
    return direct, predicted, maximum


def finite_checks():
    rows, examples = [], []
    max_integer = 0
    for n, p in PARAMETERS["systems"]:
        assert p == 1
        operators = car_system(n)
        for mode in ("double", "single"):
            fields = tuple(map(tuple, PARAMETERS[mode+"_slot_fields"]))
            length = fields[-1][1]
            kinds_list = PARAMETERS["double_slot_kinds"] if mode == "double" else ["LLLL"]
            assignments = PARAMETERS["double_slot_assignments"][str(n)] if mode == "double" else list(product(range(n), repeat=length))
            counts = Counter()
            for pairing in pairings(tuple(range(length))):
                if not legal(pairing, fields):
                    continue
                counts["pairings"] += 1
                for kinds in kinds_list:
                    for assignment in assignments:
                        actual, prediction, maximum = finite_case(n, pairing, fields, kinds, assignment, operators)
                        max_integer = max(max_integer, maximum)
                        variants = [("raw", actual, prediction)]
                        if mode == "single":
                            # Independent subset expansion for the direct CAR side.
                            shifted = {}
                            lam = Fr(p*p, n)
                            for mask in range(16):
                                value = sum(actual.get(sub, Fr(0))*(-lam)**(mask.bit_count()-sub.bit_count())
                                            for sub in range(16) if sub & mask == sub)
                                if value:
                                    shifted[mask] = value
                            variants.append(("W1", shifted, transform(prediction, lam)))
                        for label, left, right in variants:
                            assert left == right
                            counts[label+"_cases"] += 1
                            counts[label+"_coefficients"] += 16
                            counts[label+"_zero"] += sum(left.get(mask, 0) == 0 for mask in range(16))
                            counts[label+"_negative"] += sum(left.get(mask, 0) < 0 for mask in range(16))
                            counts[label+"_four_mark_nonzero"] += left.get(15, 0) != 0
                            if left.get(15, 0) and not any(e["N"] == n and e["mode"] == mode and e["label"] == label for e in examples):
                                examples.append(dict(N=n, mode=mode, label=label, pairing=pairing,
                                                     kinds=kinds, assignment=assignment,
                                                     coefficients={str(k): str(v) for k, v in sorted(left.items())}))
            assert counts["raw_four_mark_nonzero"] > 0
            if mode == "single":
                assert counts["W1_four_mark_nonzero"] > 0
            rows.append(dict(N=n, p=p, mode=mode, assignments=len(assignments), counts=dict(counts)))
    return rows, examples, max_integer


def set_partitions(items):
    if not items:
        yield ()
        return
    first, rest = items[0], items[1:]
    for partition in set_partitions(rest):
        yield ((first,),) + partition
        for i, block in enumerate(partition):
            yield partition[:i] + ((first,) + block,) + partition[i+1:]


def geometry_checks():
    rows, counts = [], Counter()
    lam = Fr(PARAMETERS["formal_lambda_for_partition_check"])
    for mode in ("double", "single", "mixed"):
        fields = tuple(map(tuple, PARAMETERS[mode+"_slot_fields"]))
        length = fields[-1][1]
        kinds_list = (PARAMETERS["double_slot_kinds"] if mode == "double" else
                      ["LDLD"] if mode == "mixed" else
                      ["".join(x) for x in product("DL", repeat=4)])
        for pairing in pairings(tuple(range(length))):
            if not legal(pairing, fields):
                continue
            pi = dict(pairing) | {b: a for a, b in pairing}
            for kinds in kinds_list:
                blocks = tuple((a, b, k) for (a, b), k in zip(fields, kinds))
                single = [local_source(length, pi, blocks, (i,)) for i in range(length)]
                singleton_product = {0: Fr(1)}
                for factor in single:
                    singleton_product = pmul(singleton_product, {a: 2*x for a, x in factor.items()})
                for index, (a, b) in enumerate(fields):
                    singleton_product = pmul(singleton_product, {0: Fr(1), 1 << index: Fr(b-a)})
                assert singleton_product == {0: Fr(1)}
                kappa = Counter({0: Fr(-comb(length, 2)-length//2)})
                for i, j in combinations(range(length), 2):
                    labels = list(range(length))
                    labels[j] = i
                    sequence = [labels[x] for slot in range(length) for x in (slot, pi[slot])]
                    numerator = {a: sign(sequence)*x for a, x in local_source(length, pi, blocks, (i, j)).items()}
                    ratio = pmul(numerator, pinverse(pmul(single[i], single[j])))
                    kappa.update(ratio)
                kappa = {a: x for a, x in kappa.items() if x}
                assert kappa == sealed_exponent(length, pairing, blocks)
                counts["exponent_cases"] += 1
                counts["all_mask_comparisons"] += 16
                counts["higher_connected_zero_coefficients"] += sum(kappa.get(mask, 0) == 0 for mask in range(16) if mask.bit_count() >= 3)
                assert all(kappa.get(mask, 0) == 0 for mask in range(16) if mask.bit_count() >= 3)
                if mode == "mixed":
                    for mark in (0, 2):
                        assert kappa.get(1 << mark, 0) == 1
                        for dmark in (1, 3):
                            assert kappa.get((1 << mark) | (1 << dmark), 0) == 0
                            counts["single_L_double_D_zero_checks"] += 1
                    shared = any(set(chord) == {0, 3} for chord in pairing)
                    assert kappa.get(5, 0) == int(shared)
                    counts["single_L_pair_checks"] += 1
                # No truncation is used when computing the exponential.
                reduced = {a: lam*x for a, x in kappa.items() if a}
                exponential, power = Counter({0: Fr(1)}), {0: Fr(1)}
                for degree in range(1, 5):
                    power = pmul(power, reduced)
                    exponential.update({a: x/factorial(degree) for a, x in power.items()})
                for mask in range(16):
                    items = tuple(i for i in range(4) if mask >> i & 1)
                    partition_value = Fr(0)
                    for partition in set_partitions(items):
                        value = Fr(1)
                        for block in partition:
                            value *= lam*kappa.get(sum(1 << i for i in block), 0)
                        partition_value += value
                    assert partition_value == exponential.get(mask, 0)
                    counts["partition_comparisons"] += 1
                counts["four_mark_moment_nonzero"] += exponential.get(15, 0) != 0
                rows.append(dict(mode=mode, pairing=pairing, kinds=kinds,
                                 kappa={str(a): str(x) for a, x in sorted(kappa.items())},
                                 four_mark_reduced_moment=str(exponential.get(15, 0))))
    assert counts["four_mark_moment_nonzero"] > 0
    return rows, dict(counts)


def main():
    if OUT.exists():
        raise FileExistsError(OUT)
    start = time.perf_counter()
    finite, examples, maximum = finite_checks()
    geometry, counts = geometry_checks()
    inputs = ["code/check_bilocal_four_sources.py", "code/memory_body_source_chords.py",
              "code/memory_size_marked_chords.py", "code/singlet_multimark.py",
              "code/singlet_coefficient.py", "code/complex_syk_dictionary.py"]
    record = dict(status="passed", schema_version=1,
                  created_utc=datetime.now(timezone.utc).isoformat(),
                  parameters=PARAMETERS, parameter_sha256=digest(PARAMETERS),
                  finite=finite, finite_nonzero_examples=examples,
                  geometry=geometry, geometry_counts=counts,
                  maximum_absolute_integer=maximum,
                  environment=dict(python=platform.python_version(), numpy=np.__version__, platform=platform.platform()),
                  input_sha256={name: sha(ROOT/name) for name in inputs},
                  elapsed_seconds=time.perf_counter()-start,
                  scope="Exact finite fixed-word CAR/source checks and fixed-diagram exponent checks. "
                        "All 16 source masks retained. The pair-Poisson limit and uniform moments have a separate proof. "
                        "No actual full-time or exact-pinched four-point limit is claimed.")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("x", encoding="utf-8") as stream:
        json.dump(record, stream, indent=2)
        stream.write("\n")
    print(json.dumps(dict(status="passed", finite=finite, geometry_counts=counts,
                          maximum_absolute_integer=maximum, elapsed_seconds=record["elapsed_seconds"],
                          output_sha256=sha(OUT))))


if __name__ == "__main__":
    main()
