"""Independent finite-CAR and exact geometric checks for new oriented marks.

Imports sealed helper functions, never their scientific main routines.
The raw nonadjacent shortcut failure is an expected retained diagnostic.
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

from complex_syk_dictionary import fermions, creators
from singlet_multimark import local as old_pair_local
from memory_oriented_marked_chords import (
    I2, C2, CD2, KINDS, mm, loss, deficit, exponent, geometric_exponent,
    finite_local, pair_site_terms,
)

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results/susskind-memory-oriented-transfer/oriented-mark-checks.json"
PARAMETERS = {
    "seed": 110711,
    "local_word_max_length": 6,
    "geometry_layouts": [
        [2, [[0, 1], [1, 2]]], [4, [[0, 2], [2, 4]]],
        [6, [[0, 2], [3, 5]]], [6, [[0, 3], [3, 6]]],
        [8, [[0, 3], [4, 7]]], [8, [[0, 4], [4, 8]]],
    ],
    "finite_systems": [
        [3, 1, 2, [[0, 1], [1, 2]]],
        [3, 1, 4, [[0, 2], [2, 4]]],
        [4, 2, 4, [[0, 2], [2, 4]]],
        [4, 1, 6, [[0, 3], [3, 6]]],
        [4, 2, 6, [[0, 2], [3, 5]]],
    ],
    "finite_words_per_system": 16,
    "kinds": list(KINDS),
    "loss_forms": ["asymmetric", "symmetric"],
    "conditional_two_slot_lengths": [4, 6, 8],
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def pairings(slots):
    if not slots:
        yield ()
    else:
        for i in range(1, len(slots)):
            for tail in pairings(slots[1:i]+slots[i+1:]):
                yield ((slots[0], slots[i]),)+tail


def allowed(pairing, blocks):
    return not any(a <= x < b and a <= y < b
                   for a, b in blocks for x, y in pairing)


def serialize_map(value):
    return {str(key): str(coefficient) for key, coefficient in sorted(value.items()) if coefficient}


def local_checks():
    counts = Counter()
    rows = []
    n = mm(CD2, C2)
    for length in range(PARAMETERS["local_word_max_length"]+1):
        for letters in product(("c", "cd"), repeat=length):
            x = I2
            for letter in letters:
                x = mm(x, C2 if letter == "c" else CD2)
            nc, na = letters.count("cd"), letters.count("c")
            commutator = tuple(u-v for u, v in zip(mm(n, x), mm(x, n)))
            assert commutator == tuple((nc-na)*v for v in x)
            for kind in KINDS:
                asym = deficit(x, nc, na, kind, "asymmetric")
                sym = deficit(x, nc, na, kind, "symmetric")
                assert asym == sym
                counts["local_deficit_comparisons"] += 1
            counts["local_words"] += 1
            counts["local_nonzero_odd_words"] += bool(any(x) and length % 2)
            if any(x):
                rows.append({"word": list(letters), "creation_count": nc, "annihilation_count": na,
                             "matrix": list(map(str, x)),
                             "D_asymmetric": list(map(str, loss(x, "D", "asymmetric"))),
                             "D_symmetric": list(map(str, loss(x, "D", "symmetric"))),
                             "D_deficit_both_forms": list(map(str, deficit(x, nc, na, "D")))})
    assert counts["local_nonzero_odd_words"] > 0
    return rows, counts


def geometry_checks():
    rows, counts = [], Counter()
    for length, raw_blocks in PARAMETERS["geometry_layouts"]:
        blocks = tuple(map(tuple, raw_blocks))
        for pairing in pairings(tuple(range(length))):
            if not allowed(pairing, blocks):
                continue
            counts["geometry_pairings"] += 1
            for kinds in product(KINDS, repeat=2):
                actual = exponent(length, pairing, blocks, kinds)
                sym = exponent(length, pairing, blocks, kinds, "symmetric")
                expected = geometric_exponent(length, pairing, blocks, kinds)
                for mask in range(4):
                    assert actual.get(mask, 0) == expected.get(mask, 0), (length, pairing, kinds, mask)
                    assert sym.get(mask, 0) == actual.get(mask, 0)
                    counts["geometry_coefficient_comparisons"] += 2
                counts["geometry_cases"] += 1
                counts["negative_single_coefficients"] += sum(actual.get(k, 0) < 0 for k in (1, 2))
                counts["negative_connected_coefficients"] += actual.get(3, 0) < 0
                rows.append({"length": length, "blocks": blocks, "pairing": pairing, "kinds": kinds,
                             "kappa": serialize_map(actual), "both_loss_forms_exact": True,
                             "geometry_exact": True})
    assert counts["negative_connected_coefficients"] > 0
    return rows, counts


def pair_site_table():
    # Independent expected classification. Column order matches KINDS products.
    same_expected = {("D", "D"): 0, ("D", "L"): 0, ("D", "Db"): 0,
                     ("L", "D"): 0, ("L", "L"): 1, ("L", "Db"): 2,
                     ("Db", "D"): 0, ("Db", "L"): 2, ("Db", "Db"): 4}
    noncross_expected = {("D", "D"): 4, ("D", "L"): 2, ("D", "Db"): 0,
                        ("L", "D"): 2, ("L", "L"): 2, ("L", "Db"): 2,
                        ("Db", "D"): 0, ("Db", "L"): 2, ("Db", "Db"): 4}
    cross_expected = {k: v-4 for k, v in noncross_expected.items()}
    rows = []
    for category, length, pa, blocks, expected in [
        ("same_shared_chord", 2, ((0, 1),), ((0, 1), (1, 2)), same_expected),
        ("distinct_shared_noncrossing", 4, ((0, 3), (1, 2)), ((0, 2), (2, 4)), noncross_expected),
        ("distinct_shared_crossing", 4, ((0, 2), (1, 3)), ((0, 2), (2, 4)), cross_expected),
    ]:
        for kinds in product(KINDS, repeat=2):
            entries = pair_site_terms(length, pa, blocks, kinds)
            # For distinct chords, discard the two within-chord label overlaps.
            relevant = [(labels, coefficients) for labels, coefficients in entries
                        if category == "same_shared_chord" or labels not in pa]
            connected = sum((v.get(3, F(0)) for _, v in relevant), F(0))
            assert connected == expected[kinds]
            rows.append({"category": category, "kinds": kinds, "connected": str(connected),
                         "signed_pair_site_terms": [
                             {"occupied_labels": labels, "kappa": serialize_map(values)}
                             for labels, values in relevant]})
    return rows


def shortcut_checks():
    counts, rows = Counter(), []
    for length in PARAMETERS["conditional_two_slot_lengths"]:
        for j in range(2, length//2+1):
            ncases = 0
            for pa in pairings(tuple(range(length))):
                if not allowed(pa, ((0, j),)):
                    continue
                for a, b in combinations(range(j), 2):
                    # The old routine asserts its alpha formula internally.
                    old_pair_local(length, pa, ((a, b),))
                    ncases += 1
            rows.append({"length": length, "source_length": j, "legal_two_slot_cases": ncases})
            counts["legal_two_slot_cases"] += ncases
    assert counts["legal_two_slot_cases"] == 446
    raw_pa, mark = ((0, 1), (2, 3)), (0, 2)
    rejected = False
    try:
        old_pair_local(4, raw_pa, (mark,))
    except AssertionError:
        rejected = True
    assert rejected
    # Direct derivative of the two-by-two determinant minor: only entry (2,1)
    # changes by -2; its paired minor has e12=+1 and diagonal sign product -1.
    raw_alpha = -1 * (1 * -2) * -1
    invalid_geometric_shortcut = 2
    assert raw_alpha == -2 and raw_alpha != invalid_geometric_shortcut
    diagnostic = {"status": "expected_failure_retained", "pairing": raw_pa, "slot_mark": mark,
                  "spanning_interval": [0, 3], "violated_condition": "source interval contains an internal pairing",
                  "direct_minor_alpha": raw_alpha, "invalid_shortcut_alpha": invalid_geometric_shortcut,
                  "old_assertion_rejected": rejected}
    return rows, diagnostic, counts


def exact_trace(matrix):
    trace = np.trace(matrix)
    assert trace.imag == 0
    scaled = trace.real*4
    assert scaled == round(scaled)
    return F(int(round(scaled)), 4*matrix.shape[0])


def global_losses(X, cs, Q, kind, form):
    eye = np.eye(X.shape[0])
    left = Q@X if form == "asymmetric" else (Q@X+X@Q)/2
    right = (len(cs)*eye-Q)@X if form == "asymmetric" else ((len(cs)*eye-Q)@X+X@(len(cs)*eye-Q))/2
    normal = left-sum(c.conj().T@X@c for c in cs)
    reverse = right-sum(c@X@c.conj().T for c in cs)
    return normal if kind == "D" else reverse if kind == "Db" else (normal+reverse)/2


def finite_checks():
    rng = np.random.default_rng(PARAMETERS["seed"])
    rows, counts = [], Counter()
    for N, p, length, raw_blocks in PARAMETERS["finite_systems"]:
        blocks = tuple(map(tuple, raw_blocks))
        cs = fermions(N)
        Q = sum(c.conj().T@c for c in cs)
        subsets = list(combinations(range(N), p))
        create = dict(zip(subsets, creators(cs, p)))
        valid = [pa for pa in pairings(tuple(range(length))) if allowed(pa, blocks)]
        eye = np.eye(1 << N)
        for case in range(PARAMETERS["finite_words_per_system"]):
            pa = valid[case % len(valid)]
            pi = dict(pa) | {b: a for a, b in pa}
            sets = ([subsets[case % len(subsets)]]*length if case < 4 else
                    [subsets[int(rng.integers(len(subsets)))] for _ in range(length)])
            ops = [create[sets[i]]@create[sets[pi[i]]].conj().T for i in range(length)]
            matrices, marked = [], []
            for start, end in blocks:
                X = eye.copy()
                for op in ops[start:end]:
                    X = X@op
                assert np.array_equal(Q@X, X@Q)
                matrices.append(X)
                by_kind = {}
                for kind in KINDS:
                    a = global_losses(X, cs, Q, kind, "asymmetric")
                    b = global_losses(X, cs, Q, kind, "symmetric")
                    assert np.array_equal(a, b)
                    by_kind[kind] = (end-start)*p*X-a
                    counts["finite_neutral_loss_comparisons"] += 1
                marked.append(by_kind)
            starts = {start: (i, end) for i, (start, end) in enumerate(blocks)}
            counts["finite_words"] += 1
            for kinds in product(KINDS, repeat=2):
                actual = []
                for mask in range(4):
                    matrix, pos = eye.copy(), 0
                    while pos < length:
                        if pos in starts:
                            i, end = starts[pos]
                            matrix = matrix@(marked[i][kinds[i]] if mask >> i & 1 else matrices[i])
                            pos = end
                        else:
                            matrix, pos = matrix@ops[pos], pos+1
                    actual.append(exact_trace(matrix))
                expected = finite_local(length, sets, pa, blocks, kinds, N)
                symmetric = finite_local(length, sets, pa, blocks, kinds, N, "symmetric")
                assert actual == expected == symmetric, (N, p, pa, sets, kinds, actual, expected)
                counts["finite_type_cases"] += 1
                counts["finite_trace_comparisons"] += 8
                counts["finite_negative_double_marks"] += actual[3] < 0
                counts["finite_zero_base_nonzero_mark"] += actual[0] == 0 and any(actual[1:])
                rows.append({"N": N, "p": p, "length": length, "blocks": blocks,
                             "pairing": pa, "sets": sets, "kinds": kinds,
                             "all_mark_subsets": list(map(str, actual)), "both_loss_forms_exact": True})
    assert counts["finite_zero_base_nonzero_mark"] > 0
    return rows, counts


def calibration_checks():
    rows = []
    for N, p in ((3, 1), (4, 2)):
        cs = fermions(N)
        Q = sum(c.conj().T@c for c in cs)
        subsets = list(combinations(range(N), p))
        create = dict(zip(subsets, creators(cs, p)))
        for A in subsets:
            for B in subsets:
                H = create[A]@create[B].conj().T
                normal = global_losses(H, cs, Q, "D", "asymmetric")
                sym = global_losses(H, cs, Q, "D", "symmetric")
                assert np.array_equal(normal, p*H) and np.array_equal(normal, sym)
                rows.append({"N": N, "p": p, "A": A, "B": B, "p_minus_D_H_exactly_zero": True})
    return rows


def main():
    if OUT.exists():
        raise FileExistsError(f"Refusing to overwrite scientific output: {OUT}")
    started = time.perf_counter()
    local_rows, counts = local_checks()
    geometry_rows, more = geometry_checks()
    counts.update(more)
    table = pair_site_table()
    shortcut_rows, diagnostic, more = shortcut_checks()
    counts.update(more)
    finite_rows, more = finite_checks()
    counts.update(more)
    calibration = calibration_checks()
    counts["pair_site_table_rows"] = len(table)
    counts["p_minus_D_H_calibrations"] = len(calibration)
    names = ["code/memory_oriented_marked_chords.py", "code/check_memory_oriented_marks.py",
             "code/complex_syk_dictionary.py", "code/singlet_multimark.py", "code/singlet_coefficient.py"]
    data = {"schema_version": 1, "status": "passed", "created_utc": datetime.now(timezone.utc).isoformat(),
            "scope": "Neutral even coupling-Wick source blocks, no internal source pairings; asym/sym loss agree with their respective incidence assignments. Not whole-field Gaussian-degree marks with internal contractions; not global odd source blocks.",
            "normalization": "Standard CAR, normalized full trace; DS graph weight C^(length/2) q^crossings with q=exp(-4 lambda), lambda=p^2/N; no DS-limit numerical approximation used here.",
            "local_identity": "(nc-D_asym)x=((nc+na)/2-D_sym)x since [n,x]=(nc-na)x; summed losses coincide on neutral X because [Q,X]=0.",
            "environment": {"python": platform.python_version(), "numpy": np.__version__, "platform": platform.platform()},
            "parameters": PARAMETERS, "parameter_sha256": canonical_hash(PARAMETERS),
            "input_sha256": {name: sha(ROOT/name) for name in names}, "counts": dict(counts),
            "local_word_rows": local_rows, "geometry_rows": geometry_rows, "pair_site_table": table,
            "conditional_two_slot_rows": shortcut_rows, "raw_nonadjacent_shortcut_diagnostic": diagnostic,
            "finite_car_rows": finite_rows, "calibration_rows": calibration,
            "elapsed_seconds": time.perf_counter()-started}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("x", encoding="utf-8") as stream:
        json.dump(data, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps({"status": data["status"], "counts": data["counts"],
                      "parameter_sha256": data["parameter_sha256"], "output_sha256": sha(OUT),
                      "elapsed_seconds": data["elapsed_seconds"]}, sort_keys=True))


if __name__ == "__main__":
    main()
