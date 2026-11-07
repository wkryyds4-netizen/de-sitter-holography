"""Independent exact source regression and no-bridge centering checks.

No project scientific module is imported. q and lambda in the finite graph
checks are independent formal variables evaluated at rational points; this
is an algebraic check, not a physical DS parameter pair or a numerical limit.
"""
from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction
from itertools import combinations, permutations
from pathlib import Path
import hashlib
import json
import platform
import time


ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results/susskind-memory-wick-descendants/source-regression-checks.json"
PARAMETERS = {
    "max_degree": 7,
    "q_exact": "1/3",
    "lambda_exact": "1/4",
    "b0_max_slots": 8,
    "b0_layouts": [
        {"d1": 2, "v1": 2, "d2": 2, "v2": 2,
         "a": a, "b": 4-a, "orientation1": "prefix", "orientation2": "suffix"}
        for a in range(5)
    ] + [
        {"d1": 3, "v1": 2, "d2": 1, "v2": 1,
         "a": 1, "b": 3, "orientation1": side, "orientation2": "prefix"}
        for side in ("prefix", "suffix")
    ],
    "rational_point_scope": "independent formal q and lambda, not q=exp(-4 lambda)",
}


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def multiply(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def qfactorial(degree):
    result = [1]
    for n in range(1, degree+1):
        result = multiply(result, [1]*n)
    return result


def evaluate(poly, q):
    return sum((Fraction(c)*q**n for n, c in enumerate(poly)), Fraction())


def cross(first, second):
    a, b = sorted(first)
    c, d = sorted(second)
    return (a < c < b < d) or (c < a < d < b)


def crossings(chords):
    return sum(cross(a, b) for a, b in combinations(chords, 2))


def kmark(chords):
    return len(chords)**2 - 4*crossings(chords)


def all_pairings(slots):
    if not slots:
        yield ()
        return
    a = slots[0]
    for i in range(1, len(slots)):
        for tail in all_pairings(slots[1:i] + slots[i+1:]):
            yield ((a, slots[i]),) + tail


def permutation_checks():
    rows, counts = [], Counter()
    for d in range(1, PARAMETERS["max_degree"]+1):
        p_d = qfactorial(d)
        weighted = {(v, side): [0]*len(p_d)
                    for v in range(1, d+1) for side in ("prefix", "suffix")}
        negative = Counter()
        minima, maxima = {}, {}
        observed = [0]*len(p_d)
        for sigma in permutations(range(d)):
            counts["permutations"] += 1
            inv = sum(sigma[i] > sigma[j] for i in range(d) for j in range(i+1, d))
            chords = tuple((i, 2*d-1-sigma[i]) for i in range(d))
            assert crossings(chords) == inv
            observed[inv] += 1
            for v, side in weighted:
                chosen = chords[:v] if side == "prefix" else chords[d-v:]
                kval = kmark(chosen)
                weighted[v, side][inv] += kval
                key = (v, side)
                minima[key] = min(minima.get(key, kval), kval)
                maxima[key] = max(maxima.get(key, kval), kval)
                negative[key] += kval < 0
                counts["weighted_mark_samples"] += 1
        assert observed == p_d
        for (v, side), lhs in weighted.items():
            p_v = qfactorial(v)
            d_v = [(v*v-4*n)*c for n, c in enumerate(p_v)]
            cleared_left = multiply(lhs, p_v)
            cleared_right = multiply(p_d, d_v)
            assert cleared_left == cleared_right, (d, v, side)
            counts["regression_cases"] += 1
            counts["negative_raw_mark_samples"] += negative[v, side]
            rows.append({"d": d, "v": v, "orientation": side,
                         "weighted_q_coefficients": lhs,
                         "cleared_polynomial": cleared_left,
                         "negative_mark_samples": negative[v, side],
                         "minimum_k": minima[v, side], "maximum_k": maxima[v, side],
                         "exact_identity": True})
    return rows, counts


def b0_checks():
    q, lam = Fraction(PARAMETERS["q_exact"]), Fraction(PARAMETERS["lambda_exact"])
    rows, counts = [], Counter()
    for layout in PARAMETERS["b0_layouts"]:
        d1, d2, v1, v2, a, b = (layout[k] for k in ("d1", "d2", "v1", "v2", "a", "b"))
        length = d1+a+d2+b
        assert length <= PARAMETERS["b0_max_slots"]
        block1 = set(range(d1))
        block2 = set(range(d1+a, d1+a+d2))
        fields = block1 | block2
        marks = []
        for block, v, side in ((block1, v1, layout["orientation1"]),
                               (block2, v2, layout["orientation2"])):
            ordered = sorted(block)
            marks.append(set(ordered[:v] if side == "prefix" else ordered[-v:]))
        means = []
        for v in (v1, v2):
            p_v = qfactorial(v)
            means.append(lam*evaluate([(v*v-4*n)*c for n, c in enumerate(p_v)], q)/evaluate(p_v, q))
        groups = defaultdict(lambda: {"count": 0, "zero": Fraction(), "one1": Fraction(),
                                      "one2": Fraction(), "two": Fraction(), "centered": Fraction()})
        negative_examples = []
        for pairing in all_pairings(tuple(range(length))):
            counts["b0_candidate_pairings"] += 1
            # No source-internal pair, and no bridge between the sources.
            if any(x in fields and y in fields for x, y in pairing):
                continue
            incident = [tuple(ch for ch in pairing if any(x in mark for x in ch)) for mark in marks]
            assert not set(incident[0]).intersection(incident[1])
            # Fix each source's external partner set; endpoint permutations are
            # summed within this group. Remaining background matching is explicit.
            partners = []
            for block in (block1, block2):
                partners.append(tuple(sorted(y if x in block else x
                                             for x, y in pairing if x in block or y in block)))
            background = tuple(ch for ch in pairing if not any(x in fields for x in ch))
            key = (partners[0], partners[1], background)
            kvals = [kmark(chords) for chords in incident]
            weight = q**crossings(pairing)
            centered = weight*(lam*kvals[0]-means[0])*(lam*kvals[1]-means[1])
            group = groups[key]
            group["count"] += 1
            group["zero"] += weight
            group["one1"] += weight*lam*kvals[0]
            group["one2"] += weight*lam*kvals[1]
            group["two"] += weight*lam*lam*kvals[0]*kvals[1]
            group["centered"] += centered
            counts["b0_pairings"] += 1
            counts["b0_negative_centered_terms"] += centered < 0
            counts["b0_nonzero_centered_terms"] += centered != 0
            if centered < 0 and len(negative_examples) < 2:
                negative_examples.append({"pairing": pairing, "k1": kvals[0], "k2": kvals[1],
                                          "q_crossings": crossings(pairing), "centered_weight": str(centered)})
        assert groups
        group_rows = []
        for key, group in sorted(groups.items()):
            assert group["one1"] == means[0]*group["zero"]
            assert group["one2"] == means[1]*group["zero"]
            assert group["two"] == means[0]*means[1]*group["zero"]
            assert group["centered"] == 0
            group_rows.append({"external_partner_sets": key[:2], "background_pairing": key[2],
                               **{k: str(v) if isinstance(v, Fraction) else v for k, v in group.items()}})
            counts["b0_groups"] += 1
        counts["b0_layouts"] += 1
        rows.append({"layout": layout, "slots": length, "kappa1": str(means[0]), "kappa2": str(means[1]),
                     "groups": group_rows, "negative_centered_examples": negative_examples,
                     "exact_factorization": True, "exact_group_centering": True})
    assert counts["b0_negative_centered_terms"] > 0
    return rows, counts


def main():
    if OUT.exists():
        raise FileExistsError(f"Refusing to overwrite scientific output: {OUT}")
    started = time.perf_counter()
    permutation_rows, counts = permutation_checks()
    b0_rows, b0_counts = b0_checks()
    counts.update(b0_counts)
    assert counts["negative_raw_mark_samples"] > 0
    payload = {
        "schema_version": 1, "status": "passed",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "environment": {"python": platform.python_version(), "platform": platform.platform(),
                        "dependencies": "Python standard library only"},
        "scope": "Independent exact finite source regression and actual no-bridge grouping; no asymptotic proof or physical rational DS point claimed.",
        "parameters": PARAMETERS, "parameter_sha256": digest(PARAMETERS),
        "input_sha256": {"code/check_memory_source_regression.py": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        "counts": dict(counts), "permutation_rows": permutation_rows, "b0_rows": b0_rows,
        "elapsed_seconds": time.perf_counter()-started,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps({"status": payload["status"], "counts": payload["counts"],
                      "parameter_sha256": payload["parameter_sha256"],
                      "output_sha256": hashlib.sha256(OUT.read_bytes()).hexdigest(),
                      "elapsed_seconds": payload["elapsed_seconds"]}, sort_keys=True))


if __name__ == "__main__":
    main()
