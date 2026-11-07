"""Exact finite-CAR covariance/coherence/charge-path checks.

Sparse Fraction arithmetic, w_N=1.  This checks algebraic identities, not
large-N estimates.  A successful artifact is exclusive-created once.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import math
import platform
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/susskind-bilocal-fourth-tail/coherence-channel-checks.json"


def choose(n, k):
    return math.comb(n, k) if 0 <= k <= n else 0


def clean(op):
    return {entry: value for entry, value in op.items() if value}


def scale(op, factor):
    return clean({entry: F(value) * factor for entry, value in op.items()})


def one_action(state, site, create):
    if bool(state & (1 << site)) == create:
        return None
    sign = -1 if (state & ((1 << site) - 1)).bit_count() % 2 else 1
    return state ^ (1 << site), sign


def action(state, created, removed):
    sign = 1
    # C_I^dagger C_J: annihilate J in increasing order, create I reversed.
    for site, create in [(j, False) for j in removed] + [(i, True) for i in reversed(created)]:
        result = one_action(state, site, create)
        if result is None:
            return None
        state, step_sign = result
        sign *= step_sign
    return state, sign


def push(op, tables):
    out = defaultdict(F)
    for table in tables:
        for (row, col), value in op.items():
            left, right = table[row], table[col]
            if left is not None and right is not None:
                out[left[0], right[0]] += value * left[1] * right[1]
    return clean(out)


def pair_sets(n, a):
    if a == 0:
        return (), ()
    if a == 1:
        return ((1,), (3,)) if n >= 4 else ((0,), (n - 1,))
    return ((0, 2), (3, 4)) if n >= 5 else ((0, 2), (1, 3))


def scalar_tests(rest):
    def occ(state, j):
        return int(bool(state & (1 << rest[j])))
    tests = [("constant", lambda state: 1)]
    if rest:
        tests.append(("one_site", lambda state: occ(state, 0)))
    if len(rest) >= 2:
        tests += [
            ("mixed_projector", lambda state: occ(state, 0) * (1 - occ(state, 1))),
            ("signed_polynomial", lambda state: 2 * occ(state, 0) - occ(state, 1)
             + 3 * occ(state, 0) * occ(state, 1) - 1),
        ]
    return tests


def make_coherence(n, a, m, created, removed, function):
    rest = tuple(i for i in range(n) if i not in created + removed)
    common = [sum(1 << i for i in sites) for sites in itertools.combinations(rest, m - a)]
    op, gauges = {}, {}
    removed_mask = sum(1 << i for i in removed)
    for state in common:
        col = state | removed_mask
        row, sign = action(col, created, removed)
        gauges[state] = sign
        value = F(sign * function(state))
        if value:
            op[row, col] = value
    return op, common, gauges


def coherence_prediction(n, p, a, m, created, removed, function, common, gauges):
    out = {}
    if m - a < p:
        return out
    created_mask = sum(1 << i for i in created)
    removed_mask = sum(1 << i for i in removed)
    for target in common:
        value = sum((F(choose(m - a - (source & ~target).bit_count(),
                              p - (source & ~target).bit_count())) * function(source)
                     for source in common), F(0))
        if value:
            out[target | created_mask, target | removed_mask] = gauges[target] * value
    return out


def charge_path(op, n, p, m, annihilation, creation):
    ans = op
    for target in range(m - 1, m - p - 1, -1):
        ans = scale(push(ans, annihilation), F(1, n - target))
    for target in range(m - p + 1, m + 1):
        ans = scale(push(ans, creation), F(1, target))
    return ans


def check_equal(actual, expected, context, counts):
    entries = set(actual) | set(expected)
    counts["sparse_entries_compared"] += len(entries)
    if actual != expected:
        bad = [(entry, str(actual.get(entry, 0)), str(expected.get(entry, 0)))
               for entry in sorted(entries) if actual.get(entry, 0) != expected.get(entry, 0)]
        raise AssertionError((context, bad[:5]))


def estimate():
    cases, direct_visits = 0, 0
    for n in range(2, 7):
        for p in (1, 2):
            for a in range(min(2, n // 2) + 1):
                for m in range(a, n - a + 1):
                    tests = 1 if n == 2 * a else 2 if n == 2 * a + 1 else 4
                    cases += tests
                    direct_visits += tests * choose(n, p) ** 2 * choose(n - 2 * a, m - a)
    return {"coherence_fixtures": cases, "direct_kraus_entry_visits_upper_bound": direct_visits,
            "algorithm": "sparse exact Fractions; cached CAR action tables; no dense products",
            "expected_wall_seconds_upper_bound": 30}


def main():
    if "--estimate" in sys.argv:
        print(json.dumps(estimate(), indent=2))
        return
    if OUT.exists():
        raise SystemExit(f"Refusing to rerun sealed result: {OUT}")
    start = time.perf_counter()
    counts = dict(coherence_identities=0, charge_path_identities=0, inactive_zero_cases=0,
                  kernel_rows=0, kernel_entries=0, tagged_flip_bounds=0,
                  sparse_entries_compared=0, nonconstant_input_fixtures=0,
                  varying_gauge_sectors=0, nonzero_outputs=0)
    fixtures = []
    gauge_witness = None
    for n in range(2, 7):
        annihilation = [[one_action(state, i, False) for state in range(1 << n)] for i in range(n)]
        creation = [[one_action(state, i, True) for state in range(1 << n)] for i in range(n)]
        for p in (1, 2):
            sets = tuple(itertools.combinations(range(n), p))
            kraus = [[action(state, i, j) for state in range(1 << n)] for i in sets for j in sets]
            for a in range(min(2, n // 2) + 1):
                created, removed = pair_sets(n, a)
                rest = tuple(i for i in range(n) if i not in created + removed)
                for m in range(a, n - a + 1):
                    cm = choose(m, p) * choose(n - m + p, p)
                    kappa = choose(m - a, p) * choose(n - m - a + p, p)
                    _, common, gauges = make_coherence(n, a, m, created, removed, lambda state: 1)
                    if len(set(gauges.values())) > 1:
                        counts["varying_gauge_sectors"] += 1
                        if gauge_witness is None:
                            gauge_witness = dict(N=n, p=p, a=a, m=m, A=created, B=removed,
                                                 common_gauges={str(k): v for k, v in gauges.items()})
                    if kappa:
                        for source in common:
                            probs = {target: F(choose(m - a - (source & ~target).bit_count(),
                                                     p - (source & ~target).bit_count()), kappa)
                                     for target in common}
                            assert sum(probs.values(), F(0)) == 1
                            assert all(value >= 0 for value in probs.values())
                            counts["kernel_rows"] += 1
                            counts["kernel_entries"] += len(probs)
                            for site in rest:
                                flip = sum((prob for target, prob in probs.items()
                                            if bool(source & (1 << site)) != bool(target & (1 << site))), F(0))
                                bound = F(p, m - a) + F(p, n - m - a + p)
                                assert flip <= bound
                                counts["tagged_flip_bounds"] += 1
                    for name, function in scalar_tests(rest):
                        op, _, _ = make_coherence(n, a, m, created, removed, function)
                        direct = push(op, kraus)
                        predicted = coherence_prediction(n, p, a, m, created, removed, function, common, gauges)
                        context = dict(N=n, p=p, a=a, m=m, function=name)
                        check_equal(direct, predicted, ("coherence", context), counts)
                        counts["coherence_identities"] += 1
                        counts["nonzero_outputs"] += bool(direct)
                        counts["nonconstant_input_fixtures"] += len({function(state) for state in common}) > 1
                        if cm:
                            path = charge_path(op, n, p, m, annihilation, creation)
                            check_equal(path, scale(direct, F(1, cm)), ("CP path", context), counts)
                            counts["charge_path_identities"] += 1
                        else:
                            assert not direct
                            counts["inactive_zero_cases"] += 1
                        fixtures.append({**context, "C_m": cm, "kappa_a": kappa,
                                         "input_nonzero_entries": len(op), "output_nonzero_entries": len(direct)})
    assert gauge_witness is not None
    assert counts["nonconstant_input_fixtures"] > 0
    assert counts["coherence_identities"] == estimate()["coherence_fixtures"]
    result = {
        "status": "passed", "scope": "exact finite-CAR algebra only; no asymptotic or Schatten estimate",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "parameters": {"N": [2, 3, 4, 5, 6], "p": [1, 2], "a": [0, 1, 2],
                       "m": "every a <= m <= N-a; a <= floor(N/2)", "covariance_weight_w_N": 1,
                       "arithmetic": "fractions.Fraction", "CAR_site_numbering": "zero-based"},
        "identities": ["Phi_m(T_AB f)=kappa_a T_AB P_a f (actual CAR gauge)",
                       "Phi_m/C_m=T_m...T_(m-p+1) S_(m-p)...S_(m-1) for m>=p",
                       "positive symmetric Markov kernel, unit row sums and tagged-site flip bound"],
        "counts": counts, "gauge_witness": gauge_witness, "fixtures": fixtures,
        "estimate": estimate(), "wall_seconds": time.perf_counter() - start,
        "environment": {"python": sys.version, "platform": platform.platform(), "dependencies": "Python standard library only"},
        "source": {"path": str(Path(__file__).relative_to(ROOT)).replace("\\", "/"),
                   "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(result, stream, indent=2, ensure_ascii=False)
        stream.write("\n")
    print(json.dumps({"output": str(OUT), "status": result["status"], "counts": counts,
                      "wall_seconds": result["wall_seconds"], "source_sha256": result["source"]["sha256"],
                      "output_sha256": hashlib.sha256(OUT.read_bytes()).hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
