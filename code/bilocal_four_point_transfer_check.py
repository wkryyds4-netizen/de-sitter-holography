"""Exact finite-CAR identities for the four-point transfer interface.

All matrix entries are fractions. This is an algebra check, not a DS fit.
The output is exclusive-create and a completed result must not be rerun.
"""
from __future__ import annotations

import hashlib
import json
import math
import platform
from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/susskind-bilocal-four-point/transfer-checks.json"


def zero(d):
    return [[F(0) for _ in range(d)] for _ in range(d)]


def diag(values):
    out = zero(len(values))
    for i, x in enumerate(values):
        out[i][i] = F(x)
    return out


def add(a, b, factor=F(1)):
    return [[x + factor * y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def scale(a, c):
    return [[F(c) * x for x in row] for row in a]


def mul(a, b):
    d = len(a)
    out = zero(d)
    for i, row in enumerate(a):
        for k, x in enumerate(row):
            if x:
                for j, y in enumerate(b[k]):
                    if y:
                        out[i][j] += x * y
    return out


def transpose(a):
    return [list(row) for row in zip(*a)]


def comm(a, b):
    return add(mul(a, b), mul(b, a), -1)


def trace(a):
    return sum((a[i][i] for i in range(len(a))), F(0)) / len(a)


def norm2(a):
    return trace(mul(transpose(a), a))


def product(items, d):
    ans = diag([1] * d)
    for a in items:
        ans = mul(ans, a)
    return ans


def annihilator(n, site):
    a = zero(1 << n)
    for state in range(1 << n):
        if state & (1 << site):
            sign = -1 if (state & ((1 << site) - 1)).bit_count() % 2 else 1
            a[state ^ (1 << site)][state] = F(sign)
    return a


def choose(n, k):
    return math.comb(n, k) if 0 <= k <= n else 0


def fixture(n, p):
    import itertools

    d = 1 << n
    identity = diag([1] * d)
    cs = [annihilator(n, i) for i in range(n)]
    ads = [transpose(a) for a in cs]
    charges = [i.bit_count() for i in range(d)]
    qop = diag(charges)
    sets = list(itertools.combinations(range(n), p))
    coupling = [[(i + 1) if i == j else ((i + j) % 3 - 1)
                 for j in range(len(sets))] for i in range(len(sets))]
    h = zero(d)
    for i, iset in enumerate(sets):
        for j, jset in enumerate(sets):
            h = add(h, product([ads[k] for k in iset]
                              + [cs[k] for k in reversed(jset)], d), F(coupling[i][j]))
    checks = 0

    def eq(a, b, name):
        nonlocal checks
        assert a == b, (n, p, name)
        checks += 1

    eq(h, transpose(h), "Hermiticity")
    eq(comm(qop, h), zero(d), "charge preservation")
    h2 = mul(h, h)

    def sandwich(x):
        ans = zero(d)
        for ad, c in zip(ads, cs):
            ans = add(ans, mul(mul(ad, x), c))
        return ans

    def body(x):
        return add(mul(qop, x), sandwich(x), -1)

    force = zero(d)
    for c in cs:
        f = comm(h, c)
        force = add(force, mul(transpose(f), f))
    eq(body(h), scale(h, p), "D H = p H")
    eq(force, add(scale(h2, 2 * p), body(h2), -1), "force body identity")

    # Unit covariance w_N=1; this deterministic charge scalar is the GUE E H^2.
    cq = diag([choose(m, p) * choose(n - m + p, p) for m in charges])
    z2 = add(h2, cq, -1)
    expectation_force = add(scale(cq, 2 * p), body(cq), -1)
    oriented = add(force, expectation_force, -1)
    eq(oriented, add(scale(z2, 2 * p), body(z2), -1), "O = (2p-D) Z2")

    lam, formal_q, target_c = F(p * p, n), F(1, 2), F(1)
    kappa = 2 * lam * (1 - formal_q) / (1 + formal_q)
    old_y = add(oriented, z2, -kappa)
    j2 = add(oriented, add(h2, scale(identity, target_c), -1), -kappa)
    eq(add(j2, old_y, -1), scale(add(scale(identity, target_c), cq, -1), kappa),
       "J2 minus old Y is charge scalar")

    # Independent adjacent-charge compression formula for m R^dagger R.
    alpha_energy = sandwich(h)
    inv_q = diag([F(1, m) if m else F(0) for m in charges])
    mg = add(sandwich(h2), mul(inv_q, mul(alpha_energy, alpha_energy)), -1)
    eq(mg, add(force, mul(scale(inv_q, p * p), h2), -1), "mG force formula")
    eq(comm(add(j2, mg, -1), h), zero(d), "J2 and mG differ by conserved field")

    reverse = zero(d)
    for ad, c in zip(ads, cs):
        reverse = add(reverse, mul(mul(c, h), ad))
    half_size_h = scale(add(add(scale(h, n), sandwich(h), -1), reverse, -1), F(1, 2))
    w1 = add(scale(h, p - lam), half_size_h, -1)
    eq(w1, transpose(w1), "local W1 Hermiticity")
    jw = comm(j2, w1)
    aabb = trace(mul(mul(j2, j2), mul(w1, w1)))
    abab = trace(product([j2, w1, j2, w1], d))
    eq([[norm2(jw)]], [[2 * (aabb - abab)]], "commutator four-point identity")
    eq([[trace(jw)]], [[F(0)]], "commutator trace")
    assert norm2(jw) >= 0
    return {"N": n, "p": p, "dimension": d, "covariance_w": 1,
            "formal_lambda": str(lam), "formal_q": str(formal_q),
            "target_C": str(target_c), "kappa": str(kappa), "coupling": coupling,
            "matrix_identity_count": checks, "scalar_entry_comparisons": (checks - 2) * d * d + 2,
            "commutator_norm_squared": str(norm2(jw)),
            "J2_squared_W1_squared": str(aabb), "J2_W1_J2_W1": str(abab)}


def main():
    if OUT.exists():
        raise FileExistsError(f"Sealed output already exists: {OUT}")
    cases = [fixture(n, p) for n, p in [(2, 1), (3, 2), (4, 2)]]
    assert any(F(x["commutator_norm_squared"]) > 0 for x in cases)
    data = {"schema": 1, "status": "passed", "created_utc": datetime.now(timezone.utc).isoformat(),
            "arithmetic": "fractions.Fraction; real Hermitian integer couplings",
            "scope": "exact finite algebra identities; no pinching eigenvector approximation; no DS fit",
            "python": platform.python_version(),
            "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "cases": cases,
            "total_matrix_identities": sum(x["matrix_identity_count"] for x in cases),
            "total_scalar_entry_comparisons": sum(x["scalar_entry_comparisons"] for x in cases),
            "nonzero_commutators": sum(F(x["commutator_norm_squared"]) > 0 for x in cases)}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("x", encoding="utf-8") as stream:
        json.dump(data, stream, indent=2)
        stream.write("\n")
    print(json.dumps({"status": data["status"], "cases": len(cases),
                      "matrix_identities": data["total_matrix_identities"],
                      "scalar_entry_comparisons": data["total_scalar_entry_comparisons"],
                      "nonzero_commutators": data["nonzero_commutators"]}))


if __name__ == "__main__":
    main()
