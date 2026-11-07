"""New exact checks for CP commutators and ordinary charged reassembly.

Imports only sparse CAR primitives from the sealed coherence driver.  The
Gaussian polynomial fixture has ONE scalar Gaussian coordinate, to test exact
number-operator algebra; it is not a sample of the full covariance ensemble.
"""
from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path

import check_bilocal_coherence_channel as car

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/susskind-bilocal-fourth-tail/charged-reassembly-checks.json"


def add(a, b, factor=F(1)):
    out = defaultdict(F, a)
    for entry, value in b.items():
        out[entry] += factor * value
    return car.clean(out)


def mul(a, b):
    out = defaultdict(F)
    rows = defaultdict(list)
    for (row, col), value in b.items():
        rows[row].append((col, value))
    for (row, inner), left in a.items():
        for col, right in rows[inner]:
            out[row, col] += left * right
    return car.clean(out)


def left_charge(op, values):
    return car.clean({(row, col): value * values[row.bit_count()] for (row, col), value in op.items()})


def interior(op, n):
    return {entry: value for entry, value in op.items() if 0 < entry[0].bit_count() < n}


def pclean(poly):
    return {degree: op for degree, op in poly.items() if op}


def padd(a, b, factor=F(1)):
    out = dict(a)
    for degree, op in b.items():
        out[degree] = add(out.get(degree, {}), op, factor)
    return pclean(out)


def pmap(poly, fn):
    return pclean({degree: fn(op) for degree, op in poly.items()})


def pcharge(poly, values):
    return pmap(poly, lambda op: left_charge(op, values))


def number(poly):
    out = {}
    for degree, op in poly.items():
        out = padd(out, {degree: car.scale(op, degree)})
        if degree >= 2:
            out = padd(out, {degree - 2: car.scale(op, -degree * (degree - 1))})
    return out


def setup(n):
    annihilation = [[car.one_action(state, i, False) for state in range(1 << n)] for i in range(n)]
    creation = [[car.one_action(state, i, True) for state in range(1 << n)] for i in range(n)]
    charges = list(range(n + 1))
    reverse = [n - m for m in charges]
    tden = [F(1, m) if m else F(0) for m in charges]
    sden = [F(1, n - m) if m < n else F(0) for m in charges]
    raw_t = lambda op: car.push(op, creation)
    raw_s = lambda op: car.push(op, annihilation)
    d = lambda op: add(left_charge(op, charges), raw_t(op), -1)
    db = lambda op: add(left_charge(op, reverse), raw_s(op), -1)
    return dict(n=n, c=annihilation, cd=creation, Q=lambda op: left_charge(op, charges),
                T=raw_t, S=raw_s, t=lambda op: left_charge(raw_t(op), tden),
                s=lambda op: left_charge(raw_s(op), sden), D=d, Db=db,
                L=lambda op: car.scale(add(d(op), db(op)), F(1, 2)))


def check(a, b, label, counts, restrict=None):
    if restrict is not None:
        a, b = restrict(a), restrict(b)
    entries = set(a) | set(b)
    counts["entry_comparisons"] += len(entries)
    if a != b:
        bad = [(key, str(a.get(key, 0)), str(b.get(key, 0)))
               for key in sorted(entries) if a.get(key, 0) != b.get(key, 0)]
        raise AssertionError((label, bad[:6]))
    counts["matrix_identities"] += 1


def pcheck(a, b, label, counts):
    for degree in set(a) | set(b):
        check(a.get(degree, {}), b.get(degree, {}), (label, degree), counts)
    counts["polynomial_identities"] += 1


def commutator_checks(n, counts):
    z = setup(n)
    Q, T, S, t, s, D, Db = [z[key] for key in ("Q", "T", "S", "t", "s", "D", "Db")]
    q_ratio = [F(m, n - m) if m < n else F(0) for m in range(n + 1)]
    r_ratio = [F(n - m, m) if m else F(0) for m in range(n + 1)]
    q_const = [F(2 * m - n, n - m) if m < n else F(0) for m in range(n + 1)]
    r_const = [F(2 * m - n, m) if m else F(0) for m in range(n + 1)]
    f = lambda m: F(m * m + 3 * m - 1)
    fv = [f(m) for m in range(n + 1)]
    backward = [f(m) - f(m - 1) for m in range(n + 1)]
    forward = [f(m) - f(m + 1) for m in range(n + 1)]
    for row in range(1 << n):
        for col in range(1 << n):
            if row.bit_count() != col.bit_count():
                continue
            counts["neutral_matrix_units"] += 1
            x = {(row, col): F(1)}
            restrict = lambda op: interior(op, n)
            check(add(Q(T(x)), T(Q(x)), -1), T(x), (n, row, col, "QT"), counts)
            check(add(Q(S(x)), S(Q(x)), -1), car.scale(S(x), -1), (n, row, col, "QS"), counts)
            check(add(T(S(x)), S(T(x)), -1), left_charge(x, [2*m-n for m in range(n+1)]),
                  (n, row, col, "TS"), counts)
            check(add(D(t(x)), t(D(x)), -1), add(t(x), t(t(x)), -1),
                  (n, row, col, "Dt"), counts, restrict)
            rhs = add(add(car.scale(s(x), -1), left_charge(t(s(x)), q_ratio)), left_charge(x, q_const), -1)
            check(add(D(s(x)), s(D(x)), -1), rhs, (n, row, col, "Ds"), counts, restrict)
            check(add(Db(s(x)), s(Db(x)), -1), add(s(x), s(s(x)), -1),
                  (n, row, col, "Dbs"), counts, restrict)
            rhs = add(add(car.scale(t(x), -1), left_charge(s(t(x)), r_ratio)), left_charge(x, r_const))
            check(add(Db(t(x)), t(Db(x)), -1), rhs, (n, row, col, "Dbt"), counts, restrict)
            check(add(D(left_charge(x, fv)), left_charge(D(x), fv), -1), left_charge(T(x), backward),
                  (n, row, col, "D f(Q)"), counts)
            check(add(Db(left_charge(x, fv)), left_charge(Db(x), fv), -1), left_charge(S(x), forward),
                  (n, row, col, "Db f(Q)"), counts)


def table_matrix(table):
    return {(output[0], state): F(output[1]) for state, output in enumerate(table) if output is not None}


def jacobi_weights(n, p, max_degree):
    weights = [[F(0) for m in range(n + 1)] for j in range(max_degree + 1)]
    for m in range(n + 1):
        cm = car.choose(m, p) * car.choose(n - m + p, p)
        asum = F(0)
        for j in range(1, max_degree + 1):
            v = j - 1
            theta = F(1) if v == 0 else F(0)
            if m >= p:
                theta = F(1)
                for k in range(p):
                    theta *= max(F(0), 1 - F(v * p, m - k))
                    theta *= max(F(0), 1 - F(v * p, n - m + p - k))
            asum += theta
            weights[j][m] = cm * asum
    return weights


def reassembly_checks(n, p, degree_max, counts):
    import itertools
    z = setup(n)
    sets = tuple(itertools.combinations(range(n), p))
    base = {}
    for i, created in enumerate(sets):
        for j, removed in enumerate(sets):
            coefficient = F(2 * (i + j + 1) + (i == j))
            table = [car.action(state, created, removed) for state in range(1 << n)]
            base = add(base, table_matrix(table), coefficient)
    assert all(base.get((col, row), 0) == value for (row, col), value in base.items())
    check(z["D"](base), car.scale(base, p), (n, p, "D H=pH"), counts)
    identity = {(state, state): F(1) for state in range(1 << n)}
    powers = [identity]
    for k in range(degree_max):
        powers.append(mul(powers[-1], base))
    c = [table_matrix(table) for table in z["c"]]
    cd = [table_matrix(table) for table in z["cd"]]
    f_left = [add(mul(a, base), mul(base, a), -1) for a in cd]
    f_right = [add(mul(base, a), mul(a, base), -1) for a in c]
    r_left = [add(mul(a, base), mul(base, a), -1) for a in c]
    r_right = [add(mul(base, a), mul(a, base), -1) for a in cd]
    def cp(op, left, right):
        ans = {}
        for a, b in zip(left, right):
            ans = add(ans, mul(mul(a, op), b))
        return ans
    gamma = [cp(x, f_left, f_right) for x in powers]
    reverse_gamma = [cp(x, r_left, r_right) for x in powers]
    for w in range(1, degree_max + 1):
        body_rhs = car.scale(powers[w], p * w)
        length_rhs = {}
        for a in range(w):
            length_rhs = add(length_rhs, mul(mul(powers[a], z["L"](base)), powers[w - a - 1]))
        for a in range(w - 1):
            for b in range(w - 1 - a):
                cc = w - 2 - a - b
                body_rhs = add(body_rhs, mul(mul(powers[a], gamma[b]), powers[cc]), -1)
                sym = car.scale(add(gamma[b], reverse_gamma[b]), F(1, 2))
                length_rhs = add(length_rhs, mul(mul(powers[a], sym), powers[cc]), -1)
        check(z["D"](powers[w]), body_rhs, (n, p, w, "D power"), counts)
        check(z["L"](powers[w]), length_rhs, (n, p, w, "L power"), counts)
    weights = jacobi_weights(n, p, degree_max)
    q = [{0: identity}, {1: base}]
    for j in range(1, degree_max):
        hq = {k + 1: mul(base, op) for k, op in q[j].items()}
        q.append(padd(hq, pcharge(q[j - 1], weights[j]), -1))
    pure = []
    for j in range(degree_max + 1):
        pure.append({j - 2 * k: car.scale(powers[j], F((-1)**k * math.factorial(j),
                         2**k * math.factorial(k) * math.factorial(j - 2*k)))
                     for k in range(j // 2 + 1)})
        pcheck(number(pure[-1]), pmap(pure[-1], lambda op: car.scale(op, j)),
               (n, p, j, "pure number eigenvalue"), counts)
    eta = [pmap(padd(a, b, -1), lambda op: car.scale(op, p)) for a, b in zip(pure, q)]
    coeffs = [[F(1) for m in range(n + 1)]]
    fixtures = []
    for degree in range(degree_max + 1):
        if degree:
            new = [[F(0) for m in range(n + 1)] for j in range(degree + 1)]
            for j, row in enumerate(coeffs):
                for m, value in enumerate(row):
                    new[j + 1][m] += value
                    if j:
                        new[j - 1][m] += value * weights[j][m]
            coeffs = new
        reconstructed = {}
        rhs = {}
        difference_term = {}
        for j, hj in enumerate(coeffs):
            reconstructed = padd(reconstructed, pcharge(q[j], hj))
            delta = [hj[m] - (hj[m-1] if m else F(0)) for m in range(n+1)]
            pure_deficit = padd(pmap(pure[j], lambda op: car.scale(op, p*j)), pmap(pure[j], z["D"]), -1)
            rhs = padd(rhs, pcharge(pure_deficit, hj))
            diffz = pcharge(pmap(pure[j], z["T"]), delta)
            difference_term = padd(difference_term, diffz)
            rhs = padd(rhs, diffz, -1)
            rhs = padd(rhs, pcharge(number(eta[j]), hj), -1)
            rhs = padd(rhs, pcharge(pmap(eta[j], z["D"]), hj), F(1, p))
            rhs = padd(rhs, pcharge(pmap(eta[j], z["T"]), delta), F(1, p))
        power_poly = {degree: powers[degree]}
        pcheck(reconstructed, power_poly, (n, p, degree, "inverse Jacobi"), counts)
        lhs = padd(pmap(number(power_poly), lambda op: car.scale(op, p)), pmap(power_poly, z["D"]), -1)
        pcheck(rhs, lhs, (n, p, degree, "five-term charged reassembly"), counts)
        counts["nonzero_charge_difference_fixtures"] += bool(difference_term)
        fixtures.append(dict(N=n, p=p, degree=degree, charge_difference_nonzero=bool(difference_term),
                             reconstruction_residual_nonzero_entries=sum(len(op) for op in padd(reconstructed, power_poly, -1).values()),
                             charged_residual_nonzero_entries=sum(len(op) for op in padd(rhs, lhs, -1).values())))
    return fixtures


def main():
    estimate = {"neutral_basis_units": sum(math.comb(2*n, n) for n in range(2, 6)),
                "superoperator_identities": 9 * sum(math.comb(2*n, n) for n in range(2, 6)),
                "polynomial_fixtures": [[2, 1], [3, 1], [3, 2], [4, 2]],
                "maximum_gaussian_degree": 6, "expected_wall_seconds_upper_bound": 30,
                "algorithm": "sparse Fraction CAR matrices and exact one-coordinate Gaussian polynomials"}
    if "--estimate" in sys.argv:
        print(json.dumps(estimate, indent=2))
        return
    if OUT.exists():
        raise SystemExit(f"Refusing to rerun sealed result: {OUT}")
    start = time.perf_counter()
    counts = dict(neutral_matrix_units=0, matrix_identities=0, entry_comparisons=0,
                  polynomial_identities=0, nonzero_charge_difference_fixtures=0)
    for n in range(2, 6):
        commutator_checks(n, counts)
    fixtures = []
    for n, p in estimate["polynomial_fixtures"]:
        fixtures += reassembly_checks(n, p, 6, counts)
    assert counts["nonzero_charge_difference_fixtures"] > 0
    dependencies = [Path(__file__), Path(car.__file__)]
    result = dict(status="passed", created_utc=datetime.now(timezone.utc).isoformat(),
                  scope="exact finite-CAR identities; no finite-N Schatten inequality or DS fit",
                  parameters=estimate, counts=counts, records=fixtures,
                  gaussian_fixture="H(g)=g A with deterministic Hermitian p-body A and one standard real Gaussian g; number=g*d_g-d_g^2",
                  jacobi_fixture="finite C_m=choose(m,p)choose(N-m+p,p), theta_v=mu_m(vp); algebra valid without an ensemble covariance assumption",
                  checks=["raw CAR superoperators", "normalized CP commutators", "charge-coefficient Leibniz rules",
                          "ordinary D/L power product rules", "pure Hermite number eigenvalues",
                          "inverse Jacobi reconstruction", "all five terms of charged ordinary reassembly"],
                  environment=dict(python=sys.version, platform=platform.platform(), dependencies="Python standard library only"),
                  arithmetic="exact fractions.Fraction",
                  source=dict(path=str(Path(__file__).relative_to(ROOT)).replace("\\", "/"),
                              sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),
                  source_hashes={str(path.relative_to(ROOT)).replace("\\", "/"): hashlib.sha256(path.read_bytes()).hexdigest()
                                 for path in dependencies}, wall_seconds=time.perf_counter() - start)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps({"output": str(OUT), "counts": counts, "wall_seconds": result["wall_seconds"],
                      "source_hashes": result["source_hashes"],
                      "output_sha256": hashlib.sha256(OUT.read_bytes()).hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
