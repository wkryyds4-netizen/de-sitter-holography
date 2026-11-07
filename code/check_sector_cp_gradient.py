"""Finite CAR diagnostics for a candidate fixed-charge gradient argument.

Integer identities are exact (small int64 fixtures, with a conservative bound).
PSD and Schatten checks use floating-point eigenvalues/SVD and are diagnostics.
This does not test or prove the marked-Wick estimate or an ordinary-series tail.
The output is exclusive: never overwrite a successful scientific run.
"""
from pathlib import Path
import hashlib
import json
import math
import platform
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/sector-cp-gradient/finite-car-checks.json"


def car_blocks(n):
    sectors = [[a for a in range(1 << n) if a.bit_count() == k]
               for k in range(n + 1)]
    blocks = {}
    for k in range(1, n + 1):
        rows = {a: i for i, a in enumerate(sectors[k - 1])}
        blocks[k] = []
        for j in range(n):
            c = np.zeros((len(rows), len(sectors[k])), dtype=np.int64)
            for col, a in enumerate(sectors[k]):
                if a & (1 << j):
                    c[rows[a ^ (1 << j)], col] = (-1) ** ((a & ((1 << j) - 1)).bit_count())
            blocks[k].append(c)
    return [len(s) for s in sectors], blocks


def sandwich(c, x, reverse=False):
    return sum((a @ x @ a.T if reverse else a.T @ x @ a) for a in c)


def gram(y, row=False):
    return sum(a @ a.T if row else a.T @ a for a in y)


def norm(x, order):
    v = np.linalg.svd(x.astype(float), compute_uv=False)
    return float(max(v) if order == "inf" else np.mean(v ** order) ** (1 / order))


def main():
    if OUT.exists():
        raise SystemExit("Saved result exists; do not rerun or overwrite.")
    seed = 20261108
    rng = np.random.default_rng(seed)
    exact_count = 0
    psd_rows = []
    norm_rows = []
    max_integer = 0

    def equal(a, b, label):
        nonlocal exact_count, max_integer
        max_integer = max(max_integer, int(np.max(np.abs(a))), int(np.max(np.abs(b))))
        assert np.array_equal(a, b), label
        exact_count += 1

    def psd(lhs, rhs, label):
        nonlocal max_integer
        max_integer = max(max_integer, int(np.max(np.abs(lhs))), int(np.max(np.abs(rhs))))
        gap = rhs - lhs
        assert np.array_equal(gap, gap.T)
        scale = max(1.0, float(np.linalg.norm(rhs.astype(float), 2)),
                    float(np.linalg.norm(lhs.astype(float), 2)))
        value = float(np.linalg.eigvalsh(gap.astype(float))[0])
        tol = 2e-12 * scale
        assert value >= -tol, (label, value, tol)
        psd_rows.append({"case": label, "min_eigenvalue": value,
                         "scale": scale, "tolerance": tol})

    for n in range(3, 7):
        dims, c = car_blocks(n)
        lcm = math.lcm(*range(1, n + 1))
        # Integer entries <=3, dimension <=20, n<=6: all intermediate
        # products here are bounded conservatively by 10**12 < 2**63.
        assert max(dims) <= 20 and lcm <= 60
        for trial in range(3):
            base = [rng.integers(-3, 4, (d, d), dtype=np.int64) for d in dims]
            for adjoint in (False, True):
                x = [a.T.copy() if adjoint else a for a in base]
                tx = [np.zeros_like(x[0])] + [
                    (lcm // k) * sandwich(c[k], x[k - 1]) for k in range(1, n + 1)]
                sx = [(lcm // (n - k)) * sandwich(c[k + 1], x[k + 1], True)
                      for k in range(n)] + [np.zeros_like(x[n])]
                y = {k: [a @ x[k] - x[k - 1] @ a for a in c[k]]
                     for k in range(1, n + 1)}
                for m in range(1, n + 1):
                    d = n - m + 1
                    for i, a in enumerate(c[m]):
                        label = f"N{n}/trial{trial}/adj{int(adjoint)}/m{m}/i{i}"
                        pt = (sum(c[m - 1][j].T @ y[m - 1][i] @ c[m][j]
                                  for j in range(n)) if m > 1 else np.zeros_like(a))
                        ps = (sum(c[m][j] @ y[m + 1][i] @ c[m + 1][j].T
                                  for j in range(n)) if m < n else np.zeros_like(a))
                        equal(m * (a @ tx[m] - tx[m - 1] @ a),
                              -lcm * pt + (lcm * x[m - 1] - tx[m - 1]) @ a,
                              "T/" + label)
                        equal(d * (a @ sx[m] - sx[m - 1] @ a),
                              -lcm * ps + a @ (sx[m] - lcm * x[m]),
                              "S/" + label)

            # Arbitrary rectangular arrays, independent of gradient structure.
            for m in range(1, n + 1):
                if m > 1:
                    y = [rng.integers(-3, 4, (dims[m - 2], dims[m - 1]), dtype=np.int64)
                         for _ in range(n)]
                    p = [sum(c[m - 1][j].T @ a @ c[m][j] for j in range(n)) for a in y]
                    psd(gram(p), (m - 1) * sandwich(c[m], gram(y)), f"T/col/{n}/{trial}/{m}")
                    psd(gram(p, True), m * sandwich(c[m - 1], gram(y, True)), f"T/row/{n}/{trial}/{m}")
                if m < n:
                    d = n - m + 1
                    y = [rng.integers(-3, 4, (dims[m], dims[m + 1]), dtype=np.int64)
                         for _ in range(n)]
                    p = [sum(c[m][j] @ a @ c[m + 1][j].T for j in range(n)) for a in y]
                    psd(gram(p), d * sandwich(c[m + 1], gram(y), True), f"S/col/{n}/{trial}/{m}")
                    psd(gram(p, True), (d - 1) * sandwich(c[m], gram(y, True), True), f"S/row/{n}/{trial}/{m}")

            for k in range(n + 1):
                for name in ("T", "S"):
                    if (name == "T" and k == 0) or (name == "S" and k == n):
                        continue
                    src = k - 1 if name == "T" else k + 1
                    denom = k if name == "T" else n - k
                    mats = c[k] if name == "T" else c[k + 1]
                    reverse = name == "S"
                    z = sandwich(mats, base[src], reverse)
                    equal(sandwich(mats, np.eye(dims[src], dtype=np.int64), reverse),
                          denom * np.eye(dims[k], dtype=np.int64), f"unital/{n}/{k}/{name}")
                    equal(np.array([np.trace(z) * dims[src]]),
                          np.array([denom * dims[k] * np.trace(base[src])]), f"trace/{n}/{k}/{name}")
                    for order in (2, 4, "inf"):
                        before = norm(base[src], order)
                        after = norm(z, order) / denom
                        assert after <= before + 1e-11 * max(1.0, before)
                        norm_rows.append({"N": n, "trial": trial, "k": k, "map": name,
                                          "order": order, "source": before, "target": after})

    assert max_integer < 10**12
    result = {
        "status": "passed", "scope": "Finite CAR identities and CP principal-part diagnostics only; no marked-Wick or ordinary-tail theorem.",
        "parameters": {"N": [3, 4, 5, 6], "trials": 3, "seed": seed,
                       "entries": "integer -3..3; non-selfadjoint neutral blocks and their transposes",
                       "arbitrary_rectangular_arrays": True},
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
        "exact_integer_equalities": exact_count, "max_observed_integer": max_integer,
        "conservative_intermediate_integer_bound": 10**12,
        "psd_checks": psd_rows, "schatten_checks": norm_rows,
        "code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "limitations": "Finite fixtures do not prove the identities generally. PSD/SVD tests use numerical tolerances; no interval certification. No degree-uniform Wick bound, CP path test or infinite-series result is asserted."
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("x", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
        f.write("\n")
    print(json.dumps({"status": result["status"], "exact": exact_count,
                      "psd": len(psd_rows), "schatten": len(norm_rows),
                      "output_sha256": hashlib.sha256(OUT.read_bytes()).hexdigest()}))


if __name__ == "__main__":
    main()
