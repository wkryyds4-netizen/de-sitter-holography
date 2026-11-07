"""Exact finite CAR checks for restricted covariances and tagged Wick networks.

New output only; never imports or reruns historical scientific mains.
The integer fixtures use coupling weight w=1. Both sides scale homogeneously
with the physical positive variance w_N, which is also recorded.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import platform
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/susskind-bilocal-fourth-tail/tagged-covariance-checks.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def choose(n: int, k: int) -> int:
    return math.comb(n, k) if 0 <= k <= n else 0


def subsets(n: int, max_size: int) -> list[tuple[int, ...]]:
    return [s for k in range(max_size + 1) for s in itertools.combinations(range(n), k)]


def annihilator(n: int, i: int) -> np.ndarray:
    a = np.zeros((1 << n, 1 << n), dtype=np.int64)
    for state in range(1 << n):
        if state & (1 << i):
            a[state ^ (1 << i), state] = (-1) ** ((state & ((1 << i) - 1)).bit_count())
    return a


def system(n: int, p: int) -> dict:
    dim = 1 << n
    cs = [annihilator(n, i) for i in range(n)]
    ps = list(itertools.combinations(range(n), p))
    creators = {}
    for indices in ps:
        c = np.eye(dim, dtype=np.int64)
        for i in indices:
            c = c @ cs[i].T
        creators[indices] = c
    labels, mats = [], []
    for left in ps:
        for right in ps:
            labels.append((frozenset(left), frozenset(right)))
            mats.append(creators[left] @ creators[right].T)
    matrices = np.asarray(mats, dtype=np.int64)
    col = np.asarray([a.T @ a for a in matrices])
    row = np.asarray([a @ a.T for a in matrices])
    offdiag = 0
    for matrix in itertools.chain(col, row):
        offdiag = max(offdiag, int(np.max(np.abs(matrix - np.diag(np.diag(matrix))))))
    lam = p * p / n
    w = (4**p / lam) * (math.factorial(p) / n**p) ** 2
    return {
        "n": n, "p": p, "d": dim, "cs": cs, "labels": labels,
        "matrices": matrices, "col": np.diagonal(col, axis1=1, axis2=2),
        "row": np.diagonal(row, axis1=1, axis2=2), "offdiag": offdiag,
        "physical_lambda_N": lam, "physical_w_N_Js_1": w,
    }


def cov_bound(n: int, p: int, m: int, ti: frozenset, tj: frozenset, side: str) -> int:
    if not 0 <= m <= n or m < p:
        return 0
    if side == "row":
        ti, tj = tj, ti
    return choose(m - len(tj), p - len(tj)) * choose(n - m + p - len(ti), p - len(ti))


def cov_exact(n: int, p: int, occupation: frozenset, ti: frozenset, tj: frozenset, side: str) -> int:
    if side == "row":
        ti, tj = tj, ti
    m = len(occupation)
    if m < p or not tj <= occupation:
        return 0
    needed = tj | (ti & occupation)
    return choose(m - len(needed), p - len(needed)) * choose(n - m + p - len(ti), p - len(ti))


def pairings(points: tuple[int, ...]):
    if not points:
        yield ()
        return
    first = points[0]
    for j in range(1, len(points)):
        other = points[j]
        rest = points[1:j] + points[j + 1:]
        for pairing in pairings(rest):
            yield ((first, other),) + pairing


def crossings(pairing) -> int:
    return sum(a < c < b < d or c < a < d < b for (a, b), (c, d) in itertools.combinations(pairing, 2))


def psd_exact(matrix: np.ndarray) -> tuple[bool, str]:
    """Rational Schur complements, including the necessary zero-pivot row test."""
    a = [[Fraction(int(v)) for v in row] for row in matrix]
    smallest = None
    for k in range(len(a)):
        pivot = a[k][k]
        smallest = pivot if smallest is None else min(smallest, pivot)
        if pivot < 0:
            return False, str(pivot)
        if pivot == 0:
            if any(a[k][j] != 0 for j in range(k + 1, len(a))):
                return False, "zero pivot with nonzero row"
            continue
        for i in range(k + 1, len(a)):
            for j in range(i, len(a)):
                a[i][j] -= a[i][k] * a[k][j] / pivot
                a[j][i] = a[i][j]
    return True, str(smallest if smallest is not None else 0)


def build_network(sysdata: dict, pairing, mode: str) -> dict:
    n, p, dim = sysdata["n"], sysdata["p"], sysdata["d"]
    length = 2 * len(pairing)
    endpoint = {}
    domains, base_tags, constraint_records = [], [], []
    for edge, (left, right) in enumerate(pairing):
        # Alternating orientation exercises both A and A^dagger opening/closing.
        trans_left = bool(edge % 2)
        restrictions = [(frozenset(), frozenset()), (frozenset(), frozenset())]
        if mode == "tagged" and edge == 0:
            restrictions = [(frozenset({0}), frozenset({0})), (frozenset({0}), frozenset())]
        if mode == "union" and edge == 0:
            restrictions = [(frozenset({0}), frozenset()), (frozenset(), frozenset({n - 1}))]
        if mode == "union" and edge == 1:
            restrictions = [(frozenset(), frozenset({0})), (frozenset({0}), frozenset({0}))]
        ti, tj = frozenset(), frozenset()
        for (ri, rj), trans in zip(restrictions, (trans_left, not trans_left)):
            bi, bj = (rj, ri) if trans else (ri, rj)
            ti, tj = ti | bi, tj | bj
        domain = []
        for k, (i_set, j_set) in enumerate(sysdata["labels"]):
            direct_ok = True
            for (ri, rj), trans in zip(restrictions, (trans_left, not trans_left)):
                ei, ej = (j_set, i_set) if trans else (i_set, j_set)
                direct_ok &= ri <= ei and rj <= ej
            union_ok = ti <= i_set and tj <= j_set
            if direct_ok != union_ok:
                raise AssertionError("two-endpoint tag union mismatch")
            if union_ok:
                domain.append(k)
        domains.append(domain)
        base_tags.append((ti, tj))
        constraint_records.append({
            "edge": list(pairing[edge]), "left_transpose": trans_left,
            "endpoint_tags": [[sorted(a), sorted(b)] for a, b in restrictions],
            "base_union_TI": sorted(ti), "base_union_TJ": sorted(tj), "domain_size": len(domain),
        })
        endpoint[left] = (edge, trans_left)
        endpoint[right] = (edge, not trans_left)
    sandwiches = [np.eye(dim, dtype=np.int64) for _ in range(length + 1)]
    shifts = [0] * (length + 1)
    sandwich_names = ["I"] * (length + 1)
    if mode in ("sandwich", "tagged", "union"):
        # Genuine CAR contractions; the intervening charge is m-1.
        sandwiches[0] = sysdata["cs"][0].T
        sandwiches[length // 2] = sysdata["cs"][0]
        shifts[0], shifts[length // 2] = 1, -1
        sandwich_names[0], sandwich_names[length // 2] = "c0_dagger", "c0"
    return {
        "endpoint": endpoint, "domains": domains, "base_tags": base_tags,
        "constraints": constraint_records, "sandwiches": sandwiches,
        "shifts": shifts, "sandwich_names": sandwich_names, "length": length,
    }


def direct_network(sysdata: dict, network: dict) -> tuple[np.ndarray, int]:
    answer = np.zeros((sysdata["d"], sysdata["d"]), dtype=np.int64)
    count = 0
    for assignment in itertools.product(*network["domains"]):
        product = network["sandwiches"][0].copy()
        for position in range(network["length"]):
            edge, trans = network["endpoint"][position]
            a = sysdata["matrices"][assignment[edge]]
            product = product @ (a.T if trans else a) @ network["sandwiches"][position + 1]
        answer += product
        count += 1
    return answer, count


def ancilla_network(sysdata: dict, network: dict) -> tuple[np.ndarray, int, int]:
    if any(not domain for domain in network["domains"]):
        return np.zeros((sysdata["d"], sysdata["d"]), dtype=np.int64), 0, 0
    tensor = network["sandwiches"][-1].copy()
    active = []
    permutes, peak_entries = 0, tensor.size
    for position in reversed(range(network["length"])):
        edge, trans = network["endpoint"][position]
        mats = [sysdata["matrices"][i].T if trans else sysdata["matrices"][i] for i in network["domains"][edge]]
        if edge not in active:
            tensor = np.stack([np.tensordot(a, tensor, axes=(1, 0)) for a in mats], axis=-2)
            active.append(edge)
        else:
            axis = 1 + active.index(edge)
            permutes += int(active[-1] != edge)
            tensor = np.moveaxis(tensor, axis, -2)
            tensor = sum((np.tensordot(a, tensor[..., i, :], axes=(1, 0)) for i, a in enumerate(mats)))
            active.remove(edge)
        tensor = np.tensordot(network["sandwiches"][position], tensor, axes=(1, 0))
        peak_entries = max(peak_entries, tensor.size)
    if active:
        raise AssertionError("unclosed ancilla labels")
    return tensor, permutes, peak_entries


def endpoint_bound_squared(sysdata: dict, network: dict, m: int) -> tuple[int, list]:
    charge = m + network["shifts"][-1]
    active = set()
    result, factors = 1, []
    for position in reversed(range(network["length"])):
        edge, trans = network["endpoint"][position]
        side = "column" if edge not in active else "row"
        actual_side = ("row" if side == "column" else "column") if trans else side
        ti, tj = network["base_tags"][edge]
        factor = cov_bound(sysdata["n"], sysdata["p"], charge, ti, tj, actual_side)
        factors.append({"position": position, "charge": charge, "operation": side, "base_covariance": actual_side, "bound": factor})
        result *= factor
        active.symmetric_difference_update({edge})
        charge += network["shifts"][position]
    return result, factors


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(f"Sealed output exists; refusing to rerun: {OUTPUT}")
    started = time.perf_counter()
    systems = {(n, p): system(n, p) for n in range(2, 6) for p in (1, 2)}
    failures = []
    cov_records = []
    cov_count = same_site_count = impossible_count = strict_count = 0
    for (n, p), data in systems.items():
        tags = [frozenset(s) for s in subsets(n, 2)]
        max_difference = max_excess = 0
        local_count = 0
        for ti, tj in itertools.product(tags, repeat=2):
            indices = [k for k, (i, j) in enumerate(data["labels"]) if ti <= i and tj <= j]
            for side in ("column", "row"):
                observed = data["col" if side == "column" else "row"][indices].sum(axis=0)
                for state in range(data["d"]):
                    occupied = frozenset(i for i in range(n) if state & (1 << i))
                    exact = cov_exact(n, p, occupied, ti, tj, side)
                    bound = cov_bound(n, p, len(occupied), ti, tj, side)
                    actual = int(observed[state])
                    delta, excess = abs(actual - exact), actual - bound
                    max_difference, max_excess = max(max_difference, delta), max(max_excess, excess)
                    cov_count += 1
                    local_count += 1
                    same_site_count += bool(ti & tj)
                    impossible_count += not indices
                    strict_count += actual < bound
                    if delta or excess > 0:
                        failures.append({"type": "covariance", "N": n, "p": p, "state": state, "TI": sorted(ti), "TJ": sorted(tj), "side": side, "actual": actual, "exact": exact, "bound": bound})
        if data["offdiag"]:
            failures.append({"type": "nonzero_covariance_offdiagonal", "N": n, "p": p, "error": data["offdiag"]})
        cov_records.append({"N": n, "p": p, "labels": len(data["labels"]), "tags_per_side": len(tags), "comparisons": local_count, "exact_formula_max_abs_error": max_difference, "upper_bound_max_positive_excess": max_excess, "individual_covariance_max_offdiagonal": data["offdiag"], "lambda_N": data["physical_lambda_N"], "physical_w_N_for_Js_1": data["physical_w_N_Js_1"]})

    configs = [(2, 1, k, mode) for k in (2, 3, 4) for mode in ("plain", "sandwich", "tagged", "union")]
    configs += [(3, 1, 3, mode) for mode in ("plain", "sandwich", "tagged", "union")]
    configs += [(4, 2, 2, mode) for mode in ("plain", "sandwich", "tagged", "union")]
    network_records = []
    totals = {"networks": 0, "crossed_networks": 0, "assignments": 0, "ancilla_permutations": 0, "sector_psd_checks": 0, "nonzero_networks": 0, "nonzero_crossed_sandwich_networks": 0}
    max_ancilla_error = max_bound_ratio = 0.0
    for n, p, k, mode in configs:
        data = systems[n, p]
        for pairing in pairings(tuple(range(2 * k))):
            network = build_network(data, pairing, mode)
            direct, assignments = direct_network(data, network)
            ancilla, permutations, peak_entries = ancilla_network(data, network)
            error = int(np.max(np.abs(direct - ancilla)))
            max_ancilla_error = max(max_ancilla_error, error)
            if error:
                failures.append({"type": "ancilla", "N": n, "p": p, "pairing": pairing, "mode": mode, "error": error})
            cross = crossings(pairing)
            nonzero = bool(np.any(direct))
            totals["networks"] += 1
            totals["crossed_networks"] += bool(cross)
            totals["assignments"] += assignments
            totals["ancilla_permutations"] += permutations
            totals["nonzero_networks"] += nonzero
            totals["nonzero_crossed_sandwich_networks"] += nonzero and bool(cross) and mode != "plain"
            sectors = []
            for m in range(n + 1):
                columns = [state for state in range(data["d"]) if state.bit_count() == m]
                rows = [state for state in range(data["d"]) if state.bit_count() == m + sum(network["shifts"])]
                block = direct[np.ix_(rows, columns)]
                bound_squared, factors = endpoint_bound_squared(data, network, m)
                gram = block.T @ block
                witness = bound_squared * np.eye(len(columns), dtype=np.int64) - gram
                passed, min_pivot = psd_exact(witness)
                totals["sector_psd_checks"] += 1
                opnorm = float(np.linalg.norm(block, ord=2)) if block.size else 0.0
                ratio = opnorm / math.sqrt(bound_squared) if bound_squared else (0.0 if opnorm == 0 else None)
                if ratio is not None:
                    max_bound_ratio = max(max_bound_ratio, ratio)
                if not passed:
                    failures.append({"type": "exact_operator_bound", "N": n, "p": p, "pairing": pairing, "mode": mode, "m": m, "bound_squared": bound_squared, "schur_pivot": min_pivot})
                sectors.append({"m": m, "input_dimension": len(columns), "bound_squared_exact": bound_squared, "operator_norm_float": opnorm, "norm_over_bound": ratio, "exact_psd_pass": passed, "minimum_schur_pivot": min_pivot, "endpoint_factors": factors, "normalized_trace_exact": str(Fraction(int(np.trace(block)), len(columns))) if block.shape[0] == block.shape[1] else None})
            network_records.append({"N": n, "p": p, "pairs": k, "pairing": pairing, "crossings": cross, "mode": mode, "endpoint_constraints": network["constraints"], "sandwiches": network["sandwich_names"], "assignments": assignments, "nonzero": nonzero, "ancilla_permutations": permutations, "peak_ancilla_tensor_entries": peak_entries, "direct_vs_ancilla_max_abs_error": error, "sectors": sectors})

    payload = {
        "stage": "susskind-bilocal-fourth-tail", "check": "tagged-covariance-and-ancilla-networks", "success": not failures,
        "scope": "Exact finite CAR fixtures only. Uniform double-scaling degree bounds require the separate analytic tagged proof.",
        "conventions": {"CAR": "{c_i,c_j^dagger}=delta_ij", "C_I_dagger": "ascending ordered product of creators", "A_IJ": "C_I_dagger C_J", "I_J_overlap": "allowed", "coupling_covariance": "E J_IJ J_KL=delta_IL delta_JK", "fixture_variance": 1, "physical_w_N": "4^p Js^2/lambda_N (p!/N^p)^2", "variance_rescaling": "a 2k-endpoint pairing and its norm bound both multiply by w_N^k", "tags": "all TI,TJ cardinality <=2 for covariance; each endpoint restriction is intersected before comparing to its union in the base edge orientation", "norm": "operator norm on fixed occupation sectors; normalized traces also recorded", "proof_of_psd": "exact Fraction Schur complements with zero-pivot row test", "sampling": "none; all finite labels/pairings in stated configurations enumerated"},
        "summary": {"covariance_comparisons": cov_count, "overlapping_TI_TJ_comparisons": same_site_count, "empty_restricted_family_comparisons": impossible_count, "strict_covariance_bound_comparisons": strict_count, **totals, "max_direct_ancilla_abs_error": max_ancilla_error, "max_operator_norm_over_bound_float": max_bound_ratio, "failure_count": len(failures)},
        "covariance_systems": cov_records, "network_configurations": configs, "networks": network_records,
        "failures": failures, "runtime_seconds": time.perf_counter() - started,
        "versions": {"python": sys.version, "numpy": np.__version__, "platform": platform.platform()},
        "input_sha256": {str(path.relative_to(ROOT)).replace("\\", "/"): sha(path) for path in (Path(__file__).resolve(), ROOT / "notes/susskind-bilocal-four-point-transfer.md", ROOT / "literature/susskind-bilocal-fourth-tail-sources.md")},
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({"output": str(OUTPUT), "success": payload["success"], "summary": payload["summary"], "runtime_seconds": payload["runtime_seconds"], "sha256": sha(OUTPUT)}, ensure_ascii=False))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
