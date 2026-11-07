"""Complete Gaussian-Hermite checks of the finite neutral marked series.

All coupling coordinates and matrix entries are retained.  Hermite dictionary
keys label unnormalised probabilists' Hermite polynomials, so their annealed
norm has the multi-index factorial weight used in memory_jacobi_bridge.
The output is exclusive-create; these small CAR identities are not a DS limit.
"""
from pathlib import Path
from fractions import Fraction
from math import sqrt
import datetime
import hashlib
import json
import platform
import time

import numpy as np
import sympy as sp

from memory_jacobi_bridge import add, norm2, raw_product, step
from memory_charged_transfer import system
from memory_wick_rank import mu


ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results/susskind-memory-neutral-marked-series"
PARAMETERS = [(3, 1, 2, 6), (4, 1, 2, 5), (4, 2, 3, 4)]
ATOL = 5e-10
RTOL = 5e-10


def scale(poly, value):
    return {key: value * mat for key, mat in poly.items()}


def degree(poly, value):
    return {key: mat for key, mat in poly.items() if sum(key) == value}


def norm(poly):
    return sqrt(max(0.0, norm2(poly)))


def compare(left, right):
    absolute = norm(add(left, right, -1))
    reference = max(norm(left), norm(right))
    assert absolute <= ATOL + RTOL * reference, (absolute, reference)
    return dict(absolute_error=absolute, reference_norm=reference,
                relative_error=absolute / reference if reference > ATOL else None,
                tolerance=ATOL + RTOL * reference)


def right_step(poly, basis):
    # Gaussian coordinates and Hermite polynomials are real, and every B_i is
    # Hermitian.  Adjoints therefore turn the imported left step into X H.
    adjoint = {key: mat.conj().T for key, mat in poly.items()}
    return {key: mat.conj().T for key, mat in step(adjoint, basis).items()}


def check(N, p, m, max_degree):
    started = time.time()
    sys = system(N, p, m)
    basis, dim = sys["up"], sys["dp"]
    zero = (0,) * len(basis)
    identity = {zero: np.eye(dim, dtype=complex)}
    theta_exact = [mu(N, p, m, j * p) for j in range(max_degree + 1)]
    theta = list(map(float, theta_exact))
    a_exact = [sum(theta_exact[:j], Fraction(0)) for j in range(max_degree + 1)]
    a = list(map(float, a_exact))
    # Row vectorisation: vec_row(B X B) = (B tensor B^T) vec_row(X).
    phi_matrix = sum(np.kron(mat, mat.T) for mat in basis)

    def phi(poly):
        return {key: (phi_matrix @ mat.reshape(-1)).reshape(dim, dim)
                for key, mat in poly.items()}

    def theta_number(poly):
        return {key: theta[sum(key)] * mat for key, mat in poly.items()}

    covariance_error = compare(phi(identity), identity)
    ordinary, wick, jacobi = [identity], [identity], [identity]
    for j in range(max_degree):
        ordinary.append(step(ordinary[-1], basis))
        wick.append(step(wick[-1], basis, "creation"))
        lower = jacobi[j - 1] if j else {}
        jacobi.append(add(step(jacobi[-1], basis), lower, -a[j]))
    errors = [add(w, q, -1) for w, q in zip(wick, jacobi)]
    eta = [scale(err, p) for err in errors]

    # Monomial-to-Jacobi conversion uses rational coefficients and the scalar
    # recurrence x Q_j = Q_{j+1} + a_j Q_{j-1}; it does not fit sample values.
    polynomial_rows = []
    coefficients = [Fraction(1)]
    for n in range(max_degree + 1):
        reconstructed, conserved, correction = {}, {}, {}
        for j, exact in enumerate(coefficients):
            value = float(exact)
            reconstructed = add(reconstructed, jacobi[j], value)
            conserved = add(conserved, jacobi[j], value * theta[j])
            marked_error = add(scale(errors[j], theta[j]),
                               theta_number(errors[j]), -1)
            correction = add(correction, marked_error, value)
        first = add(phi(ordinary[n]), theta_number(ordinary[n]), -1)
        right = add(first, correction)
        left = add(phi(ordinary[n]), conserved, -1)
        polynomial_rows.append(dict(
            degree=n, ordinary_Hermite_matrices=len(ordinary[n]),
            highest_chaos_matrices=len(wick[n]),
            jacobi_coefficients=list(map(str, coefficients)),
            monomial_reconstruction=compare(ordinary[n], reconstructed),
            exact_neutral_decomposition=compare(left, right),
            conserved_energy_polynomial_norm=norm(conserved),
            rank_minus_chaos_term_norm=norm(first),
            lower_chaos_correction_norm=norm(correction),
            total_residual_norm=norm(left),
            forbidden_top_error_norm=norm({key: mat for key, mat in errors[n].items()
                                          if sum(key) > n - 2})))
        assert polynomial_rows[-1]["forbidden_top_error_norm"] < ATOL
        next_coefficients = [Fraction(0)] * (n + 2)
        for j, value in enumerate(coefficients):
            next_coefficients[j + 1] += value
            if j:
                next_coefficients[j - 1] += a_exact[j] * value
        coefficients = next_coefficients

    source_rows = []
    for j in range(1, max_degree):
        annihilation = step(wick[j], basis, "annihilation")
        residual = add(annihilation, wick[j - 1], -a[j])
        ordinary_source, wick_source = {}, {}
        for v in range(j):
            # This is literally W_{j-1} sum_v [Phi(H^v)-theta_v H^v]
            # H^{j-1-v}, using ordinary Hermite multiplication by H on the
            # right before taking the indicated coupling-chaos projection.
            word = add(phi(ordinary[v]), ordinary[v], -theta[v])
            for _ in range(j - 1 - v):
                word = right_step(word, basis)
            ordinary_source = add(ordinary_source, word)
            # Independent highest-chaos reconstruction with the imported
            # raw (Wick) product.  This excludes contractions by definition.
            block = add(phi(wick[v]), wick[v], -theta[v])
            wick_source = add(wick_source, raw_product(block, wick[j - 1 - v]))
        projected_source = degree(ordinary_source, j - 1)
        forced = add(step(eta[j], basis), eta[j - 1], -a[j])
        forced = add(forced, projected_source, -p)
        source_rows.append(dict(
            j=j, source_chaos_degree=j - 1, source_norm=norm(residual),
            annihilation_vs_ordinary_Wick_source=compare(residual, projected_source),
            ordinary_vs_raw_product_Wick_source=compare(projected_source, wick_source),
            forced_eta_recurrence=compare(eta[j + 1], forced),
            retained_unprojected_source_norm=norm(ordinary_source)))

    H, Z2 = ordinary[1], wick[2]
    G1 = add(phi(H), H, -theta[1])
    K = sum((phi_matrix @ mat.reshape(-1)).reshape(dim, dim) @ mat
            for mat in basis)
    delta_matrix = K - theta[1] * np.eye(dim)
    delta_scalar = np.trace(delta_matrix) / dim
    delta = {zero: delta_matrix}
    scalarity = compare(delta, {zero: delta_scalar * np.eye(dim)})
    explicit3 = scale(G1, -p)
    full4 = add(step(G1, basis), right_step(G1, basis))
    full4 = add(full4, add(phi(Z2), Z2, -theta[2]))
    full4 = add(full4, delta, -1)
    explicit4 = scale(full4, -p)
    low_orders = dict(
        formula_eta3="-p * (Phi(H) - theta_1 H)",
        formula_eta4="-p * (H G1 + G1 H + (Phi-theta_2) Z2 - delta I)",
        delta_definition="sum_i Phi(B_i) B_i - theta_1 I; C=1",
        delta_real=float(delta_scalar.real), delta_imag=float(delta_scalar.imag),
        contraction_scalar_identity=scalarity,
        eta3_explicit_identity=compare(eta[3], explicit3),
        eta4_explicit_identity=compare(eta[4], explicit4),
        eta4_scalar_identity=compare(degree(eta[4], 0), scale(delta, -p)),
        eta3_norm=norm(eta[3]), eta4_norm=norm(eta[4]),
        eta4_chaos_zero_norm=norm(degree(eta[4], 0)),
        eta4_chaos_two_norm=norm(degree(eta[4], 2)),
        error_if_eta3_energy_term_omitted=norm(add(eta[3], phi(H), p)),
        error_if_eta4_explicit_scalar_correction_omitted=norm(scale(delta, p)))
    assert low_orders["eta4_chaos_zero_norm"] > 1e-5
    return dict(N=N, p=p, m=m, C=1, charge_sector_dimension=dim,
                Gaussian_coordinates=len(basis), max_degree=max_degree,
                theta_exact=list(map(str, theta_exact)),
                a_exact=list(map(str, a_exact)),
                covariance_normalization=covariance_error,
                polynomial_decomposition=polynomial_rows,
                Wick_sources=source_rows, low_orders=low_orders,
                elapsed_seconds=time.time() - started)


if __name__ == "__main__":
    started = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    destination = OUT / "coefficient-checks.json"
    if destination.exists():
        raise FileExistsError(destination)
    rows = []
    for args in PARAMETERS:
        row = check(*args)
        rows.append(row)
        print(json.dumps(dict(system=args, status="passed",
                              elapsed_seconds=row["elapsed_seconds"],
                              eta4_scalar_norm=row["low_orders"]["eta4_chaos_zero_norm"])),
              flush=True)
    inputs = ["code/check_neutral_marked_series.py", "code/memory_jacobi_bridge.py",
              "code/memory_charged_transfer.py", "code/memory_wick_rank.py",
              "code/cubic_linear_dynamics.py", "code/complex_syk_dictionary.py"]
    parameters = dict(systems=PARAMETERS, absolute_tolerance=ATOL,
                      relative_tolerance=RTOL, covariance_C=1,
                      arithmetic="complex128 matrix coefficients; exact rational Jacobi scalars")
    parameter_bytes = json.dumps(parameters, sort_keys=True, separators=(",", ":")).encode()
    result = dict(
        created_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        status="passed", parameters=parameters,
        parameter_sha256=hashlib.sha256(parameter_bytes).hexdigest(), systems=rows,
        environment=dict(python=platform.python_version(), numpy=np.__version__,
                         sympy=sp.__version__, platform=platform.platform()),
        scope=("Complete Gaussian-Hermite coefficients for finite CAR systems; "
               "no Gaussian samples, coordinate truncation, energy pinching or discarded "
               "conserved terms. Checks finite polynomial neutral decomposition, "
               "ordinary-word Wick sources, forced eta recurrence and eta3/eta4 including "
               "nonzero scalar chaos. Floating-point checks are not exact-arithmetic "
               "certificates, uniform tail proofs or double-scaling convergence claims."),
        input_sha256={name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                      for name in inputs}, elapsed_seconds=time.time() - started)
    with destination.open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(dict(status="passed", systems=len(rows),
                          polynomial_cases=sum(len(row["polynomial_decomposition"]) for row in rows),
                          source_cases=sum(len(row["Wick_sources"]) for row in rows),
                          output_sha256=hashlib.sha256(destination.read_bytes()).hexdigest(),
                          elapsed_seconds=result["elapsed_seconds"])), flush=True)
