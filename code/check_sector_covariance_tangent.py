"""Exact second-order path identities and scalar remainder diagnostics."""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
import hashlib
import json
import math
import platform

import sector_covariance_tangent as t

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/susskind-sector-neutral-transfer/tangent-checks.json'


def main():
    if OUT.exists():
        raise SystemExit('Successful scientific output exists; do not overwrite.')
    counts = {'matrix_units': 0, 'path_equalities': 0, 'norm_bounds': 0,
              'nonzero_remainders': 0, 'nonzero_first_terms': 0,
              'zero_reference_cases': 0}
    rows = []
    for n in (2, 3, 4):
        before = dict(counts)
        for row, col in product(range(1 << n), repeat=2):
            if row.bit_count() != col.bit_count():
                continue
            x = {(row, col): F(1)}
            counts['matrix_units'] += 1
            for p in range(1, min(2, n)+1):
                lhs_full = t.path_apply(x, n, ['T']*p+['S']*p)
                for m in range(p, n+1):
                    _, den, _ = t.path_data(n, p, m, F(0))
                    for b in (F(0), F(1, 2), F(min(den))):
                        rhs, first, rem, coeff = t.tangent_expansion(x, n, p, m, b)
                        assert t.restrict(lhs_full, m) == rhs, (n, row, col, p, m, b)
                        counts['path_equalities'] += 1
                        counts['nonzero_remainders'] += bool(rem)
                        counts['nonzero_first_terms'] += bool(first)
                        counts['zero_reference_cases'] += b == min(den)
                        second_norms = []
                        for a, c in product(('T', 'S'), repeat=2):
                            y = t.deficit(t.deficit(x, n, b, c), n, b, a)
                            second_norms += [t.norm2_squared(y, n, k) for k in range(n+1)]
                        bound = F(math.comb(2*p, 2), min(den)**2)**2 * max(second_norms)
                        assert t.norm2_squared(rem, n, m) <= bound
                        counts['norm_bounds'] += 1
        rows.append({'N': n, 'counts': {k: counts[k]-before[k] for k in counts}})

    scalars = []
    for p in (4, 8, 16, 32, 64):
        n = 4*p*p
        for z, v in product((-2, 0, 2), (1, 2, 3, 5)):
            m = n//2+z*p
            _, den, _ = t.path_data(n, p, m, v*p)
            if min(den) <= v*p:
                continue
            d = t.scalar_tangent(n, p, m, v)
            for key in ('zeta_T', 'zeta_S'):
                assert abs(d[key]-d['target_each']) <= d['coefficient_error_bound']
            scalars.append({'N': n, 'p': p, 'central_offset': z, 'v': v,
                            'exact': {k: str(a) for k, a in d.items()},
                            'display': {k: float(a) for k, a in d.items()},
                            'p_scaled_zeta_errors': {
                                k: float(p*abs(d[k]-d['target_each'])) for k in ('zeta_T', 'zeta_S')}})
    assert counts['nonzero_remainders'] > 0 and counts['zero_reference_cases'] > 0
    result = {'status': 'passed', 'counts': counts, 'rows': rows, 'scalar_rows': scalars,
              'parameters': {'N_exact': [2, 3, 4], 'p_exact_max': 2,
                             'b': ['0', '1/2', 'minimum path denominator'],
                             'scalar_p': [4, 8, 16, 32, 64], 'scalar_lambda': '1/4'},
              'environment': {'python': platform.python_version(), 'arithmetic': 'fractions.Fraction'},
              'input_sha256': {
                  'code/check_sector_covariance_tangent.py': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  'code/sector_covariance_tangent.py': hashlib.sha256(Path(t.__file__).read_bytes()).hexdigest()},
              'scope': 'Full neutral matrix-unit basis for exact path/first/second-deficit identity, including theta=0, and exact S2 remainder upper bounds; scalar finite-parameter coefficient bounds. Not a numerical simulation of the original Gaussian ensemble or a proof of asymptotic field convergence.'}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open('x', encoding='utf-8') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    print(json.dumps({'status': 'passed', 'counts': counts, 'scalar_rows': len(scalars),
                      'sha256': hashlib.sha256(OUT.read_bytes()).hexdigest()}))


if __name__ == '__main__':
    main()
