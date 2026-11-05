"""Read-only provenance and exact-endpoint checks, not an analytic proof checker."""
from pathlib import Path
from fractions import Fraction
import hashlib
import json


def audit(root=None):
    root = Path(root) if root else Path(__file__).resolve().parents[1]
    prior_path = root / 'results/stage4/residual-certificate/evidence.json'
    chain_path = root / 'results/stage4/residual-certificate/tensor-chain-audit.json'
    prior = json.loads(prior_path.read_text(encoding='utf-8'))
    chain = json.loads(chain_path.read_text(encoding='utf-8'))
    errors = []
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    for name, expected in prior['hashes'].items():
        if digest(root / name) != expected:
            errors.append('prior bound file changed: ' + name)
    endpoint = lambda obj, side: Fraction(*map(int, obj[side]))
    separation_lower = endpoint(chain['separation'], 'lower')
    error_upper = endpoint(chain['total_step_error_upper'], 'upper')
    saved_lower = endpoint(chain['word_norm_lower_bound'], 'lower')
    threshold = Fraction(4489, 100000)
    if not chain['passed'] or chain['failed']:
        errors.append('prior chain records a failure')
    if not threshold < saved_lower <= separation_lower - error_upper:
        errors.append('conservative lower endpoint does not support threshold')
    if prior['word_norm_lower_bound'] != chain['word_norm_lower_bound']:
        errors.append('prior evidence and chain bound differ')
    return dict(
        status='passed' if not errors else 'failed', errors=errors,
        prior_bound_files_checked=len(prior['hashes']),
        threshold=[str(threshold.numerator), str(threshold.denominator)],
        neighborhood_threshold=[str((threshold/2).numerator), str((threshold/2).denominator)],
        interval_radius=None,
        scope='Hash and exact rational endpoint checks only. No overlap reintegration, '
              'no formal proof checking, no quantitative parameter neighborhood.',
        hashes={p.relative_to(root).as_posix(): digest(p)
                for p in (prior_path, chain_path, Path(__file__).resolve())})


if __name__ == '__main__':
    result = audit()
    print(json.dumps(result, indent=2))
    raise SystemExit(result['status'] != 'passed')
