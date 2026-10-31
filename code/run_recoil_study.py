"""Generate reproducible stage-3 recoil data and a scientific figure."""
from __future__ import annotations

import csv
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import mpmath as mp

from recoil_model import (SOURCE_URL, SOURCE_EQUATIONS, recoil_hyper,
    recoil_series, recoil_angle, recoil_conformal, conformal_boundary,
    time_branch_q)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results' / 'stage3'
DPS = 60
TOL = '1e-45'
NU_GRID = [('0', mp.mpf(0)), ('0.5', mp.mpf('0.5')),
           ('1', mp.mpf(1)), ('2', mp.mpf(2)), ('i/2', mp.mpc(0, '0.5'))]
Q_GRID = ['0.01', '0.1', '0.3', '0.6', '0.8', '0.95']


def dec(x):
    return mp.nstr(x, 55)


def complex_record(x):
    return {'real': dec(mp.re(x)), 'imag': dec(mp.im(x))}


def write_csv(path, rows):
    with path.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    mp.mp.dps = DPS
    OUT.mkdir(parents=True, exist_ok=True)
    baseline = []
    for label, nu in NU_GRID:
        for q_str in Q_GRID:
            q = mp.mpf(q_str)
            h = recoil_hyper(nu, q)
            s = recoil_series(nu, q, TOL)
            a = recoil_angle(nu, q)
            baseline.append({
                'nu': label, 'q': q_str,
                'hyper_real': dec(mp.re(h)), 'hyper_imag': dec(mp.im(h)),
                'series_real': dec(mp.re(s['value'])),
                'angle_real': dec(mp.re(a)),
                'hyper_series_abs_error': dec(abs(h-s['value'])),
                'hyper_angle_abs_error': dec(abs(h-a)),
                'series_tail_bound': dec(s['tail_bound']),
                'series_terms': s['terms'],
                'conformal_abs_error': dec(abs(h-recoil_conformal(q))) if label == 'i/2' else '',
            })
    branches = []
    for q_str in ['1.05', '1.2', '2', '5']:
        q = mp.mpf(q_str)
        for sign_t in [-1, 1]:
            for contour in [-1, 1]:
                for eps_str in ['1e-3', '1e-5', '1e-7', '1e-9']:
                    t, qc = time_branch_q(q, '0.01', sign_t, contour, eps_str)
                    h = recoil_hyper(mp.mpc(0, '0.5'), qc)
                    elementary = recoil_conformal(qc)
                    boundary = conformal_boundary(q, contour*sign_t)
                    branches.append({
                        'q_real': q_str, 'time_sign': sign_t,
                        'contour': 'F12_upper_t' if contour == 1 else 'F14_lower_t',
                        'q_cut_side': contour*sign_t, 'time': dec(t),
                        'epsilon': eps_str, 'q_real_regulated': dec(mp.re(qc)),
                        'q_imag_regulated': dec(mp.im(qc)),
                        'hyper_real': dec(mp.re(h)), 'hyper_imag': dec(mp.im(h)),
                        'boundary_real': dec(mp.re(boundary)),
                        'boundary_imag': dec(mp.im(boundary)),
                        'elementary_abs_error': dec(abs(h-elementary)),
                        'boundary_abs_error': dec(abs(h-boundary)),
                    })
    write_csv(OUT/'baseline.csv', baseline)
    write_csv(OUT/'branches.csv', branches)
    report = {
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'source': {'url': SOURCE_URL, 'equations': SOURCE_EQUATIONS},
        'scope': 'Leading heavy-light pure-recoil dS3 formula only; not finite-clock dynamics.',
        'normalization': 'F/F0 in the fixed-average-mass, Boltzmann-weighted representation.',
        'precision_dps': DPS, 'series_absolute_tail_tolerance': TOL,
        'parameters': {'nu': [x[0] for x in NU_GRID], 'q': Q_GRID,
            'branch_nu': 'i/2', 'branch_heavy_to_observer_mass_ratio': '0.01',
            'branch_epsilon': ['1e-3', '1e-5', '1e-7', '1e-9']},
        'baseline': baseline, 'branches': branches,
        'max_hyper_series_error': dec(max(mp.mpf(x['hyper_series_abs_error']) for x in baseline)),
        'max_hyper_angle_error': dec(max(mp.mpf(x['hyper_angle_abs_error']) for x in baseline)),
        'max_branch_elementary_error': dec(max(mp.mpf(x['elementary_abs_error']) for x in branches)),
        'branch_warning': 'Both continuations are diagnostics; neither is selected as the physical state.',
        'environment': {'python': sys.version, 'platform': platform.platform(), 'mpmath': mp.__version__},
        'code_sha256': {name: sha(ROOT/'code'/name) for name in
                        ['recoil_model.py', 'run_recoil_study.py', 'validate_recoil.py']},
    }
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    report['environment']['matplotlib'] = matplotlib.__version__
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.1), constrained_layout=True)
    for label, _ in NU_GRID:
        rows = [row for row in baseline if row['nu'] == label]
        axes[0].plot([float(row['q']) for row in rows],
                     [float(row['hyper_real']) for row in rows], 'o-', label=r'$\nu=$'+label)
    axes[0].set(xlabel=r'$q=\nu^{\prime}\cosh(t)/(2m)$', ylabel=r'$\mathcal{F}/\mathcal{F}_0$',
                title='Pure recoil, before the cut', yscale='log')
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=.25)
    eps = ['1e-3', '1e-5', '1e-7', '1e-9']
    for q_str in ['1.05', '1.2', '2', '5']:
        rows = [row for row in branches if row['q_real']==q_str and row['time_sign']==1
                and row['contour']=='F12_upper_t']
        axes[1].loglog([float(row['epsilon']) for row in rows],
                       [float(row['boundary_abs_error']) for row in rows], 'o-', label='q='+q_str)
    axes[1].set(xlabel=r'time regulator $\epsilon$', ylabel='absolute error to upper boundary',
                title=r'Conformal light field, $\nu=i/2$')
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=.25, which='both')
    fig.savefig(OUT/'recoil-baseline.png', dpi=180)
    plt.close(fig)
    (OUT/'study.json').write_text(json.dumps(report, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    print(json.dumps({k: report[k] for k in ['max_hyper_series_error','max_hyper_angle_error',
        'max_branch_elementary_error']}, indent=2))
    print(f'Saved {len(baseline)} baseline rows and {len(branches)} branch diagnostics to {OUT}')


if __name__ == '__main__':
    main()
