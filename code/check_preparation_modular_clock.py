"""One-shot independent Jacobi/product-integral comparison; no old drivers run."""
from pathlib import Path
import hashlib, json, platform, time
import numpy as np
import scipy
from preparation_modular_clock import (parameters, angular, jacobi, spectrum,
                                      correlation, sum_rules, product_log_error, window_rates)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/susskind-preparation-modular-clock/clock-checks.json'
if OUT.exists():
    raise FileExistsError(OUT)
start = time.monotonic()
cc = lambda z: [float(complex(z).real), float(complex(z).imag)]
rows, windows, summaries = [], [], []
for q in (.25, .5):
    specs = {f'{method}{L}': spectrum(q, fun(q, L))
             for method, fun in [('jacobi', jacobi), ('angular', angular)] for L in (128, 192)}
    j, a = specs['jacobi192'], specs['angular192']
    for s in (0., .5, 1., 2., 4., .5-.5j, 1-.5j, 2-.5j, .5-1j, 1-1j, 2-1j):
        vals = {name: correlation(spec, s) for name, spec in specs.items()}
        cross = abs(vals['jacobi192'][0]-vals['angular192'][0])
        refine = max(abs(vals[f'{method}192'][0]-vals[f'{method}128'][0])
                     for method in ('jacobi', 'angular'))
        assert cross < 2e-10 and refine < 2e-10, (q, s, cross, refine)
        rows.append(dict(q=q, s=cc(s), values={k: cc(v[0]) for k,v in vals.items()},
                         cross_absolute=cross, refinement_absolute=refine,
                         cross_relative=cross/max(abs(vals['angular192'][0]), 1e-300),
                         cancellation=max(v[1] for v in vals.values())))
    kms = []
    for t in (.5, 1., 2.):
        err = abs(correlation(a, t-1j)[0]-correlation(a, -t)[0])
        assert err < 1e-12
        kms.append(dict(t=t, residual=err))
    sr = {name: sum_rules(spec) for name, spec in specs.items()}
    for key in ('mass', 'first_frequency', 'second_frequency', 'relative_entropy'):
        assert abs(sr['jacobi192'][key]-sr['angular192'][key]) < 2e-10, (q,key)
    assert abs(sr['angular192']['mass']-1) < 1e-12
    assert sr['angular192']['first_frequency'] > 0 and sr['angular192']['relative_entropy'] > 0
    lam, aa, bb, R = parameters(q)
    assert aa > 0 and abs(aa+bb-1)<1e-14
    plateau = (1+bb*(1+q))/(1+lam)
    modular_h0 = plateau+lam/(1+lam)
    assert abs(modular_h0-(1+np.sqrt(q)*(1+q))) < 1e-14
    summaries.append(dict(q=q, lambda_=lam, a=aa, b=bb, radius=R,
                          bandwidth=float(np.log((aa+bb*R*R)/aa)), sum_rules=sr, kms_checks=kms,
                          modular_h_plateau=plateau, modular_h_zero=modular_h0,
                          product_log_tail=product_log_error(q,q,192)))
    for interval in ((.8,1.2),(-1.2,-.8)):
        coarse = window_rates(q, interval, (-.2,.2), 32)
        fine = window_rates(q, interval, (-.2,.2), 64)
        floor = (1-np.sqrt(q)+np.sqrt(q)*.8**2)/(1-np.sqrt(q)+np.sqrt(q)*.2**2)
        assert fine['ratio'] >= floor > 1
        assert max(abs(fine[k]-coarse[k]) for k in fine) < 1e-12
        windows.append(dict(q=q, output_interval=interval, input_interval=(-.2,.2),
                            coarse=coarse, fine=fine, analytic_ratio_lower=floor))
    assert abs(windows[-1]['fine']['ratio']-windows[-2]['fine']['ratio']) < 1e-12

inputs = ['code/preparation_modular_clock.py', 'code/check_preparation_modular_clock.py']
result = dict(status='passed', rows=rows, summaries=summaries, energy_windows=windows,
              environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                               arithmetic='float64/complex128'),
              parameters=dict(q=[.25,.5], quadrature_or_Jacobi=[128,192], products=192,
                              window_orders=[32,64], physical_lambda='-log(q)/4'),
              scope='Exact modular correlation and energy-filtered charged KMS witness. Independent representations plus refinement are diagnostics, not total-error certificates. Product-tail bound excludes quadrature/roundoff. KMS exchange checks are implementation sign checks, not independent proofs.',
              seconds=time.monotonic()-start,
              input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
OUT.parent.mkdir(parents=True,exist_ok=True)
with OUT.open('x',encoding='utf-8') as f:
    json.dump(result,f,indent=2); f.write('\n')
print(json.dumps(dict(status='passed', rows=len(rows), windows=len(windows),
                     max_cross=max(r['cross_absolute'] for r in rows),
                     max_refinement=max(r['refinement_absolute'] for r in rows),
                     sha256=hashlib.sha256(OUT.read_bytes()).hexdigest(), seconds=result['seconds'])))
