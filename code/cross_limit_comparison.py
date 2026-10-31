"""Projected flat-energy comparison only; does not evaluate the exact 6j limit.

Independent survival-function and ordered-path calculations of the static
boundary commutator. Gauss refinements are empirical, not interval bounds.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import platform
import sys
sys.dont_write_bytecode = True
import numpy as np
import scipy
from scipy.special import roots_legendre
from semiclassical_clock import stationary_transition

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/stage3/cross-limit'
NU, SIGMA, WIDTH = .6, 1., .2

def rule(lo, hi, order, breaks=()):
    if hi <= lo:
        return np.empty(0), np.empty(0)
    # Subdivision controls smooth kernels; physical discontinuity surfaces
    # are supplied explicitly, rather than sampled on an indicator grid.
    edges = sorted(set([float(lo), float(hi)] +
        [float(x) for x in breaks if lo < x < hi] +
        list(np.linspace(lo, hi, int(np.ceil(hi-lo))+1))))
    x, w = roots_legendre(order)
    return (np.concatenate([(a+b)/2+(b-a)*x/2 for a,b in zip(edges[:-1],edges[1:])]),
            np.concatenate([(b-a)*w/2 for a,b in zip(edges[:-1],edges[1:])]))

def k(x):
    x = np.asarray(x)
    return stationary_transition(x, NU)*np.exp(-(SIGMA*x)**2)

def survival(q, start):
    t = np.clip((np.asarray(q)-start)/WIDTH, 0., 1.)
    return 1-t+2*np.sin(2*np.pi*t)/(3*np.pi)-np.sin(4*np.pi*t)/(12*np.pi)

def survival_moments(start, order, radius):
    x, wx = rule(0., radius, order, [start,start+WIDTH,radius-start,radius-start-WIDTH])
    gain = np.dot(wx,k(x))
    lk = np.dot(wx,k(-x)*survival(x,start))
    mixed, losses, diff = [], [], []
    for xx in x:
        y, wy = rule(0.,radius,order,[start,start+WIDTH,xx+start,xx+start+WIDTH])
        shifted = survival(np.maximum(0.,y-xx),start)
        base = survival(y,start)
        mixed.append(np.dot(wy,k(-y)*shifted))
        diff.append(np.dot(wy,k(-y)*(shifted-base)))
        z, wz = rule(0.,radius,order,[start-xx,start+WIDTH-xx])
        losses.append(np.dot(wz,k(-z)*survival(xx+z,start)))
    ll = np.dot(wx,k(-x)*np.asarray(losses))
    mix = np.dot(wx,k(x)*np.asarray(mixed))
    n = gain*gain + gain*lk + mix + ll
    cross = gain*gain + 2*gain*lk + ll
    floor = 2*np.dot(wx,k(x)*np.asarray(diff))
    return dict(N=float(n),X=float(cross),floor_survival=float(floor),
                floor_subtracted=float(2*(n-cross)),gain_integral=float(gain))

def direct_path_floor(start, order, radius):
    r, wr = rule(start,start+WIDTH,order)
    vals = []
    for rr in r:
        x, wx = rule(0.,radius,order,[radius-rr])
        losses=[]
        for xx in x:
            y, wy = rule(rr,min(rr+xx,radius),order)
            losses.append(np.dot(wy,k(-y)))
        vals.append(2*np.dot(wx,k(x)*np.asarray(losses)))
    density=8*np.sin(np.pi*(r-start)/WIDTH)**4/(3*WIDTH)
    return float(np.dot(wr,density*np.asarray(vals)))

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    rows=[]
    checks=[]
    for start in [0.,.25,.5,1.,2.,3.]:
        runs=[]
        for order,radius in [(12,6.),(20,6.),(28,8.)]:
            row=survival_moments(start,order,radius)
            row.update(order=order,radius=radius)
            runs.append(row)
        direct=direct_path_floor(start,24,8.)
        final=runs[-1]
        rel=abs(direct-final['floor_survival'])/max(abs(direct),1e-300)
        refinement=abs(runs[-1]['floor_survival']-runs[-2]['floor_survival'])/abs(final['floor_survival'])
        checks.append(dict(name=f'independent_path_vs_survival_{start}',value=rel,limit=2e-8,passed=rel<2e-8))
        checks.append(dict(name=f'quadrature_refinement_{start}',value=refinement,limit=2e-8,passed=refinement<2e-8))
        checks.append(dict(name=f'positive_floor_below_twice_norm_{start}',passed=bool(0<direct<2*final['N'])))
        rows.append(dict(packet_start=start,runs=runs,direct_path_floor=direct,relative_comparison=rel))
    # Verify the packet normalization separately on its original variable.
    r,wr=rule(0.,WIDTH,32)
    norm=float(np.dot(wr,8*np.sin(np.pi*r/WIDTH)**4/(3*WIDTH)))
    checks.append(dict(name='packet_normalization',value=abs(norm-1),limit=1e-13,passed=abs(norm-1)<1e-13))
    paths=['code/cross_limit_comparison.py','code/semiclassical_clock.py','code/quantized_observer.py',
           'manuscript/crosslimit.tex','sources/kolchmeyer-liu-2411.08090v1.html']
    output=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
        kind='comparison-model numerical study; no exact finite-M 6j convergence claim',
        status='passed' if all(c['passed'] for c in checks) else 'failed',
        parameters=dict(nu=NU,sigma=SIGMA,packet_width=WIDTH,packet='8/(3w) sin^4(pi(r-r0)/w)',
                        variable='exact proper-energy excess r=m-M; not modular p=2*pi*r'),
        versions=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform()),
        rows=rows,checks=checks,
        error_scope='Gauss order and finite-radius comparisons are empirical, not a certified total error. Positive formula avoids cancellation in N-X.',
        hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})
    (OUT/'comparison.json').write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=output['status'],checks=len(checks),rows=[dict(start=r['packet_start'],floor=r['direct_path_floor'],relative_difference=r['relative_comparison']) for r in rows]),indent=2))
    if output['status']!='passed':
        raise SystemExit(1)

if __name__=='__main__': main()
