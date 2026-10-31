"""Save phase, reflection and norm checks; no scalar unitarity test."""
import hashlib
import json
import platform
from pathlib import Path
import mpmath as mp
import sympy
from ds3_crossing_dictionary import (
    channel_data, period, physical_scalar, partial_wave)
from ds3_group_pairing import triple_group
from global_racah import scalar_6j, shadow_coefficient

ROOT = Path(__file__).resolve().parents[1]


def encode(x):
    if isinstance(x, dict): return {k: encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [encode(v) for v in x]
    if isinstance(x, mp.mpc): return dict(re=str(x.real), im=str(x.imag))
    if isinstance(x, mp.mpf): return str(x)
    return x


def run():
    mp.mp.dps = 45
    rows, checks = [], []
    def check(name, lhs, rhs, tol='1e-35'):
        err = abs(lhs-rhs)/max(abs(rhs), mp.mpf('1e-40'))
        checks.append(dict(name=name, lhs=lhs, rhs=rhs, error=err,
                           tolerance=tol, passed=bool(err < mp.mpf(tol))))
    for k, (labels, bs, ds) in enumerate([
        (['.4', '.7', '.93', '.6'], '.65', '.57'),
        (['1.1', '.43', '.82', '1.3'], '.37', '1.17'),
        (['.31', '1.2', '.66', '.91'], '1.4', '.29'),
    ]):
        e = list(map(mp.mpf, labels)); b, d = mp.mpf(bs), mp.mpf(ds)
        a, q, c, p = e
        ns, nt, ws, wt = channel_data(b, d, e)
        raw = scalar_6j(b, d, e)['value']
        u = physical_scalar(b, d, e)
        reverse = physical_scalar(d, b, [c, q, a, p])
        check(f'{k}: three-point physical quotient',
              abs(period(a,q,b))**2/triple_group([a,q,b]), 8*mp.pi**3)
        check(f'{k}: s-channel raw norm', abs(ns)**2/ws, 32*mp.pi**7/b**2)
        check(f'{k}: t-channel raw norm', abs(nt)**2/wt, 32*mp.pi**7/d**2)
        check(f'{k}: physical U from quotient',
              32*mp.pi**3*mp.sqrt(ws*wt)*raw/(mp.conj(ns)*nt), u)
        check(f'{k}: Hermitian channel exchange', u, mp.conj(reverse))
        check(f'{k}: real phase convention', u, mp.re(u))
        h = (1+1j*b)/2
        hi = [(1+1j*x)/2 for x in e]
        nsh = period(a,q,-b)*period(c,p,b)
        check(f'{k}: shadow coefficient ratio', ns/nsh,
              shadow_coefficient(1-h,hi[2],hi[3])/
              shadow_coefficient(1-h,hi[0],hi[1]))
        fs, fh = partial_wave(b,e), partial_wave(-b,e)
        for j, z in enumerate([mp.mpc('.2','.37'),mp.mpc('1.2','.45')]):
            check(f'{k}: reflected normalized wave {j}', fs(z)/ns, fh(z)/nsh)
            check(f'{k}: half-plane reflection {j}', fs(z), fs(mp.conj(z)))
        rows.append(dict(external=e,b=b,d=d,ns=ns,nt=nt,ws=ws,wt=wt,
                         raw=raw,u=u,reverse=reverse))
    paths = ['code/ds3_crossing_dictionary.py',
             'code/validate_ds3_crossing_dictionary.py',
             'code/global_racah.py','code/ds3_local_vertex.py',
             'code/ds3_two_matter.py','code/ds3_group_pairing.py']
    result = dict(dps=mp.mp.dps,python=platform.python_version(),
                  mpmath=mp.__version__,sympy=sympy.__version__,
                  hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                  rows=rows,checks=checks,
                  passed=sum(c['passed'] for c in checks),total=len(checks))
    out=ROOT/'results/stage4/ds3-crossing-dictionary/diagnostics.json'
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(encode(result),indent=2),encoding='utf-8')
    print(f"{result['passed']}/{result['total']} checks passed")
    for ch in checks:
        if not ch['passed']:print(encode(ch))
    if result['passed']!=result['total']:raise SystemExit(1)


if __name__=='__main__':run()
