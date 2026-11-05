"""Read-only check of saved calibration against charge-sector identities."""
from pathlib import Path
from fractions import Fraction
from math import comb, factorial
import hashlib
import json
import numpy as np

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-dictionary'


def audit(root=ROOT):
    dest=root/'results/susskind-dictionary'
    r=json.loads((dest/'calibration.json').read_text(encoding='utf-8'))
    errors=[]
    rows=[]
    for name,sha in r['code_sha256'].items():
        if hashlib.sha256((root/name).read_bytes()).hexdigest()!=sha:
            errors.append('code hash: '+name)
    if not r['passed'] or not all(c['passed'] and c['error']<=c['tolerance'] for c in r['checks']):
        errors.append('saved diagnostics')
    for row in r['covariance_rows']:
        n,p=row['N'],row['p']
        v=4**p*Fraction(n,p*p)*Fraction(factorial(p),n**p)**2
        f=[v*comb(m,p)*comb(n-m+p,p) if m>=p else Fraction(0) for m in range(n+1)]
        w=sum(Fraction(comb(n,m),2**n)*f[m] for m in range(n+1))
        a=sum(Fraction(comb(n,m),2**n)*f[m]*Fraction(p*(n-p+1),n-m+1) for m in range(n+1))
        err=max(abs(float(w)-row['H2_formula']),abs(float(a)-row['A2_formula']))
        if err>1e-12:
            errors.append(f'sector normalization {n}-{p}')
        rows.append(dict(N=n,p=p,H2_sector=float(w),A2_sector=float(a),error=err))
    for s in r['samples']:
        path=dest/s['file']
        if hashlib.sha256(path.read_bytes()).hexdigest()!=s['sha256']:
            errors.append('sample hash: '+s['file'])
        with np.load(path,allow_pickle=False) as data:
            h,a=data['H'],data['A2']
            w,om=data['weight'],data['omega']
            d=len(h)
            o=a-np.trace(a)/d*np.eye(d)
            ho=h@o-o@h
            direct_g=float(np.trace(o@o).real/d)
            direct_k=float(np.trace(ho.conj().T@ho).real/d)
            err=max(abs(w.sum()-direct_g)/max(1,direct_g),
                    abs(np.dot(w,om**2)-direct_k)/max(1,direct_k))
            if err>2e-11 or np.any(w<0):
                errors.append('spectral moments: '+s['file'])
    return dict(status='passed' if not errors else 'failed',errors=errors,
                sample_hashes_checked=len(r['samples']),code_hashes_checked=len(r['code_sha256']),
                sector_rows=rows,scope='Saved finite-system moments and independent sector normalization; no large-N singlet claim.')


if __name__=='__main__':
    result=audit()
    with (OUT/'audit.json').open('x',encoding='utf-8') as f:
        json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result))
    if result['errors']:
        raise SystemExit(1)
