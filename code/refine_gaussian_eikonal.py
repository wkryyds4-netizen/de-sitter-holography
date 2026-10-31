"""Higher-order fresh quadratures and exact symbolic bookkeeping."""
import hashlib,json
from pathlib import Path
import numpy as np
import sympy as sp
from gaussian_eikonal import smeared

def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/gaussian-eikonal';dest=base/'refinement.json'
    if dest.exists():raise RuntimeError('Preserve refinement')
    data=json.loads((base/'study.json').read_text());rows=[];checks=[]
    for row in data['rows']:
        new=smeared(row['m'],row['gamma'],row['t'],row['sigma'],48,64)
        diff=abs(complex(*new['value'])-complex(*row['fine']['value']))
        rows.append(dict(gamma=row['gamma'],sigma=row['sigma'],u=row['u'],value=new,difference=diff))
        checks.append(dict(label='32/48 quadrature',passed=bool(diff<2e-7),difference=diff,tolerance=2e-7))
        print(row['gamma'],row['sigma'],row['u'],diff,flush=True)
    mat=sp.Matrix([[sp.Rational(1,2),0,-sp.Rational(1,2),0],[0,sp.Rational(1,2),0,-sp.Rational(1,2)],
                   [sp.Rational(1,2),-sp.Rational(1,2),sp.Rational(1,2),-sp.Rational(1,2)]])
    checks.append(dict(label='Gaussian covariance and independence',passed=mat*mat.T==sp.diag(sp.Rational(1,2),sp.Rational(1,2),1)))
    z,n,r=sp.symbols('z n r')
    A=lambda k:1-k*z/4+k*(3*k+2)*z*z/32
    logq=sp.series(n*n*z/2+2*(sp.log(A(n+1))-sp.log(A(1))),z,0,3).removeO().expand()
    checks.append(dict(label='small-width coefficient',passed=sp.simplify(logq-n*(n-1)*z/2-n*(n+3)*z*z/8)==0))
    diagnostics=json.loads((base/'diagnostics.json').read_text())
    checks.append(dict(label='52 independent diagnostics retained',passed=diagnostics['status']=='passed' and len(diagnostics['checks'])==52 and all(x['passed'] for x in diagnostics['checks'])))
    out=dict(status='passed' if all(x['passed'] for x in checks) else 'failed',rows=rows,checks=checks,
             scope='Higher quadrature order and exact covariance/small-width algebra. No interval or omitted-eikonal error certification.',
             hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in ['code/refine_gaussian_eikonal.py','code/gaussian_eikonal.py','results/stage4/gaussian-eikonal/study.json','results/stage4/gaussian-eikonal/diagnostics.json']})
    dest.write_text(json.dumps(out,indent=2),encoding='utf-8');print(out['status'],len(checks))
    if out['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
