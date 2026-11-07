"""Companion equal-time Gaussian-centered Gram, using two cross pairings."""
from pathlib import Path
from fractions import Fraction
from math import comb,factorial,exp
import json,hashlib
from singlet_coefficient import local_polynomial,power_coefficients

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-singlet-double-scale'


def centered_coefficients(p):
    cc=[0]*(4*p+1)
    for pair in [((0,2),(1,3)),((0,3),(1,2))]:
        row=power_coefficients(p,4,local_polynomial(4,pair,((0,1),(2,3))))
        cc=[a+b for a,b in zip(cc,row)]
    return cc


def main():
    target=OUT/'centered-gram.json'
    if target.exists():raise FileExistsError(target)
    old=json.loads((ROOT/'results/susskind-singlet-moments/moments.json').read_text())
    rows=[]
    for p in range(1,6):
        cc=centered_coefficients(p)
        def evaluate(n):return Fraction(sum(comb(n,j)*v for j,v in enumerate(cc) if j<=n),2**(4*p))
        if p==1:
            for n in range(1,12):assert evaluate(n)==Fraction(n*(n+1)**2,4)
        if p==2:assert evaluate(4)==Fraction(old['records'][-1]['total']['ensemble_centered_gram'])
        samples=[]
        dynamics=json.loads((OUT/f'coefficients-p{p}.json').read_text())['samples']
        for dy in dynamics:
            n=dy['N'];lam=Fraction(dy['lambda_value'])
            normalized=evaluate(n)*2**(4*p)*Fraction(factorial(p),n**p)**4
            pred=4*float(lam)+4*float(lam)**2*(1+exp(-4*float(lam)))
            curvature=Fraction(dy['normalized_mu2'])/normalized/lam
            samples.append(dict(N=n,lambda_value=lam,normalized_centered_gram=normalized,
                prediction=pred,ratio=float(normalized)/pred,curvature_over_Js_squared=curvature))
        rows.append(dict(p=p,coefficients=cc,samples=samples))
    payload=dict(status='passed',definition='E tr(A-E[A]|Q)^2, retaining any sample-dependent conserved part',
        double_scaled='4 Js^4 (1/lambda+1+exp(-4lambda))',rows=rows,
        input_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in
            ['code/run_singlet_centered_moment.py','code/singlet_coefficient.py',
             'results/susskind-singlet-moments/moments.json']})
    with target.open('x') as f:json.dump(payload,f,indent=2,default=str);f.write('\n')
    print('Centered Gram recorded; exact p=1 and N=4,p=2 controls passed')


if __name__=='__main__':main()
