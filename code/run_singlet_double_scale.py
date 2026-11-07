"""One-shot integer coefficients, graph table and independent checks."""
from pathlib import Path
from fractions import Fraction
from math import comb,factorial,exp
import hashlib,json,platform,time
import numpy as np
import sympy as sp
from singlet_coefficient import moment_coefficients,evaluate,overlap_table
from check_singlet_coefficient import check_words
from singlet_wick_moments import compute

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-singlet-double-scale'


def save(name,obj):
    with (OUT/name).open('x',encoding='utf-8') as f:json.dump(obj,f,indent=2,default=str);f.write('\n')


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    if any(OUT.glob('*.json')):raise FileExistsError('Refusing to replace existing stage data')
    checks=check_words();save('word-checks.json',checks)
    print(f"Direct CAR: {checks['scalar_comparisons']} comparisons passed",flush=True)
    rows=[];n=sp.symbols('N')
    for p in range(1,6):
        start=time.perf_counter();cc=moment_coefficients(p)
        if p<=3:assert cc==moment_coefficients(p,'all15')
        assert all(evaluate(nn,p,cc)==0 for nn in range(p+2))
        if p==1:assert not any(cc)
        expr=sp.factor(sum(sp.prod(n-i for i in range(j))*v/sp.factorial(j)
                           for j,v in enumerate(cc))/2**(6*p))
        samples=[]
        for lam in [Fraction(1,4),Fraction(1,2),Fraction(1)]:
            nn=int(p*p/lam)
            raw=evaluate(nn,p,cc)
            scaled=raw*2**(6*p)*Fraction(factorial(p),nn**p)**6
            q=exp(-4*float(lam))
            prediction=8*float(lam)*(1-q*q)+32*float(lam)**2*(q-q*q)
            samples.append(dict(N=nn,lambda_value=lam,normalized_mu2=scaled,
                                double_scaled_prediction=prediction,ratio=float(scaled)/prediction))
        row=dict(p=p,coefficients=cc,factored_raw_mu2=str(expr),samples=samples,
                 elapsed_seconds=time.perf_counter()-start)
        rows.append(row);save(f'coefficients-p{p}.json',row)
        print(f'p={p} complete, {row["elapsed_seconds"]:.2f} s',flush=True)
    # A different N from the preceding stage, evaluated by full Gaussian tensors.
    independent,_=compute(5,2)
    assert independent['total']['mu2']==evaluate(5,2,rows[1]['coefficients'])==2700
    save('independent-tensor.json',independent)
    graph=overlap_table();grouped={c:[0,0] for c in range(4)}
    for r in graph:
        grouped[r['crossings']][0]+=r['sign']*r['delta']
        grouped[r['crossings']][1]+=r['sign']*r['alpha']*r['beta']
    assert grouped=={0:[8,0],1:[0,32],2:[-8,-32],3:[0,0]}
    save('overlap-graphs.json',dict(rows=graph,lambda_and_lambda_squared_coefficients=grouped))
    sources=['code/singlet_coefficient.py','code/check_singlet_coefficient.py',
             'code/run_singlet_double_scale.py','code/singlet_wick_moments.py','code/complex_syk_dictionary.py',
             'results/susskind-singlet-moments/moments.json']
    save('validation.json',dict(status='passed',environment=dict(python=platform.python_version(),
         numpy=np.__version__,sympy=sp.__version__),
         input_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sources},
         exact_claims=['720 marked word cases, 2880 exact scalar comparisons',
           'p=1,2,3 all-15 and cross-only coefficient arrays agree',
           'p=1 identically zero; N<=p+1 controls vanish for p=1..5',
           'N=5,p=2 independent full covariance tensor gives 2700',
           '24 marked overlap graphs produce the claimed double-scaled expression'],
         asymptotic_scope='Finite-p ratios are checks of approach only; limit derived analytically from Poisson overlaps.'))
    print('All stage checks passed',flush=True)


if __name__=='__main__':main()
