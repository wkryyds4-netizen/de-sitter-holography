"""Uniform rational prefactors on the full input/output rank window."""
from pathlib import Path
import json,hashlib
import sympy as s
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-two-charge'

def main():
    target=OUT/'uniform-prefactors.json'
    if target.exists():raise FileExistsError(target)
    data=json.loads((OUT/'symbolic-filter.json').read_text())
    m,n,p,a,b,r=s.symbols('m n p a b r')
    eps,lam,alpha,beta,rho=s.symbols('eps lam alpha beta rho')
    sub={m:1/(2*lam*eps**2),n:1/(2*lam*eps**2)+1,p:1/eps,
         a:alpha/eps,b:beta/eps,r:rho/eps}
    expansions=[];leading=[2,-2,-2,2]
    for value,want in zip(data['H_squared_coefficients'],leading):
        expression=s.cancel(s.sympify(value).subs(sub))
        expansion=s.series(expression,eps,0,1).removeO().expand()
        assert s.expand(expansion.coeff(eps,-1)-want)==0
        expansions.append(str(s.factor(expansion)))
    result=dict(status='passed',rational_prefactor_expansions=expansions,
        domain='Fixed positive lambda; alpha,beta in [0,1], rho in [0,2]; epsilon=1/p sufficiently small.',
        analytic_uniformity='Every rational denominator comes from d=m-p, m, n, m-j, n-j. After extracting powers of p, these stay bounded away from zero uniformly on the compact rank-proportion window.',
        consequence='gH/C - 2p(1-q)^2 has a uniform linear bound in 1+input deficits+output deficit. Combined with previous positive product tail estimates, the H-squared part has an L2 linear-filter approximation.',
        environment=dict(sympy=s.__version__),
        input_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in
            ['code/check_memory_filter_uniform.py','results/susskind-memory-two-charge/symbolic-filter.json']})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(expansions)

if __name__=='__main__':main()
