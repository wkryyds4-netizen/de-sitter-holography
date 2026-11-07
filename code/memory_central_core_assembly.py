"""Uniform-rank rational coefficients and central core assembly algebra."""
from pathlib import Path
import hashlib,json,platform
import sympy as s
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-central-norm'


def main():
    target=OUT/'core-assembly.json'
    if target.exists():raise FileExistsError(target)
    source='results/susskind-memory-two-charge/symbolic-filter.json'
    old=json.loads((ROOT/source).read_text())
    m,n,p,a,b,r=s.symbols('m n p a b r')
    e,lam,x,alpha,beta,rho,q=s.symbols('e lam x alpha beta rho q')
    sub={m:1/(2*lam*e**2)+x/e,n:1/(2*lam*e**2)-x/e+1,
         p:1/e,a:alpha/e,b:beta/e,r:rho/e}
    expected=[2/e+lam*(1+4*x),-2/e+lam*(6-11*alpha-4*x),
              -2/e+lam*(6-11*beta-4*x),2/e+lam*(5*alpha+5*beta+6*rho-13+4*x)]
    rows=[]
    for val,want in zip(old['H_squared_coefficients'],expected):
        ex=s.series(s.cancel(s.sympify(val).subs(sub)),e,0,1).removeO()
        assert s.factor(ex-want)==0
        rows.append(str(s.factor(ex)))
    K0=lam*(1-10*q+9*q*q+16*lam*q-32*lam*q*q)
    centered=s.factor(K0+16*lam*lam*(q*q-q)+4*lam*x*(1-q)**2)
    assert s.factor(centered-(lam*(1-10*q+9*q*q)-16*lam*lam*q*q+4*lam*x*(1-q)**2))==0
    bb=1+q;bo=2*lam*(1-q);aa=4*lam+4*lam**2*(1+q);ll=2*lam*(1+q)
    gram=s.Matrix([[bb,0,bo],[0,ll,0],[bo,0,aa]])
    y=s.factor(aa-bo*bo/bb)
    assert s.factor(y-(4*lam+16*lam**2*q/(1+q)))==0
    assert s.factor(gram.det()-bb*ll*y)==0
    names=['code/memory_central_core_assembly.py',source]
    target.write_text(json.dumps(dict(status='passed',uniform_H_prefactors=rows,
        centered_core_B_finite_coefficient=str(centered),
        gram_B_Linput_O=[[str(v) for v in row] for row in gram.tolist()],
        Y_variance_over_C_squared=str(y),gram_determinant=str(s.factor(gram.det())),
        analytic_domain='Bounded x and positive compact lambda, a,b<=p,r<=a+b; scaled rational denominators uniformly nonzero.',
        strong_A_filter_bound='||sum(gA-C*gamma)W2A_ab;r||^2 <= K*C^4*(p^-2+R^2/p^2+p^4/R^6); R=p^(3/4).',
        scope='Coefficient identities and Gram algebra; proof of strong replacement uses general coherent norm, factorial input tail and central sixth output moment.',
        environment=dict(python=platform.python_version(),sympy=s.__version__),
        input_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names}),indent=2)+'\n')
    print('central uniform H coefficients and core assembly passed')


if __name__=='__main__':main()
