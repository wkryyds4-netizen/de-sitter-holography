"""Exact endpoint certificate for Q4 versus the reverse fermionic force.

The old forward certificate is reused as an algorithm, with only the
independently verified opposite-orientation local minor changed. All
nonzero polynomials are retained; no zero result is presumed.
"""
from pathlib import Path
import hashlib, inspect, json, platform
import sympy as s
import quartic_cross_kernel as old
from quartic_orbits import q

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-oriented-kernel'


def main():
    target=OUT/'reverse-cross.json'
    if target.exists():raise FileExistsError(target)
    source=inspect.getsource(old.cross_certificate)
    before='4*int(pi[A[a]]==B[e] and pi[A[b]]==B[c])'
    after='-4*int(pi[A[a]]==B[c] and pi[A[b]]==B[e])'
    assert source.count(before)==1
    modified=source.replace(before,after)
    env=vars(old).copy()
    exec(compile(modified,'<reverse endpoint certificate>','exec'),env)
    result=env['cross_certificate'](2)
    # Symbolic projection and positive spectral matrix; C denotes sigma^2.
    l,C=s.symbols('lam C',positive=True)
    z=1+q;c=2*l*(1-q)/z
    w=4*C**2*l+16*C**2*l**2*q/z
    x=-4*C**2*l*q+16*C**2*l**2*q/z
    v=2*l*(1-q)+16*l*l*q/z;d=2*l*z
    rho=v/(v+d)
    checks=[w+x-2*C**2*v,w-x-2*C**2*d,
            (s.Rational(1,2)-rho)*w+x/2,
            (s.Rational(1,2)-rho)**2*w+w/4+(s.Rational(1,2)-rho)*x-C**2*v*d/(v+d)]
    assert all(s.factor(t)==0 for t in checks)
    vnew=C**2*l*z*(1-q*q+8*l*q)/(z+4*l*q)
    assert s.factor(vnew-C**2*v*d/(v+d))==0
    v2=(8*l*C*q*q)**2*vnew
    small=s.limit(v2.subs({C:1/l,q:s.exp(-4*l)}),l,0,dir='+')
    assert small==1024
    prior=json.loads((ROOT/'results/susskind-memory-quartic/quartic.json').read_text())
    vv=s.sympify(prior['residual_variance'].replace('lambda','lam'),locals={'q':q,'lam':l})
    v4=C**4*vv
    ratio=s.limit((v2/v4/l).subs({C:1/l,q:s.exp(-4*l)}),l,0,dir='+')
    assert ratio==s.Rational(20,3)
    names=['code/memory_oriented_cross.py','code/quartic_cross_kernel.py',
           'code/quartic_orbits.py','code/memory_oriented_marks.py',
           'results/susskind-memory-oriented-kernel/oriented-marks.json',
           'results/susskind-memory-quartic/quartic.json']
    data=dict(status='passed' if not result['counterexamples'] else 'counterexamples',
              certificate=result,algorithm_change=dict(before=before,after=after),
              generated_function_sha256=hashlib.sha256(modified.encode()).hexdigest(),
              spectral_matrix=dict(diagonal=str(w),off_diagonal=str(x)),
              symbolic_identities=len(checks)+1,Vnew=str(vnew),V2=str(v2),V4=str(v4),
              small_lambda_V2_Js8=str(small),small_lambda_V2_over_lambda_V4=str(ratio),
              environment=dict(python=platform.python_version(),sympy=s.__version__),
              input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],result['templates'],'templates',result['words'],'words',flush=True)


if __name__=='__main__':main()
