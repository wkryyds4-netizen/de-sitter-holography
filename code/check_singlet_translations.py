"""Exact differential identities and independent high-precision quadratures."""
from pathlib import Path
import hashlib,json,platform,time
import sympy as sp
import mpmath as mp
from singlet_translations import kernel,integral,two_mass_density
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-translations'

def symbolic():
    p,d,m,mu=sp.symbols('p d m mu',positive=True);f=sp.Function('f')(p)
    Q=lambda z:p*z
    B=lambda z:sp.I*(p*sp.diff(z,p)+z/2)
    S=lambda z:-p*sp.diff(z,p,2)-sp.diff(z,p)+(d-sp.Rational(1,2))**2*z/p
    E=lambda z:(Q(z)+S(z))/2
    P=lambda z:(Q(z)-S(z))/2
    pp=lambda z:mu*p*z
    pm=lambda z:m*m*z/(mu*p)
    comm=lambda a,b:a(b(f))-b(a(f))
    exprs={
      'sl2_BE':comm(B,E)-sp.I*P(f),'sl2_BP':comm(B,P)-sp.I*E(f),
      'sl2_EP':comm(E,P)-sp.I*B(f),
      'sl2_Casimir':E(E(f))-P(P(f))-B(B(f))-d*(d-1)*f,
      'B_Pplus':comm(B,pp)-sp.I*pp(f),'B_Pminus':comm(B,pm)+sp.I*pm(f),
      'translations_commute':comm(pp,pm),'mass_shell':pp(pm(f))-m*m*f,
      'lowest_weight':E(p**(d-sp.Rational(1,2))*sp.exp(-p))-d*p**(d-sp.Rational(1,2))*sp.exp(-p)}
    checks=[]
    for name,z in exprs.items():
        value=sp.simplify(sp.expand(z));assert value==0,name
        checks.append(dict(name=name,residual=str(value)))
    # Exact seed action of Lin--Stanford J_ij at any fixed positive h.
    q=sp.symbols('q',positive=True);r=q*q;h=sp.symbols('h',positive=True)
    gram=sp.Matrix([[1,r],[r,1]]);v=sp.Matrix([1,-1])/(2*sp.sqrt(h))
    var=(v.T*gram*v)[0]
    assert sp.simplify(var-(1-q*q)/(2*h))==0
    seed=dict(d=2,Jij=-(1+q),E=(1+q)/q,B=0,P=0)
    symbolic_seed={k:str(v) for k,v in seed.items()}
    symbolic_seed.update(boost_liouvillian_variance=str(var),variance_limit=str(sp.limit(var.subs(q,sp.exp(-h)),h,0)),
        normalized_sl2_relative_defect='1 exactly',chord_Casimir='q**2 + q**(-1)',
        rescaled_Casimir_limit=str(sp.limit((sp.exp(-2*h)+sp.exp(h)-1-sp.exp(-h))/h**2,h,0)))
    # Two independent gamma momenta give x~Beta(4,4), s=1/[x(1-x)].
    x=sp.symbols('x',positive=True)
    density=140*x**3*(1-x)**3
    norm=sp.integrate(density,(x,0,1));mean=sp.integrate(density/(x*(1-x)),(x,0,1))
    assert norm==1 and mean==sp.Rational(14,3)
    return dict(operator_checks=checks,finite_q_seed=symbolic_seed,two_particle=dict(beta_normalization=str(norm),mean_mass_squared_over_m2=str(mean)))

def main():
    target=OUT/'study.json'
    if target.exists():raise FileExistsError(target)
    OUT.mkdir(parents=True,exist_ok=True);start=time.monotonic();exact=symbolic();rows=[]
    specs=[('0.5','0.5','.2','-.1','0','0'),('1','1','.2','-.1','.5','0'),
        ('2','2','.2','-.1','-.5','0'),('1','1','-.3','.1','.25','0'),
        ('1','1','.2','.2','.5','0'),('1','1','.2','-.2','.5','0'),
        ('1','1','.3','-.1','.5','.2'),('2','2','.3','-.1','0','.2')]
    for spec in specs:
        m,mu,a,b,eta,eps=spec
        mp.mp.dps=45;v,tail,amp=integral(m,mu,a,b,eta,eps=eps);ref=kernel(m,mu,a,b,eta,eps=eps)
        err=abs(v-ref);assert err<tail+mp.mpf('1e-30'),(spec,err,tail)
        mp.mp.dps=65;fine,tail2,amp2=integral(m,mu,a,b,eta,eps=eps)
        change=abs(fine-v);assert change<mp.mpf('1e-30'),(spec,change)
        rows.append(dict(parameters=dict(zip(['m','mu','a','b','eta','eps'],spec)),value=str(ref),
            integral=str(v),absolute_difference=str(err),relative_difference=str(err/abs(ref)),
            analytic_tail_bound=str(tail),precision_change=str(change),cancellation_upper_bound=str(amp/abs(v))))
        print('kernel checked',spec,flush=True)
    boost=[]
    for eta in ['0','.5','1','2']:
        ref=mp.sech(mp.mpf(eta)/2)**4
        value=kernel(1,1,0,0,eta)
        assert abs(value-ref)<mp.mpf('1e-55')
        boost.append(dict(eta=eta,value=str(value)))
    # Normalize the invariant-mass density with s=4/(1-y^2), removing the threshold root.
    mass_norm=mp.quad(lambda y:(1-y*y)**3/(mp.beta(4,4)*64),[0,1])
    mass_mean=mp.quad(lambda y:4*(1-y*y)**2/(mp.beta(4,4)*64),[0,1])
    assert abs(mass_norm-1)<mp.mpf('1e-55') and abs(mass_mean-mp.mpf(14)/3)<mp.mpf('1e-55')
    masses=[dict(m=str(m),spacetime_value=str(kernel(m,m,'.5',0,0))) for m in [mp.mpf('.5'),mp.mpf(1),mp.mpf(2)]]
    finite=[dict(h=str(h),liouvillian_variance=str(-mp.expm1(-2*h)/(2*h)),length_boost_variance='0') for h in map(mp.mpf,['.1','.01','.001'])]
    names=['code/singlet_translations.py','code/check_singlet_translations.py']
    data=dict(status='passed',scope='Exact candidate representation identities, finite-q seed diagnosis, and high-precision kernel checks; no microscopic spatial dictionary claimed',
        symbolic=exact,kernels=rows,boost=boost,two_particle_density=dict(normalization=str(mass_norm),mean=str(mass_mean)),
        mass_comparison=masses,finite_q_comparison=finite,
        error_scope='Explicit absolute tail bounds; interior quadrature checked at 45/65 digits but not interval-certified. No finite-h dynamical approximation is used in candidate kernel.',
        environment=dict(python=platform.python_version(),sympy=sp.__version__,mpmath=mp.__version__),
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in names},seconds=time.monotonic()-start)
    target.write_text(json.dumps(data,indent=2)+'\n');print(json.dumps(dict(status='passed',seconds=data['seconds'],mass_comparison=masses)))
if __name__=='__main__':main()
