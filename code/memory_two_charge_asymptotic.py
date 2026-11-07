"""Symbolic cancellation and fixed-deficit expansion of the two-charge filter."""
from pathlib import Path
import hashlib,json,platform
import sympy as s
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-two-charge'

def main():
    target=OUT/'symbolic-filter.json'
    if target.exists():raise FileExistsError(target)
    m,n,p,a,b,r=s.symbols('m n p a b r');C,Cm,jc,ka,kb,kr,ea,eb,er=s.symbols('C Cm jc ka kb kr ea eb er')
    d=m-p;za=(m-a)*(n-a);zb=(m-b)*(n-b)
    kapminus=lambda k,e:k+(m-n-p)*e
    va,vb,vr=ka-n*ea,kb-n*eb,kr-n*er
    ua,ub=d*(jc+va),d*(jc+vb)
    D4=(0,3*Cm+kapminus(ka,ea)+kapminus(kb,eb)+kapminus(kr,er))
    D13=(d*(4*Cm+kapminus(ka,ea)+kapminus(kb,eb)-d*(ea+eb)),
        2*(jc+vr)+ka+kb+(m-n)*(ea+eb))
    D22=(2*m*Cm+ua+ub+d*(ka+kb)+2*d*kr,
        2*C-p*(ka+kb)/d-er*(za+zb)/d)
    D121=(ua+ub+m*Cm,kr);D40=(3*C+ka+kb+kr,0)
    raw=[D4[i]-2*D13[i]+s.Rational(5,4)*D22[i]+s.Rational(7,2)*D121[i]+(4*p-3*m)*D40[i] for i in range(2)]
    gp=(8*p-2*m)*Cm+s.Rational(19,2)*d*jc+(12*p-9*m)*C+m*(ka+kb)+(3*p-m)*kr/2-s.Rational(11,4)*d*n*(ea+eb)
    gm=3*Cm+s.Rational(5,2)*C-4*jc-(4*m+p)*(ka+kb)/(4*d)+kr/2-(m-n+p)*(ea+eb)+(m+3*n-p-5*(za+zb)/(4*d))*er
    assert all(s.cancel(x-y)==0 for x,y in zip(raw,[gp,gm]))
    h=s.cancel(gp+(m-2*p)*gm)
    h=h.subs({Cm:C*d*(n+p)/(m*n),jc:C*d/m,ea:ka*p/za,eb:kb*p/zb,er:kr*p/((m-r)*(n-r))})
    hcoeff=[s.cancel(s.diff(h,x)) for x in [C,ka,kb,kr]]
    eps,lam,u,v,S,q=s.symbols('eps lam u v S q',positive=True)
    sub={m:1/(2*lam*eps**2),n:1/(2*lam*eps**2)+1,p:1/eps,a:1/eps-u,b:1/eps-v,r:2/eps-u-v-S}
    expansions=[]
    for value in hcoeff:
        ex=s.series(s.cancel(value.subs(sub)),eps,0,1).removeO().expand()
        expansions.append(s.factor(ex));print('H coefficient',ex,flush=True)
    kexp=[1,q*(1+eps*(4*lam*u-4*lam**2)),q*(1+eps*(4*lam*v-4*lam**2)),q**2*(1+eps*(4*lam*(u+v+S)-16*lam**2))]
    combined=s.expand(sum(x*y for x,y in zip(expansions,kexp)))
    lead=s.factor(combined.coeff(eps,-1));finite=s.factor(combined.coeff(eps,0))
    assert s.expand(lead-2*(1-q)**2)==0
    assert s.diff(finite,u)==s.diff(finite,v)
    assert s.factor(s.diff(finite,u)-8*lam*(q*q-q))==0
    assert s.factor(s.diff(finite,S)-8*lam*q*q)==0
    gmnorm=gm.subs({Cm:C*d*(n+p)/(m*n),jc:C*d/m,ea:ka*p/za,eb:kb*p/zb,er:kr*p/((m-r)*(n-r))})/C
    gmnorm=gmnorm.subs({ka:C*q,kb:C*q,kr:C*q*q})
    glimit=s.factor(s.limit(s.cancel(gmnorm.subs(sub)),eps,0))
    assert s.factor(glimit-(3-4*q+q*q)/2)==0
    names=['code/memory_two_charge_asymptotic.py','code/memory_two_charge_filter.py']
    result=dict(status='passed',simplified_g_plus=str(gp),simplified_g_minus=str(gm),
        H_squared_coefficients=list(map(str,hcoeff)),half_filling_rational_coefficient_expansions=list(map(str,expansions)),
        gA_over_C_limit=str(glimit),gH_over_C_leading=str(lead),gH_over_C_finite=str(finite),
        constant_finite_part=str(s.factor(finite.subs({u:0,v:0,S:0}))),
        scope='Fixed deficits at exact half filling; coefficient asymptotics only, not an annealed operator limit.',
        environment=dict(python=platform.python_version(),sympy=s.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print('leading',lead,'finite',finite,'gA',glimit,flush=True)

if __name__=='__main__':main()
