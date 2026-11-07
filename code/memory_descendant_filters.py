"""Exact lower-chaos filters for the five raw quartic single-mark words."""
from pathlib import Path
import hashlib,json,platform
import sympy as s

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-descendants'


def symbolic_words():
    m,n,p,a,b,r=s.symbols('m n p a b r')
    C,Cm,jc,ka,kb,kr,ea,eb,er=s.symbols('C Cm jc ka kb kr ea eb er')
    d=m-p;za=(m-a)*(n-a);zb=(m-b)*(n-b)
    km=lambda k,e:k+(m-n-p)*e
    ua,ub=d*(jc+ka-n*ea),d*(jc+kb-n*eb)
    words=s.Matrix([
        [0,3*Cm+km(ka,ea)+km(kb,eb)+km(kr,er)],
        [d*(4*Cm+km(ka,ea)+km(kb,eb)-d*(ea+eb)),
         2*(jc+kr-n*er)+ka+kb+(m-n)*(ea+eb)],
        [2*m*Cm+ua+ub+d*(ka+kb)+2*d*kr,
         2*C-p*(ka+kb)/d-er*(za+zb)/d],
        [ua+ub+m*Cm,kr],
        [3*C+ka+kb+kr,0]])
    # Rows E0,M0,S0,C0,Z0; columns R4,{h,R3},{h2,R2},hR2h,h4.
    raw=s.Matrix([[0,0,s.Rational(1,2),0,2*p-m],
                  [0,0,0,1,2*p-m],
                  [0,s.Rational(1,2),-s.Rational(1,2),-1,m-p],
                  [1,-1,0,1,0],
                  [0,0,0,0,1]])
    filters=raw*words
    substitution={Cm:C*d*(n+p)/(m*n),jc:C*d/m,
                  ea:ka*p/za,eb:kb*p/zb,er:kr*p/((m-r)*(n-r))}
    pair=[]
    for gp,gm in filters.tolist():
        pair.append([s.cancel((gp+(m-2*p)*gm).subs(substitution)),s.cancel(gm.subs(substitution))])
    return (m,n,p,a,b,r,C,ka,kb,kr),raw,pair


def main():
    OUT.mkdir(exist_ok=True);target=OUT/'symbolic-filters.json'
    if target.exists():raise FileExistsError(target)
    vars,raw,pairs=symbolic_words();m,n,p,a,b,r,C,ka,kb,kr=vars
    e,lam,u,v,S,q=s.symbols('e lam u v S q')
    alpha,beta,rho=s.symbols('alpha beta rho')
    fixed={m:1/(2*lam*e**2),n:1/(2*lam*e**2)+1,p:1/e,a:1/e-u,b:1/e-v,r:2/e-u-v-S}
    broad={m:1/(2*lam*e**2),n:1/(2*lam*e**2)+1,p:1/e,a:alpha/e,b:beta/e,r:rho/e}
    kap=[1,q*(1+e*(4*lam*u-4*lam**2)),q*(1+e*(4*lam*v-4*lam**2)),
         q*q*(1+e*(4*lam*(u+v+S)-16*lam**2))]
    rows=[]
    for name,(hh,aa) in zip(['E','M','S','C','Z'],pairs):
        coeff=[[s.cancel(s.diff(expr,z)) for z in [C,ka,kb,kr]] for expr in [hh,aa]]
        limits=[];uniform=[];fixed_coeff=[]
        for cc in coeff:
            ff=[s.series(s.cancel(c.subs(fixed)),e,0,1).removeO().expand() for c in cc]
            uu=[s.series(s.cancel(c.subs(broad)),e,0,1).removeO().expand() for c in cc]
            assert all(not any(c.coeff(e,i)!=0 for i in range(-5,-1)) for c in ff+uu)
            limits.append(s.expand(sum(c*k for c,k in zip(ff,kap))))
            uniform.append(list(map(str,uu)));fixed_coeff.append(list(map(str,ff)))
        Hlead=s.factor(limits[0].coeff(e,-1));Hfinite=s.factor(limits[0].coeff(e,0))
        Alead=s.factor(limits[1].coeff(e,-1));Alimit=s.factor(limits[1].coeff(e,0))
        assert Alead==0
        input_coeff=s.factor(s.diff(Hfinite,u));output_coeff=s.factor(s.diff(Hfinite,S))
        assert s.factor(s.diff(Hfinite,v)-input_coeff)==0
        constant=s.factor(Hfinite.subs({u:0,v:0,S:0})+2*lam*input_coeff)
        assert s.factor(Hfinite-Hfinite.subs({u:0,v:0,S:0})-input_coeff*(u+v)-output_coeff*S)==0
        bridge=s.factor(input_coeff/(4*lam*q))
        row=dict(name=name,H_prefactors=list(map(str,coeff[0])),A_prefactors=list(map(str,coeff[1])),
                 fixed_prefactors=fixed_coeff,uniform_prefactors=uniform,
                 B_leading=str(Hlead),B_finite_centered=str(constant),O_coefficient=str(Alimit),
                 JhF_coefficient=str(bridge),S_coefficient=str(output_coeff))
        rows.append(row);print(name,[row[k] for k in ['B_leading','B_finite_centered','O_coefficient','JhF_coefficient','S_coefficient']],flush=True)
    old=json.loads((ROOT/'results/susskind-memory-quartic/quartic.json').read_text())
    # Match symbols in the saved fourth-chaos projection.
    projection=[s.sympify(x.replace('lambda','lam'),locals={'lam':lam,'q':q}) for x in old['projection']]
    vec=[-projection[0],-projection[1],-projection[2],s.Integer(1),-projection[3]]
    combined={k:str(s.factor(sum(c*s.sympify(row[k]) for c,row in zip(vec,rows))))
              for k in ['B_leading','B_finite_centered','O_coefficient','JhF_coefficient','S_coefficient']}
    core=[s.Rational(1,2),s.Rational(1,2),-2,1,0]
    assert s.factor(sum(c*s.sympify(row['S_coefficient']) for c,row in zip(core,rows))-8*lam*q*q)==0
    names=['code/memory_descendant_filters.py','code/memory_two_charge_filter.py',
           'results/susskind-memory-quartic/quartic.json']
    result=dict(status='passed',raw_word_matrix=[[str(x) for x in row] for row in raw.tolist()],
                rows=rows,raw_Q4_combination=combined,raw_Q4_coefficients=list(map(str,vec)),
                environment=dict(python=platform.python_version(),sympy=s.__version__),
                input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print('Raw Q4 lower chaos',combined,flush=True)


if __name__=='__main__':main()
