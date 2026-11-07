"""Central-charge expansion of exact raw quartic filters and cubic rank slope."""
from pathlib import Path
import hashlib,json,platform,time
import sympy as s
from memory_descendant_filters import symbolic_words

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-central-raw'

def main():
    OUT.mkdir(exist_ok=True);target=OUT/'symbolic-descendants.json'
    if target.exists():raise FileExistsError(target)
    start=time.time();vs,raw,pairs=symbolic_words();m,n,p,a,b,r,C,ka,kb,kr=vs
    e,lam,x,u,v,S,q=s.symbols('e lam x u v S q')
    alpha,beta,rho=s.symbols('alpha beta rho')
    fixed={m:1/(2*lam*e**2)+x/e,n:1/(2*lam*e**2)-x/e+1,p:1/e,
           a:1/e-u,b:1/e-v,r:2/e-u-v-S}
    broad={m:fixed[m],n:fixed[n],p:1/e,a:alpha/e,b:beta/e,r:rho/e}
    kap=[1,q*(1+e*(4*lam*u-4*lam**2)),q*(1+e*(4*lam*v-4*lam**2)),
         q*q*(1+e*(4*lam*(u+v+S)-16*lam**2))]
    old=json.loads((ROOT/'results/susskind-memory-descendants/symbolic-filters.json').read_text())
    rows=[]
    for name,(hh,aa),prior in zip(['E','M','S','C','Z'],pairs,old['rows']):
        coeff=[[s.cancel(s.diff(expr,z)) for z in [C,ka,kb,kr]] for expr in [hh,aa]]
        limits=[];uniform=[]
        for cc in coeff:
            ff=[s.series(s.cancel(c.subs(fixed)),e,0,1).removeO().expand() for c in cc]
            uu=[s.series(s.cancel(c.subs(broad)),e,0,1).removeO().expand() for c in cc]
            assert all(not any(c.coeff(e,i)!=0 for i in range(-5,-1)) for c in ff+uu)
            limits.append(s.expand(sum(c*k for c,k in zip(ff,kap))));uniform.append(list(map(str,uu)))
        lead=s.factor(limits[0].coeff(e,-1));finite=s.factor(limits[0].coeff(e,0));oa=s.factor(limits[1].coeff(e,0))
        assert limits[1].coeff(e,-1)==0
        inp=s.factor(s.diff(finite,u));out=s.factor(s.diff(finite,S))
        assert s.factor(s.diff(finite,v)-inp)==0
        assert s.factor(finite-finite.subs({u:0,v:0,S:0})-inp*(u+v)-out*S)==0
        centered=s.factor(finite.subs({u:0,v:0,S:0})+2*lam*inp)
        row=dict(name=name,B_leading=str(lead),B_finite_centered=str(centered),O_coefficient=str(oa),
                 L_input_coefficient=str(inp),JhF_coefficient=str(s.factor(inp/(4*lam*q))),
                 S_coefficient=str(out),uniform_prefactors=uniform)
        for key in ['B_leading','O_coefficient','JhF_coefficient','S_coefficient']:
            assert s.factor(s.sympify(row[key])-s.sympify(prior[key]))==0,(name,key)
        row['central_B_shift']=str(s.factor(centered-s.sympify(prior['B_finite_centered'])))
        rows.append(row);print(name,'central B shift',row['central_B_shift'],flush=True)
    vec=[s.sympify(z) for z in old['raw_Q4_coefficients']]
    shift=s.factor(sum(c*s.sympify(row['central_B_shift']) for c,row in zip(vec,rows)))
    corevec=[s.Rational(1,2),s.Rational(1,2),-2,1,0]
    coreshift=s.factor(sum(c*s.sympify(row['central_B_shift']) for c,row in zip(corevec,rows)))
    assert s.factor(coreshift-4*lam*x*(1-q)**2)==0
    # Independent rank-dependent part of the exact cubic W1 coefficient.
    pp=1/e;NN=pp*pp/lam;mm=NN/2+x*pp;nu=pp-u
    A=mm-nu;B=NN-mm-nu+1;al=nu*(NN-nu-pp+1)/(NN-2*nu+1)
    term=(mm-pp)*(1-pp/A)*(1+pp/B)-mm+3*(al*(1+pp/B)+(pp-al)*(1-pp/A))
    es=s.series(s.cancel(term),e,0,1).removeO().expand()
    dep=s.expand(q*(1+e*(4*lam*u-4*lam**2))*es)
    slope=s.factor(dep.coeff(e,0)-dep.coeff(e,0).subs(u,0))
    assert s.factor(slope-8*lam*q*u)==0
    R1=u*(NN-2*pp+u+1)*(NN-mm+pp)/((mm-pp+1)*pp*pp)
    ratio=s.factor(s.limit(R1,e,0));assert s.factor(ratio-u/lam)==0
    names=['code/memory_central_descendants.py','code/memory_descendant_filters.py',
           'results/susskind-memory-descendants/symbolic-filters.json','code/cubic_linear_channels.py','code/cubic_linear_limit.py']
    data=dict(status='passed',rows=rows,raw_Q4_B_shift=str(shift),core_B_shift=str(coreshift),
        cubic_rank_coefficient_difference_limit=str(slope),p_minus_1_covariance_ratio_limit=str(ratio),
        scope='Exact rational prefactor expansion at general central charge. All leading dynamic coefficients equal half filling; finite corrections lie in B. Cubic centered rank slope is independent of x. Uniform tail arguments and strong replacements are given in the note, not inferred from these expansions alone.',
        elapsed_seconds=time.time()-start,environment=dict(python=platform.python_version(),sympy=s.__version__),
        input_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names})
    target.write_text(json.dumps(data,indent=2)+'\n')

if __name__=='__main__':main()
