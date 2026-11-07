"""Exact charge trace combinatorics and symbolic quotient-space assembly."""
from pathlib import Path
from fractions import Fraction as F
from math import comb,prod
import hashlib,json,platform
import sympy as s
from memory_two_charge_filter import filters
from cubic_linear_dynamics import weights

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-descendants'


def main():
    target=OUT/'assembly-certificate.json'
    if target.exists():raise FileExistsError(target)
    q,lam,L,t=s.symbols('q lam L t')
    leading=s.expand(-(t*(t-1)+(L-t)*(L-t-1))+L*(L-1)/2)
    assert s.expand(leading-(L-(2*t-L)**2)/2)==0
    checks=0
    for N in [10,14,20]:
        m=N//2
        for ll in range(m+1):
            for tt in range(ll+1):
                direct=F(comb(N-ll,m-tt),comb(N,m))*2**ll
                falling=F(2**ll*prod(m-j for j in range(tt))*prod(m-j for j in range(ll-tt)),
                          prod(N-j for j in range(ll)))
                assert direct==falling
                checks+=1
    # Independent mean from diagonal rank filters vs full-CAR second contraction.
    finite=json.loads((OUT/'finite-descendants.json').read_text());meanchecks=0;skipped=[]
    for row in finite['rows']:
        N,p,m=row['N'],row['p'],row['m'];C=comb(m,p)*comb(N-m+p,p)
        if p>min(m-1,N-m):
            skipped.append(dict(N=N,p=p,m=m,reason='Lower-only ranks outside the common-rank filter domain.'))
            continue
        om=weights(N,m,p);omm=weights(N,m-1,p)
        Cm=comb(m-1,p)*comb(N-m+1+p,p)
        means=[F(0)]*5
        for a in range(min(p,m,N-m)+1):
            w=filters(N,p,m,a,a,0)['words']
            coeff=[tuple(w['D22'][i]/2+(2*p-m)*w['D40'][i] for i in range(2)),
                   tuple(w['D121'][i]+(2*p-m)*w['D40'][i] for i in range(2)),
                   tuple(w['D13'][i]/2-w['D22'][i]/2-w['D121'][i]+(m-p)*w['D40'][i] for i in range(2)),
                   tuple(w['D4'][i]-w['D13'][i]+w['D121'][i] for i in range(2)),w['D40']]
            for j,(gp,gm) in enumerate(coeff):
                means[j]+=(gp*C*om[a]+gm*m*Cm*omm[a])/2
        assert means==list(map(F,row['means'])),(row,means)
        meanchecks+=5
    desc=json.loads((OUT/'symbolic-filters.json').read_text())
    assert all(s.sympify(row['S_coefficient'])==0 for row in desc['rows'] if row['name']!='C')
    old=json.loads((ROOT/'results/susskind-memory-quartic/quartic.json').read_text())
    parse=lambda x:s.sympify(x.replace('lambda','lam'),locals={'q':q,'lam':lam})
    G=s.Matrix([[parse(x) for x in row] for row in old['gram']])
    ids=[0,1,2,4];nuisance=G.extract(ids,ids)
    det=s.factor(nuisance.det());assert det!=0
    vec=s.Matrix(list(map(parse,old['residual_vector'])))
    assert all(s.factor(x)==0 for x in nuisance*vec.extract(ids,[0])+G.extract(ids,[3]))
    V4=s.factor(parse(old['residual_variance'])/lam**4)
    Vnew=(1+q)*(1-q*q+8*lam*q)/(lam*(1+q+4*lam*q))
    V2=64*q**4*Vnew
    small2=s.limit(V2.subs(q,s.exp(-4*lam)),lam,0)
    small4=s.limit((lam*V4).subs(q,s.exp(-4*lam)),lam,0)
    assert small2==1024 and small4==s.Rational(768,5)
    paths=['code/check_memory_descendant_assembly.py','code/memory_two_charge_filter.py',
           'code/cubic_linear_dynamics.py','results/susskind-memory-descendants/finite-descendants.json',
           'results/susskind-memory-descendants/symbolic-filters.json',
           'results/susskind-memory-quartic/quartic.json']
    data=dict(status='passed',fixed_charge_trace_identities=checks,rank_mean_identities=meanchecks,rank_mean_skipped_systems=skipped,
              half_filled_log_trace_ratio_leading=str(s.factor(leading)),
              fourth_chaos_nuisance_determinant=str(det),V4_over_Js8=str(V4),
              V2_over_Js8=str(s.factor(V2)),small_lambda_V2=str(small2),small_lambda_lambda_V4=str(small4),
              trace_bridge_scope='Exact occupancy-count ratio checked at finite N. Analytic bridge in manuscript uses balanced singleton modes for Wick blocks with no internal Gaussian chords, bounded intersection count, and the established absolute exponential overlap tail; numeric identities alone do not prove the limit.',
              environment=dict(python=platform.python_version(),sympy=s.__version__),
              input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in paths})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print('passed',checks,'trace identities,',meanchecks,'rank means; determinant',det,flush=True)


if __name__=='__main__':main()
