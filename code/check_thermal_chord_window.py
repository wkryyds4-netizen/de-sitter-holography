"""One-shot verification of the thermal endpoint window; refuses overwrite."""
from pathlib import Path
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib, json, math, platform, sys
import numpy as np
import scipy, sympy as sp
from scipy.stats import poisson
from thermal_chord_window import (gram_factory, heat_vector, coherent,
                                  endpoint_actions, observables)

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/susskind-thermal-chord-window/window-checks.json'
if OUT.exists(): raise FileExistsError(OUT)
OUT.parent.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()

counts=dict(gram_symmetry=0, adjoints=0, endpoint_gram=0, generator_actions=0,
            endpoint_norms=0, boost_actions=0, moment_comparisons=0)
for q in (F(1,4),F(4,9)):
  for d in (1,2,3):
    r=q**d; gram=gram_factory(q,r)
    qi=lambda n:(1-q**n)/(1-q)
    fac=lambda n:math.prod(qi(j) for j in range(1,n+1))
    def lower(v,side):
      l,k=v;out={}
      if l:out[l-1,k]=qi(l)*(1 if side=='L' else r*q**k)
      if k:out[l,k-1]=qi(k)*(r*q**l if side=='L' else 1)
      return out
    def raisev(v,side):return (v[0]+(side=='L'),v[1]+(side=='R'))
    def apply(v,i,j):
      return {raisev(u,i):c for u,c in lower(v,j).items()}
    def inner(x,y):
      return sum(sp.conjugate(a)*b*sp.Rational(gram(u,v))
                 for u,a in x.items() for v,b in y.items())
    for n in range(8):
      states=[(n-j,j) for j in range(n+1)]
      for u in states:
        for v in states:
          assert gram(u,v)==gram(v,u);counts['gram_symmetry']+=1
          for side in ('L','R'):
            lhs=gram(raisev(u,side),raisev(v,side))
            rhs=sum(c*gram(u,w) for w,c in lower(raisev(v,side),side).items())
            assert lhs==rhs;counts['adjoints']+=1
      assert gram((n,0),(n,0))==fac(n);counts['endpoint_gram']+=1
      if n:
        assert gram((n,0),(n-1,1))==r*fac(n)
        assert gram((n-1,1),(n-1,1))/fac(n)==r*r+(1-r*r)/qi(n)
        counts['endpoint_gram']+=2
      e={(n,0):sp.Integer(1)}
      c=sp.Rational(q)**sp.Rational(n+d,2)
      jj={}
      for i,j in [('L','L'),('L','R'),('R','L'),('R','R')]:
        z={u:sp.Rational(v) for u,v in apply((n,0),i,j).items()}
        z[n,0]=z.get((n,0),0)-sp.Rational(qi(n+d))
        jj[i+j]=z
      keys=set().union(*[set(v) for v in jj.values()])
      actions={
        'B':{u:(jj['LL'].get(u,0)-jj['RR'].get(u,0))/(2*c) for u in keys},
        'E':{u:-(jj['LL'].get(u,0)+jj['RR'].get(u,0))/(2*c) for u in keys},
        'P':{u:sp.I*(jj['LR'].get(u,0)-jj['RL'].get(u,0))/2 for u in keys}}
      nn,rr=map(sp.Rational,(qi(n),r))
      predicted={'B':{(n,0):nn/(2*c)},'P':{(n,0):sp.I*rr*nn/2},
                 'E':{(n,0):(sp.Rational(qi(n+d))-nn/2)/c}}
      if n:
        predicted['B'][n-1,1]=-rr*nn/(2*c)
        predicted['P'][n-1,1]=-sp.I*nn/2
        predicted['E'][n-1,1]=-rr*nn/(2*c)
      for name in actions:
        assert all(sp.simplify(actions[name].get(u,0)-predicted[name].get(u,0))==0
                   for u in keys|set(predicted[name]));counts['generator_actions']+=1
      # Independent nonorthogonal Gram versus orthogonal e/f formulas.
      bg=nn*(1-rr**2)/(2*c)
      bf2=rr**2*nn*(1-rr**2)/(4*c**2)
      eg=(sp.Rational(qi(d))-(1-rr)**2*nn/2)/c
      norms={'B':bg**2+bf2,'E':eg**2+bf2,'P':nn*(1-rr**2)/4}
      for name,expected in norms.items():
        assert sp.simplify(inner(actions[name],actions[name])/fac(n)-expected)==0
        counts['endpoint_norms']+=1
      # sqrt(epsilon) Btr: both degree-changing components, in raw basis.
      raw={(n+1,0):sp.Rational(1,2),(n,1):sp.Rational(-1,2)}
      if n:raw[n-1,0]=(1-rr)*nn/2
      expected=((1-rr)**2*(sp.Rational(qi(n+1))+nn)+(1-rr**2))/4
      assert sp.simplify(inner(raw,raw)/fac(n)-expected)==0
      counts['boost_actions']+=1

# Independent Wick pairings versus q-Jacobi moments; no inherited driver.
def pairs(items):
  if not items:yield [];return
  a=items[0]
  for j in range(1,len(items)):
    b=items[j]
    for rest in pairs(items[1:j]+items[j+1:]):yield [(a,b)]+rest
moment_rows=[]
for q in (F(1,4),F(4,9)):
  for k in range(1,7):
    terms=[]
    for pp in pairs(list(range(2*k))):
      cr=sum(a<c<b<e or c<a<e<b for j,(a,b) in enumerate(pp) for c,e in pp[j+1:])
      terms.append(q**cr)
    vec={0:F(1)}
    for step in range(2*k):
      nxt={}
      for n,v in vec.items():
        nxt[n+1]=nxt.get(n+1,0)+v
        if n:nxt[n-1]=nxt.get(n-1,0)+v*(1-q**n)/(1-q)
      vec=nxt
    assert vec[0]==sum(terms);counts['moment_comparisons']+=1
    moment_rows.append(dict(q=str(q),k=k,pairings=len(terms),moment=str(vec[0])))

rows=[]
for m in (8.,32.,128.,512.):
  epsilon=m**-3;delta=2*epsilon*m*m
  cutoff=math.ceil(2*m+14*math.sqrt(2*m)+50)
  fine=cutoff+64
  pa,ma=heat_vector(epsilon,m,cutoff,'taylor')
  pb,mb=heat_vector(epsilon,m,fine,'taylor')
  pc,mc=heat_vector(epsilon,m,fine,'taylor-longdouble')
  pd,md=heat_vector(epsilon,m,fine,'expm')
  pad=np.pad(pa,(0,64))
  cross=float(np.linalg.norm(pb-pc));refine=float(np.linalg.norm(pb-pad))
  print(json.dumps(dict(diagnostic_m=m,precision=cross,refine=refine,
                       split_expm_difference=float(np.linalg.norm(pb-pd)))),flush=True)
  assert cross<2e-11 and refine<2e-11
  chi=coherent(m,fine)
  distance=float(np.linalg.norm(pb-chi))
  bound=math.sqrt(2*(1-math.sqrt(1-delta)))
  assert distance<=bound+1e-12
  n=np.arange(fine+1);p=pb**2
  tail=poisson.sf(cutoff,2*m)
  cutoff_norm_bound=math.sqrt(2*tail/(1-delta))
  assert mb['scaled_norm_squared']>=1-delta-1e-11
  assert mb['scaled_norm_squared']<=1+1e-11
  assert np.all(p <=chi**2/(1-delta)+2e-14)
  measures={}
  for d in (1.,2.,3.):
    ob=observables(epsilon,d,pb)
    oa=observables(epsilon,d,pad)
    oc=observables(epsilon,d,pc)
    acts_b=endpoint_actions(epsilon,d,pb)
    acts_c=endpoint_actions(epsilon,d,pc)
    ob['representation_action_difference_max']=max(float(np.linalg.norm(acts_b[k]-acts_c[k]))
                                                for k in ('B','P','E','Btr'))
    ob['cutoff_boost_norm_difference']=abs(ob['boost_sum_norm']-oa['boost_sum_norm'])
    ob['relative_boost_representation_difference']=abs(ob['boost_sum_norm']-oc['boost_sum_norm'])/ob['boost_sum_norm']
    assert ob['representation_action_difference_max']<2e-10
    measures[str(int(d))]=ob
  rows.append(dict(m=m,epsilon=epsilon,beta=2*math.sqrt(epsilon*m),delta=delta,
     cutoff=cutoff,fine_cutoff=fine,taylor=mb,taylor_longdouble=mc,
     split_expm_diagnostic=md,split_expm_norm_difference=float(np.linalg.norm(pb-pd)),
     longdouble_mantissa_bits=int(np.finfo(np.longdouble).nmant),
     representation_norm_difference=cross,cutoff_norm_difference=refine,
     exact_arithmetic_cutoff_norm_bound=cutoff_norm_bound,
     poisson_2m_tail=float(tail),coherent_distance=distance,coherent_distance_bound=bound,
     number_mean=float(n@p),number_variance=float(((n-n@p)**2)@p),
     concentration_mse=float(((n/m-1)**2)@p),seed_overlap=float(pb[0]),
     observables=measures))

result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),
 scope='Endpoint one-insertion tracial module; d=2 is the actual charged curvature channel. No selected rho insertion, operator convergence, full sl2 commutator-norm or mass/finite-G claim.',
 conventions=dict(epsilon='4 lambda = -log(q)',K='H/(2 Js)',Btr='(K_L-K_R)/2',
    state='normalized exp(-beta K_L/2)|d;0,0>',Bhat='B/sqrt(1-q^(n+d))',
    window='m=beta^2/(4 epsilon) -> infinity, epsilon m^2 -> 0'),
 exact_counts=counts,moment_rows=moment_rows,rows=rows,
 environment=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,
                  sympy=sp.__version__,platform=platform.platform()),
 input_sha256={p:sha(p) for p in ['code/thermal_chord_window.py','code/check_thermal_chord_window.py',
             'notes/susskind-translations.md','notes/susskind-preparation-modular-clock.md']},
 limitations='Positive Taylor tail and Jacobi cutoff bounds are exact-arithmetic bounds; floating-point/cancellation assessed by 53/64-bit significand comparison and cutoff refinement, not interval certified. Split scipy expm is a diagnostic with known lost-tail amplification, not a high-accuracy reference. One-vector matching is not operator equality.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',exact_counts=counts,rows=len(rows),sha256=sha(OUT.relative_to(ROOT)),
 max_cross=max(z['representation_norm_difference'] for z in rows),
 d2=[dict(m=z['m'],**z['observables']['2']) for z in rows])))
