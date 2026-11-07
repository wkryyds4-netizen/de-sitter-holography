"""One-shot exact frame and moving finite-word diagnostics; no old drivers."""
from pathlib import Path
from fractions import Fraction as F
from datetime import datetime,timezone
import hashlib,json,math,platform,sys
import numpy as np
import scipy,sympy as sp,mpmath as mp
from chord_module_frame import raw_basis,raw_ladder,coherent_descendant,act,word,residuals
from thermal_chord_window import gram_factory
ROOT=Path(__file__).resolve().parents[1]
STAGE='susskind-chord-module-frame';OUT=ROOT/f'results/{STAGE}/frame-checks.json'
if OUT.exists():raise FileExistsError(OUT)
OUT.parent.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
counts=dict(gram=0,ladder=0,laguerre=0)
for q in (F(1,4),F(4,9)):
 for d in (1,2,3):
  r=q**d;g=gram_factory(q,r);qi=lambda z:(1-q**z)/(1-q)
  def inner(u,v):return sum(a*b*g(x,y) for x,a in u.items() for y,b in v.items())
  def add(vec,more,c):
   for u,a in more.items():vec[u]=vec.get(u,0)+c*a
   return {u:a for u,a in vec.items() if a}
  for total in range(8):
   states=[raw_basis(q,r,total-k,k) for k in range(total+1)]
   for k,v in enumerate(states):
    n=total-k
    norm=math.prod(qi(i) for i in range(1,n+1))*math.prod(qi(i)*(1-r*r*q**(i-1)) for i in range(1,k+1))
    for z,w in enumerate(states):
     assert inner(v,w)==(norm if z==k else 0);counts['gram']+=1
    for side,creation in [('L',False),('L',True),('R',False),('R',True)]:
     expected={}
     if creation:
      expected=add(expected,raw_basis(q,r,n+1,k),1 if side=='L' else r*q**k)
      if side=='R':expected=add(expected,raw_basis(q,r,n,k+1),1)
     else:
      if n:expected=add(expected,raw_basis(q,r,n-1,k),qi(n)*(1 if side=='L' else r*q**k))
      if k and side=='R':expected=add(expected,raw_basis(q,r,n,k-1),qi(k)*(1-r*r*q**(k-1)))
     assert raw_ladder(q,r,v,side,creation)==expected;counts['ladder']+=1
x=sp.symbols('x')
for d in (1,2,3):
 for k in range(7):
  pol=lambda j:sp.assoc_laguerre(j,2*d-1,2*x) if j>=0 else sp.Integer(0)
  assert sp.expand(2*x*pol(k)-(2*k+2*d)*pol(k)+(k+1)*pol(k+1)+(k+2*d-1)*pol(k-1))==0
  counts['laguerre']+=1

# Bound all j<=2 input tails by (n+m+1)^4/m^2 times a Poisson law.
stirling=[[1],[0,1],[0,1,1],[0,1,3,1],[0,1,7,6,1]]
def chernoff(mu,k):
 return math.exp(-mu+k*(1+math.log(mu/k))) if k>mu else 1.
def input_tail_bound(m,M):
 return sum(math.comb(4,a)*(m+1)**(4-a)*sum(stirling[a][b]*m**b*chernoff(m,M+1-b)
                   for b in range(a+1)) for a in range(5))/m**2

rows=[]
for m in (8.,32.,128.,512.):
 eps=m**-3;cut=math.ceil(m+14*math.sqrt(m+1)+80);fine=cut+48;d=2.
 for j in (0,1,2):
  cv=coherent_descendant(m,j,cut);fv=coherent_descendant(m,j,fine)
  assert abs(float(np.dot(fv,fv))-1)<3e-12
  h=np.zeros_like(fv)
  qn=-np.expm1(-eps*np.arange(fine+1))/(-math.expm1(-eps))
  h[1:]+=np.sqrt(qn[1:])*fv[:-1];h[:-1]+=np.sqrt(qn[1:])*fv[1:]
  eta_h=float(np.linalg.norm(h+2*math.sqrt(m)*fv)/(2*math.sqrt(m)))
  s=math.exp(-eps/2);a=-math.expm1(-eps/2)-s*eps/4
  wbound=math.sqrt(eps/(-math.expm1(-eps)))
  eta=eta_h+wbound/(2*math.sqrt(m))+math.sqrt(a)/(2*math.sqrt(s*m))
  for k in (0,1,2,3):
   v=np.zeros((fine+9,9),complex);v[:fine+1,k]=fv
   u=np.zeros_like(v);u[:cut+1,k]=cv
   res=residuals(eps,d,v)
   ref=max(float(np.linalg.norm(word(eps,d,w,v)-word(eps,d,w,u)))
           for w in ('B','P','E','T','Q','BE','PT','EBP','EEEE'))
   assert ref<1e-9
   rows.append(dict(m=m,epsilon=eps,d=d,j=j,k=k,cutoff=cut,fine_cutoff=fine,
     input_norm=float(np.linalg.norm(v)),input_tail_squared_bound=input_tail_bound(m,cut),
     word_refinement_absolute=ref,residuals=res,
     selected_vector_bound=dict(eta_h=eta_h,W_operator_bound=wbound,
          eta=eta,normalized_vector_distance_upper=2*eta)))
 print(json.dumps(dict(m=m,max_comm=max(r['residuals']['comm_BP']['absolute'] for r in rows if r['m']==m))),flush=True)

# Independent 50-digit sparse direct-coefficient implementation at the largest
# displacement, for a nontrivial common and relative descendant j=2,k=3.
mp.mp.dps=50
m=mp.mpf(512);eps=m**-3;d=mp.mpf(2);q=mp.exp(-eps)
cut=next(r['fine_cutoff'] for r in rows if r['m']==512)
def qmp(z):return -mp.expm1(-eps*z)/(-mp.expm1(-eps))
c=[mp.exp(-m/2)]
for n in range(1,cut+1):c.append(-mp.sqrt(m/n)*c[-1])
vec={(n,3):v*((n-m)**2-n)/(mp.sqrt(2)*m) for n,v in enumerate(c)}
def opmp(name,v):
 if name=='Q':
  out=opmp('E',v)
  for key,z in opmp('P',v).items():out[key]=out.get(key,0)-z
  return out
 out={}
 def put(n,k,z):
  if n>=0 and k>=0:out[n,k]=out.get((n,k),0)+z
 for (n,k),z in v.items():
  R=q**(d+k);cn=q**((n+k+d)/2);root=mp.sqrt(1-cn*cn)
  up=mp.sqrt(qmp(k+1)*(1-q**(2*d+k)))
  down=mp.sqrt(qmp(k)*(1-q**(2*d+k-1)))
  if name=='T':
   put(n+1,k,z*(1-R)*mp.sqrt(qmp(n+1))/(2*mp.sqrt(eps)))
   if n:put(n-1,k,z*(1-R)*mp.sqrt(qmp(n))/(2*mp.sqrt(eps)))
   put(n,k+1,-z*up/(2*mp.sqrt(eps)))
   if k:put(n,k-1,-z*down/(2*mp.sqrt(eps)))
  elif name in ('B','E'):
   if name=='B':diag=((1-R*R)*qmp(n)-down*down)/(2*cn*root)
   else:diag=(qmp(n+k+d)-(1+R*R)*qmp(n)/2-down*down/2)/cn
   put(n,k,z*diag);den=2*cn*(root if name=='B' else 1)
   if n:put(n-1,k+1,-z*R*up*mp.sqrt(qmp(n))/den)
   if k:put(n+1,k-1,-z*q**(d+k-1)*down*mp.sqrt(qmp(n+1))/den)
  elif name=='P':
   if n:put(n-1,k+1,-mp.j*z*up*mp.sqrt(qmp(n))/(2*root))
   if k:put(n+1,k-1,mp.j*z*down*mp.sqrt(qmp(n+1))/(2*root))
  else:raise ValueError(name)
 return out

precision=[]
base=np.zeros((cut+9,9),complex);base[:cut+1,3]=coherent_descendant(float(m),2,cut)
for letters in ('B','P','E','T','Q','BE','EBP','EEEE'):
 v=vec
 for name in reversed(letters):v=opmp(name,v)
 reference=word(float(eps),2,letters,base)
 err=mp.sqrt(mp.fsum(abs(z-mp.mpc(reference[n,k]))**2 for (n,k),z in v.items()))
 norm=mp.sqrt(mp.fsum(abs(z)**2 for z in v.values()))
 assert err/norm<mp.mpf('2e-11')
 precision.append(dict(word=letters,absolute=str(err),relative=str(err/norm),norm=str(norm)))

result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),
 conventions=dict(frame='aL=a_q tensor I; aR=a_q tensor r q^K + I tensor b',
  b='b|k>=sqrt([k](1-r^2 q^(k-1)))|k-1>',moving='F_q,d (D(-sqrt(m)) tensor I)',
  limit='E->L0, Bhat->L1, Phat->i(Lplus-Lminus)/2, Btr->-L1, Q=E-Phat>=0 in limit'),
 exact_counts=counts,rows=rows,precision_50_decimal=precision,
 input_sha256={p:sha(p) for p in ['code/chord_module_frame.py','code/check_chord_module_frame.py',
  'code/thermal_chord_window.py','notes/susskind-thermal-chord-window.md','notes/susskind-slater-selected-state.md']},
 environment=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__,
  mpmath=mp.__version__,precision_decimal=50,platform=platform.platform()),
 scope='Exact one-matter frame, moving finite-word and commutator diagnostics for actual d=2. Full rho vector insertion controlled analytically using uniform W norm, not numerically replaced by energy-only density.',
 limitations='Input Poisson tail bound is exact-arithmetic. Word tail convergence uses manuscript envelopes; refinement and 50-decimal comparison are diagnostics, not interval certificates. Finite Q matrices need not be positive; positivity is proved for the limit. No mass or multiparticle translation law claimed.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',exact_counts=counts,rows=len(rows),sha256=sha(OUT.relative_to(ROOT)),
 max_precision_relative=max(float(r['relative']) for r in precision))))
