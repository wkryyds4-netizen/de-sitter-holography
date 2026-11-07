"""One-shot exact ordered Gram/ladders and joint moving-word check."""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from datetime import datetime,timezone
import math,json,hashlib,platform
import numpy as np
import scipy,mpmath as mp
from chord_module_frame import coherent_descendant
from multimatter_moving_frame import raw_image,exact_ladder,exact_inner,action,word,norm
ROOT=Path(__file__).resolve().parents[1]
STAGE='susskind-multimatter-moving-frame'
OUT=ROOT/f'results/{STAGE}/checks.json'
if OUT.exists():raise FileExistsError(OUT)
OUT.parent.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def basis(n):return [(l,m,n-l-m) for l in range(n+1) for m in range(n-l+1)]
counts=dict(gram=0,outer_ladder=0)
for q in (F(1,4),F(4,9)):
 for ds in ((2,2),(1,2)):
  def lower(v,side):
   l,m,r=v
   fac=(1,q**(ds[0]+l),q**(sum(ds)+l+m)) if side==0 else (q**(sum(ds)+m+r),q**(ds[1]+r),1)
   out={}
   for j,z in enumerate(v):
    if z:
     u=list(v);u[j]-=1;out[tuple(u)]=(1-q**z)/(1-q)*fac[j]
   return out
  @lru_cache(None)
  def gram(u,v):
   if sum(u)!=sum(v):return F(0)
   for i in (0,2):
    if u[i]:
     z=list(u);z[i]-=1
     return sum(c*gram(tuple(z),t) for t,c in lower(v,i).items())
   for i in (0,2):
    if v[i]:
     z=list(v);z[i]-=1
     return sum(c*gram(t,tuple(z)) for t,c in lower(u,i).items())
   return math.prod((1-q**i)/(1-q) for i in range(1,u[1]+1))
  for degree in range(5):
   images={z:raw_image(q,ds,z) for z in basis(degree)}
   for u,vu in images.items():
    for v,vv in images.items():
     assert exact_inner(q,ds,vu,vv)==gram(u,v),(q,ds,u,v)
     counts['gram']+=1
    for side in (0,2):
     expected={}
     for z,c in lower(u,side).items():
      for t,a in raw_image(q,ds,z).items():expected[t]=expected.get(t,0)+c*a
     expected={t:a for t,a in expected.items() if a}
     assert exact_ladder(q,ds,vu,side)==expected
     counts['outer_ladder']+=1
print(json.dumps(counts),flush=True)

stirling=[[1],[0,1],[0,1,1],[0,1,3,1],[0,1,7,6,1]]
def tail(m,M):
 def bound(k):return math.exp(-m+k*(1+math.log(m/k))) if k>m else 1
 return sum(math.comb(4,a)*(m+1)**(4-a)*sum(stirling[a][b]*m**b*bound(M+1-b)
            for b in range(a+1)) for a in range(5))/m**2
letters=('E','B','P','T','Q','BE','PT','EBP')
rows=[]
for m in (8.,32.,128.,512.):
 eps=m**-3;cut=math.ceil(m+14*math.sqrt(m+1)+80);fine=cut+48
 for j,k,l in ((0,0,0),(1,1,0),(2,1,2)):
  ds=(2.,2.)
  v=np.zeros((fine+9,8,8),complex);v[:fine+1,k,l]=coherent_descendant(m,j,fine)
  u=np.zeros_like(v);u[:cut+1,k,l]=coherent_descendant(m,j,cut)
  assert abs(norm(v)-1)<3e-12
  errors={};ref=0
  for s in letters:
   got=word(eps,ds,s,v);expected=word(eps,ds,s,v,True)
   errors[s]=dict(absolute=norm(got-expected),relative=norm(got-expected)/norm(expected))
   ref=max(ref,norm(got-word(eps,ds,s,u)))
  comms={}
  for a,b,c in (('B','P','E'),('T','Q','Q')):
   x=word(eps,ds,a+b,v);y=word(eps,ds,b+a,v);z=1j*action(eps,ds,c,v)
   defect=norm(x-y-z)
   comms[a+b]=dict(absolute=defect,relative=defect/norm(z),cancellation=(norm(x)+norm(y)+norm(z))/defect)
  assert ref<1e-9
  rows.append(dict(m=m,epsilon=eps,ds=ds,j=j,k1=k,k2=l,cutoff=cut,fine_cutoff=fine,
    input_norm=norm(v),input_tail_squared_bound=tail(m,cut),word_refinement_absolute=ref,
    word_errors=errors,commutators=comms))
 print(json.dumps(dict(m=m,Q=rows[-1]['word_errors']['Q'],TQ=rows[-1]['commutators']['TQ'])),flush=True)

# Independent 60-decimal direct outer-ladder composition, avoiding the
# cancellation-free rearrangement used by numpy. Includes commutator itself.
mp.mp.dps=60;m=mp.mpf(128);eps=m**-3;q=mp.exp(-eps);ds=(mp.mpf(2),mp.mpf(2))
cut=next(r['fine_cutoff'] for r in rows if r['m']==128)
co=[mp.exp(-m/2)]
for i in range(1,cut+1):co.append(-mp.sqrt(m/i)*co[-1])
vec={(i,1,2):x*((i-m)**2-i)/(mp.sqrt(2)*m) for i,x in enumerate(co)}
def qn(z):return -mp.expm1(-eps*z)/(-mp.expm1(-eps))
def add(*parts):
 out={}
 for factor,v in parts:
  for z,c in v.items():out[z]=out.get(z,mp.mpf(0))+factor*c
 return out
def ladder(v,side,up=False):
 out={}
 for state,c in v.items():
  for i in range(side+1):
   z=state[i]
   if not up and not z:continue
   x=z+1 if up else z
   a=mp.sqrt(qn(x))
   if i:a*=mp.sqrt(1-q**(2*ds[i-1]+x-1))
   a*=q**sum(ds[t-1]+state[t] for t in range(i+1,side+1))
   target=list(state);target[i]+=1 if up else -1;target=tuple(target)
   out[target]=out.get(target,0)+a*c
 return out
def op(name,v):
 if name=='Q':return add((1,op('E',v)),(-1,op('P',v)))
 if name=='T':return {z:c/(2*mp.sqrt(eps)) for z,c in add((1,ladder(v,0)),(1,ladder(v,0,True)),(-1,ladder(v,2)),(-1,ladder(v,2,True))).items()}
 ll=ladder(ladder(v,0),0,True);rr=ladder(ladder(v,2),2,True)
 if name=='E':
  base=add((1,{z:qn(sum(z)+4)*c for z,c in v.items()}),(-mp.mpf('.5'),ll),(-mp.mpf('.5'),rr))
  return {z:c/q**((sum(z)+4)/mp.mpf(2)) for z,c in base.items()}
 if name=='B':base=add((1,ll),(-1,rr))
 elif name=='P':base=add((mp.j,ladder(ladder(v,2),0,True)),(-mp.j,ladder(ladder(v,0),2,True)))
 else:raise ValueError(name)
 return {z:c/(2*mp.sqrt(1-q**(sum(z)+4))*(q**((sum(z)+4)/mp.mpf(2)) if name=='B' else 1)) for z,c in base.items()}
def mw(s,v):
 for a in reversed(s):v=op(a,v)
 return v
base=np.zeros((cut+9,8,8),complex);base[:cut+1,1,2]=coherent_descendant(128,2,cut)
precision=[]
for name in ('E','B','P','T','Q','BE','BP_defect','TQ_defect'):
 if name.endswith('_defect'):
  a,b=name[:2];c='E' if a=='B' else 'Q'
  v=add((1,mw(a+b,vec)),(-1,mw(b+a,vec)),(-mp.j,op(c,vec)))
  nv=word(float(eps),(2,2),a+b,base)-word(float(eps),(2,2),b+a,base)-1j*action(float(eps),(2,2),c,base)
 else:v=mw(name,vec);nv=word(float(eps),(2,2),name,base)
 err=mp.sqrt(mp.fsum(abs(z-mp.mpc(nv[key]))**2 for key,z in v.items()))
 scale=mp.sqrt(mp.fsum(abs(z)**2 for z in v.values()))
 assert err/scale<mp.mpf('1e-10'),(name,err/scale)
 precision.append(dict(name=name,absolute=str(err),relative=str(err/scale),reference_norm=str(scale)))
data=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),exact_counts=counts,
 rows=rows,precision_60_decimal=precision,
 environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__,platform=platform.platform()),
 input_sha256={p:sha(p) for p in ['code/multimatter_moving_frame.py','code/check_multimatter_moving_frame.py',
 'code/chord_module_frame.py','notes/susskind-chord-module-frame.md','notes/susskind-two-matter.md']},
 scope='Actual ordered direct Gram and outer ladders, followed by joint moving-frame additive LS limit. Same-species crossed metric is not replaced by the ordered norm.',
 limitations='Fixed word and fixed internal levels only; numerical refinements are not interval bounds. Input norm tail is not full word tail; analytic envelopes in proof. No mass, opposite momentum or finite-G four-point result.')
with OUT.open('x',encoding='utf-8') as f:json.dump(data,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',sha256=sha(OUT.relative_to(ROOT)),max_precision_relative=max(float(z['relative']) for z in precision))))
