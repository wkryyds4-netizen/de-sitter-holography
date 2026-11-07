"""Exact two-insertion identity checked against full selected and h/W boundaries."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product,combinations
import hashlib,json,platform,time
from selected_bilocal_covariance import *
from neutral_tower_projection import pair_integral
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/susskind-selected-bilocal-covariance/exact-checks.json'
if OUT.exists():raise FileExistsError(OUT)
OUT.parent.mkdir(parents=True,exist_ok=True)
start=time.monotonic();rows=[];counts={}
def check(kind,lhs,rhs,**meta):
    row=dict(kind=kind,actual=str(lhs),expected=str(rhs),passed=lhs==rhs,**meta)
    rows.append(row);counts[kind]=counts.get(kind,0)+1
    if lhs!=rhs:
        OUT.with_name('failed-exact-checks.json').write_text(json.dumps(rows,indent=2)+'\n')
        raise AssertionError(row)

for gamma,lam in [(F(1,2),F(1,3)),(F(2,3),F(1,5))]:
 q=gamma**4;m=JointModel(q,lam)
 fields={f'L{n}_{v}':(*centered_prefix(n,v,m),n) for n in range(1,4) for v in range(1,n+1)}
 for d in (2,3):fields[f'J{d}']=(*centered_prefix(d,d,m,'D'),d)
 for j,(field,phi) in enumerate(neutral_fields(4,m)):
  if j:fields[f'N{j}']=(field,phi,j)
 selected_count=0
 for nx,ny in product(fields,repeat=2):
  x,px,dx=fields[nx];y,py,dy=fields[ny]
  if dx+dy>6:continue
  gaps=[(0,0)] if (dx+dy)%2==0 else [(1,0),(0,1)]
  for a,b in gaps:
   if dx+dy+a+b>6:continue
   direct=field_pair(x,y,a,b,gamma,lam,'selected')
   tilted=tilted_pair(x,y,a,b,gamma,lam)
   correction=chain_correction(px,py,a,b,q,lam,gamma**2)
   check('full_selected',direct-tilted,correction,x=nx,y=ny,a=a,b=b,q=str(q))
   selected_count+=1
 # Separate boundary mechanisms on mixed and lower-chaos cases.
 for nx,ny,a,b in [('L2_1','L2_2',0,0),('L3_2','L1_1',0,0),
                  ('J2','L2_2',0,0),('J3','L1_1',0,0),
                  ('N3','N1',0,0),('N2','N2',0,0),
                  ('N1','L2_1',1,0),('N1','L2_1',0,1)]:
  x,px,_=fields[nx];y,py,_=fields[ny]
  odd=boundary_pair(x,y,a,b,'hW',gamma,lam)+boundary_pair(x,y,a,b,'Wh',gamma,lam)
  even=boundary_pair(x,y,a,b,'WW',gamma,lam)-lam*field_pair(x,y,a,b,gamma,lam)
  check('odd_boundary',odd,F(0),x=nx,y=ny,a=a,b=b,q=str(q))
  check('even_boundary',even,chain_correction(px,py,a,b,q,lam,F(1)),x=nx,y=ny,a=a,b=b,q=str(q))
 # Nontrivial missed residual: charged curvature cannot be replaced by its zero projection.
 j2,_,_=fields['J2']
 norm=field_pair(j2,j2,0,0,gamma,lam,'selected')
 check('curvature_marginal',norm,4*lam+16*lam*lam*q/(1+q),q=str(q))
 assert norm>0
 # A residual Gram matrix in the energy-tilted positive state.
 names=['J2','L2_1','L2_2'];gram=[]
 for nx in names:
  rr=[]
  for ny in names:
   x,px,_=fields[nx];y,py,_=fields[ny]
   projected=F(0)
   for (a,b),v in px.items():
    for (c,d),w in py.items():
     from neutral_selected_response import pair_moment
     projected+=lam*v*w*((1-gamma**2)*pair_moment(a+c,b+d,q)+gamma**2*pair_moment(a+c,b+d+2,q))
   rr.append(tilted_pair(x,y,0,0,gamma,lam)-projected)
  gram.append(rr)
 import sympy as sp
 mat=sp.Matrix([[sp.Rational(v.numerator,v.denominator) for v in rr] for rr in gram])
 assert mat==mat.T
 minors=[]
 for k in (1,2,3):
  for ids in combinations(range(3),k):
   det=mat.extract(ids,ids).det();assert det>=0,(q,ids,det)
   minors.append(dict(indices=ids,value=str(det)))
 rows.append(dict(kind='positive_residual_Gram',q=str(q),names=names,matrix=[[str(v) for v in rr] for rr in gram],principal_minors=minors,passed=True))
 counts['positive_residual_Gram']=counts.get('positive_residual_Gram',0)+1
 print(json.dumps({'q':str(q),'counts':counts,'selected_cases_this_parameter':selected_count}),flush=True)

inputs=['code/selected_bilocal_covariance.py','code/check_selected_bilocal_covariance.py','code/neutral_tower_projection.py',
        'code/slater_higher_wick_sources.py','code/slater_composite_sources.py',
        'code/memory_joint_regression.py','code/memory_wick_descendant_chords.py','code/neutral_selected_response.py']
result=dict(status='passed',counts=counts,rows=rows,seconds=time.monotonic()-start,
 environment=dict(python=platform.python_version(),sympy=sp.__version__,arithmetic='Fraction/SymPy Rational'),
 scope='Finite two-field identities with true selected boundary, independent even/odd boundary expansion, actual lower-chaos fields and nonzero residual Gram. General source proof is in the manuscript.',
 parameters='Independent rational q/lambda identities, not the physical exponential curve.',
 input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({'status':'passed','counts':counts,'seconds':result['seconds'],'sha256':hashlib.sha256(OUT.read_bytes()).hexdigest()}),flush=True)
