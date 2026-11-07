"""Exact spectral polynomial/Taylor checks and fixed-lambda time quadrature."""
from pathlib import Path
from fractions import Fraction as F
from math import factorial,sqrt,log,exp
from functools import lru_cache
import hashlib,json,platform,time
import numpy as np
import scipy
from scipy.linalg import eigh_tridiagonal
from numpy.polynomial.legendre import leggauss
from neutral_tower_projection import *
from slater_higher_wick_sources import moment
from memory_wick_descendant_chords import descendants
from memory_joint_regression import JointModel
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/susskind-selected-bilocal-mean/time-checks.json'
if OUT.exists():raise FileExistsError(OUT)

def mulpoly(a,b):
 out={}
 for (i,j),v in a.items():
  for (k,l),w in b.items():out[(i+k,j+l)]=out.get((i+k,j+l),F(0))+v*w
 return {k:v for k,v in out.items() if v}

def qfac(n,q):
 out=F(1)
 for k in range(1,n+1):out*=qi(k,q)
 return out

def cky(a,b,q):
 """Independently normalized explicit CKY sum, Eq. (1.1), r=q."""
 out={}
 for k in range(min(a,b)+1):
  coef=(-1)**k*q**(k*(k-1)//2+k)*qfac(a,q)*qfac(b,q)/(qfac(a-k,q)*qfac(b-k,q)*qfac(k,q))
  for i,v in hermite(a-k,q).items():
   for j,w in hermite(b-k,q).items():out=add(out,{(i,j):coef*v*w})
 return out

def exact_time(gamma,lam):
 q=gamma**4;model=JointModel(q,lam);_,V,_=descendants(4)
 _,Xi,_,reg=tower(4,q,lam)
 @lru_cache(None)
 def actual(j,a,b):
  total=F(0)
  for term,c in model.neutral(V[j]).items():
   marks=tuple((x+a,y+a,k) for x,y,k in term.marks)
   wr=tuple((x+a,y+a) for x,y in term.wicks)
   ds=[None] if not term.filters else range(term.length%2,term.length+1,2)
   for d in ds:
    weight=F(1);cr=[]
    if d is not None:
     for x,y,p,sel in term.filters:
      if sel>=0 and d!=sel:weight=0
      weight*=q**(p*d);cr.append((x+a,y+a,d))
    total+=c*weight*moment(a+term.length+b,marks,gamma,lam,wick_regions=wr,chaos_regions=tuple(cr),mode='selected')
  return total-sum((c*moment(a+k+b,(),gamma,lam,mode='selected') for k,c in energy_polynomial(reg[j],q).items()),F(0))
 rows=[]
 for m in (0,2,4):
  for ad in (0,2) if m<=2 else (0,):
   lhs=rhs=F(0)
   for j in range(1,min(4,m)+1):
    for k in range(j,m+1,2):
     c=power_coeff(k,q).get(j,F(0))/factorial(k)
     for a in range(m-k+1):
      b=m-k-a;coef=2*c*F((-1)**(m//2+a+b),2**(a+b)*factorial(a)*factorial(b))
      for z in range(ad+1):
       w=F((-1)**(ad-z)*factorial(ad),factorial(z)*factorial(ad-z))
       lhs+=coef*w*actual(j,a+z,b+ad-z)
       poly={(x+a+z,y+b+ad-z):v for (x,y),v in to_symbol(Xi[j],q).items()}
       rhs+=coef*w*selected_symbol(poly,q,lam,gamma)
   assert lhs==rhs,(m,ad,lhs,rhs)
   rows.append(dict(m=m,ad=ad,value=str(lhs),independent=str(rhs)))
 f2=next(F(r['value']) for r in rows if r['m']==2 and r['ad']==0)
 b4=-f2/12
 expected=F(2,3)*gamma**2*lam**2*q*(1-q*q)
 assert b4==expected
 f4=next(F(r['value']) for r in rows if r['m']==4)
 ad2=next(F(r['value']) for r in rows if r['m']==2 and r['ad']==2)
 return dict(q=str(q),lambda_=str(lam),rows=rows,bilocal_t4=str(b4),bilocal_t6=str(-f4/30+ad2/1440))

def degree_bound(q,lam,t,J):
 r=1/sqrt(1-q);a=q*r*abs(t);b=3*r*abs(t);force=4*lam*q/(1-q)
 assert a<J+1 and b<J+2
 neutral=sqrt(lam)*(4*lam*a**(J+1)/factorial(J)/(1-a/(J+1))+
                    force*b**(J+1)/factorial(J+1)/(1-b/(J+2)))
 return t*t*sqrt(1+q*(1+q)*(1+lam)**2)*neutral

def numerical(q,t,T,J,L,ns):
 lam=-log(q)/4;R=2/sqrt(1-q)
 e,v=eigh_tridiagonal(np.zeros(L),np.sqrt((1-q**np.arange(1,L))/(1-q)))
 v*=np.where(v[0]<0,-1,1)[None,:];root=v[0]
 mass=root[:,None]*root[None,:]*(v.T@(q**np.arange(L)[:,None]*v))
 x=e[:,None];y=e[None,:];delta=x-y;plus=x+y
 P={(0,0):np.ones((L,L))}
 for b in range(1,J):P[0,b]=y*P[0,b-1]-(float(qi(b-1,q))*P[0,b-2] if b>1 else 0)
 for a in range(1,J):
  for b in range(J-a):P[a,b]=x*P[a-1,b]-(float(qi(a-1,q))*P[a-2,b] if a>1 else 0)-(q**a*float(qi(b,q))*P[a-1,b-1] if b else 0)
 _,Xi,_,_=tower(J,q,lam)
 X=[np.zeros((L,L))]+[sum(float(c)*P[ab] for ab,c in Xi[j].items()) for j in range(1,J+1)]
 # Stable c_j: start its scalar Taylor sum at order j, not an eigensolver cancellation.
 M=80;ac=np.zeros((M+1,M+1));ac[0,0]=1
 for k in range(M):
  for j in range(k+1):
   ac[k+1,j+1]+=ac[k,j]
   if j:ac[k+1,j-1]+=ac[k,j]*(1-q**j)/(1-q)
 z,w=leggauss(ns);nodes=(z+1)*t/2;weights=w*t/2;answer=0j;absolute=0.
 for s,weight in zip(nodes,weights):
  c=[sum((1j*s)**k*ac[k,j]/factorial(k) for k in range(j,M+1,2)) for j in range(J+1)]
  psi=sum(c[j]*X[j] for j in range(1,J+1))
  integrand=mass*plus*np.exp(1j*T*delta-.5j*s*plus)*psi*(t-s)*np.sinc((t-s)*delta/(2*np.pi))
  answer+=weight*integrand.sum();absolute+=abs(weight)*np.abs(integrand).sum()
 answer*=2*sqrt(q)*lam;absolute*=2*sqrt(q)*lam
 return dict(q=q,lambda_=lam,t=t,T=T,J=J,L=L,time_nodes=ns,value=[float(answer.real),float(answer.imag)],
             cancellation=float(absolute/max(abs(answer),1e-300)),degree_tail=degree_bound(q,lam,t,J),
             c_scalar_Taylor_remainder=(R*t)**(M+1)/factorial(M+1),mass_error=float(abs(mass.sum()-1)))

start=time.monotonic();symbol_count=0
for q in (F(1,4),F(2,3)):
 for a in range(8):
  for b in range(8-a):
   assert symbol(a,b,q)==cky(a,b,q),(a,b,q)
   assert symbol(a,b,q)=={(y,x):v for (x,y),v in symbol(b,a,q).items()}
   symbol_count+=1
exact=[exact_time(F(1,2),F(1,3)),exact_time(F(2,3),F(1,5))]
print(json.dumps({'phase':'exact','symbols':symbol_count,'time_cases':len(exact)}),flush=True)
rows=[]
for q in (.25,.5):
 for t,T in [(.1,0),(.25,0),(.5,0),(1.,0),(.5,.5),(1.,.5)]:
  coarse=numerical(q,t,T,24,48,24);fine=numerical(q,t,T,32,64,32)
  difference=abs(complex(*coarse['value'])-complex(*fine['value']))
  assert difference<2e-11,(coarse,fine,difference)
  fine['coarse']=coarse;fine['refinement_difference']=difference
  lam=fine['lambda_'];b4=2/3*sqrt(q)*lam**2*q*(1-q*q)
  if T==0:fine['relative_to_t4']=fine['value'][0]/(b4*t**4)
  rows.append(fine)
 print(json.dumps({'phase':'time','q':q,'points':len(rows)}),flush=True)
inputs=['code/neutral_tower_projection.py','code/check_selected_bilocal_mean.py','code/slater_higher_wick_sources.py',
        'code/memory_wick_descendant_chords.py','code/memory_joint_regression.py','code/neutral_selected_response.py']
result=dict(status='passed',symbol_checks=symbol_count,exact_time=exact,numerical=rows,
 scope='Exact polynomial and Taylor identities; quadrature refinements are diagnostics, not interval-certified total errors. Analytic degree tails assume exact arithmetic.',
 environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,dtype='float64/complex128'),
 seconds=time.monotonic()-start,input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({'status':'passed','seconds':result['seconds'],'sha256':hashlib.sha256(OUT.read_bytes()).hexdigest()}),flush=True)
