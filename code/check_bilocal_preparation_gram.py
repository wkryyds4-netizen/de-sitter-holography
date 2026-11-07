"""Finite-time preparation contribution, refinement and short-time calibration."""
from pathlib import Path
from math import sqrt
import hashlib,json,platform,time
import numpy as np
import scipy
from bilocal_preparation_gram import PreparationGram
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/susskind-selected-bilocal-covariance/gram-checks.json'
if OUT.exists():raise FileExistsError(OUT)
def cc(z):return [float(z.real),float(z.imag)]
start=time.monotonic();diagonal=[];matrices=[];short=[];old_mean=[]
prior=json.loads((ROOT/'results/susskind-selected-bilocal-mean/time-checks.json').read_text())
for q in (.25,.5):
 k=PreparationGram(q,32,64);co=PreparationGram(q,24,48)
 saved={};co_saved={}
 for t in (.25,.5,1.):
  base=k.base_block(t,32);base0=co.base_block(t,24)
  for T in (0.,.5,1.,2.,4.):
   a=k.state(base,t,T);b=co.state(base0,t,T);saved[t,T]=a;co_saved[t,T]=b
   fine=k.excess(a,a);rough=co.excess(b,b)
   diff=abs(fine['connected']-rough['connected'])
   assert diff<3e-12,(q,t,T,diff)
   reference=k.excess(saved[t,0.],saved[t,0.])
   diagonal.append(dict(q=q,lambda_=k.lam,t=t,T=T,
       uncentered=cc(fine['uncentered']),connected_excess=cc(fine['connected']),
       coarse_connected_excess=cc(rough['connected']),mean=cc(a['mean']),
       total_variance_shift_from_T0=cc(fine['connected']-reference['connected']),
       refinement_absolute=diff,refinement_relative=diff/max(abs(fine['connected']),1e-300),
       propagated_degree_bound=fine['propagated_degree_bound'],
       beta_norm=a['beta_norm'],analytic_beta_bound=a['analytic_beta_bound'],
       subtraction_condition=fine['subtraction_condition'],mean_cancellation=a['mean_cancellation']))
   for old in prior['numerical']:
    if old['q']==q and old['t']==t and old['T']==T:
     error=abs(a['mean']-complex(*old['value']));assert error<3e-14
     old_mean.append(dict(q=q,t=t,T=T,difference=error))
 # General unequal-time Gram, retaining complex phases and centered subtraction.
 labels=[(.5,0.),(.5,.5),(.5,1.),(1.,0.),(1.,.5),(1.,1.)]
 G=np.zeros((6,6),complex);D=G.copy();maxdiff=0.
 for i,lab in enumerate(labels):
  for j,lbc in enumerate(labels):
   row=k.excess(saved[lab],saved[lbc]);other=co.excess(co_saved[lab],co_saved[lbc])
   G[i,j]=row['uncentered'];D[i,j]=row['connected']
   maxdiff=max(maxdiff,abs(row['connected']-other['connected']))
 assert np.max(abs(G-G.conj().T))<1e-14
 eig=np.linalg.eigvalsh(G);assert eig.min()>-1e-13
 matrices.append(dict(q=q,labels=labels,uncentered=[[cc(z) for z in row] for row in G],
                      connected_excess=[[cc(z) for z in row] for row in D],
                      uncentered_eigenvalues=eig.tolist(),
                      connected_excess_eigenvalues=np.linalg.eigvalsh(D).tolist(),
                      max_refinement_absolute=float(maxdiff)))
 # The t^7 term is present away from center zero; do not impose even-t variance.
 T=.5;lead=16/9*sqrt(q)*k.lam**4*q*q*(1+q)
 c=k.C2(T).real;dc=k.C2(T,1).real
 for t in (.1,.05,.025):
  a=k.state(k.base_block(t),t,T);value=k.excess(a,a)['connected'].real
  pred=lead*(c*t**6+(1-q)/4*dc*t**7)
  scaled=(value-pred)/t**8
  assert abs(scaled)<.1
  short.append(dict(q=q,t=t,T=T,connected=value,predicted_through_t7=pred,
                    C2=c,C2_derivative=dc,residual_over_t8=scaled,ratio_to_leading=value/(lead*c*t**6)))
 print(json.dumps({'q':q,'diagonal_rows':len(diagonal),'uncentered_min_eigenvalue':float(eig.min()),
                   'connected_excess_min_eigenvalue':float(np.linalg.eigvalsh(D).min())}),flush=True)

inputs=['code/bilocal_preparation_gram.py','code/check_bilocal_preparation_gram.py',
        'code/neutral_tower_projection.py','code/neutral_selected_response.py',
        'results/susskind-selected-bilocal-mean/time-checks.json']
result=dict(status='passed',diagonal=diagonal,matrices=matrices,short_time=short,previous_means=old_mean,
 parameters=dict(coarse=[24,48,24],fine=[32,64,32],scalar_Taylor_cutoff=80,physical_lambda='-log(q)/4'),
 environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,dtype='float64/complex128'),
 scope='Preparation contribution only. The common-center variance difference is for the full bilocal because its stationary energy-tilted baseline cancels. Absolute total variances are not numerically evaluated here. Analytic degree bounds exclude quadrature and roundoff.',
 seconds=time.monotonic()-start,input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({'status':'passed','seconds':result['seconds'],'sha256':hashlib.sha256(OUT.read_bytes()).hexdigest()}),flush=True)
