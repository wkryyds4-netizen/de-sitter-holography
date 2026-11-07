"""Fresh odd-source response, using only existing pure helpers."""
from pathlib import Path
from datetime import datetime,timezone
from fractions import Fraction as F
from math import cosh,sqrt
import hashlib,json,platform
import numpy as np
import scipy
from neutral_selected_response import SpectralKernel,pair_moment
from slater_composite_transfer import moment
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/susskind-neutral-selected-response/mixed-source-checks.json'
if OUT.exists():raise FileExistsError(OUT)
exact=[]
for gamma,lam in [(F(1,2),F(1,3)),(F(2,3),F(1,5))]:
    q=gamma**4;count=0
    for a in range(9):
        for b in range(9-a):
            word1='H'*(a+1)+'W'+'H'*b
            word2='H'*a+'W'+'H'*(b+1)
            actual=moment('-'+word1+'+',gamma,lam)+moment('-'+word2+'+',gamma,lam)
            expected=-gamma**2*lam*(pair_moment(a+2,b,q)+2*pair_moment(a+1,b+1,q)+pair_moment(a,b+2,q))
            assert actual==expected,(a,b,actual,expected)
            count+=1
    exact.append(dict(gamma=str(gamma),lambda_value=str(lam),identities=count))
parent='results/susskind-neutral-selected-response/response-checks.json'
assert hashlib.sha256((ROOT/parent).read_bytes()).hexdigest()=='475535bab4df2947c1de0d886ae42672e478becd612ffb66338afc14b09af843'
old=json.loads((ROOT/parent).read_text());rows=[]
for row in old['rows']:
    eps=row['epsilon'];points=[]
    models=[SpectralKernel(eps,L) for L in row['dimensions']]
    for t in (0.,.25,.5,1.,1.5):
        values=[]
        for model in models:
            a=model.V@(model.root*np.exp(1j*t*model.E/sqrt(eps)))
            b=model.V@(model.root*model.E*np.exp(1j*t*model.E/sqrt(eps)))
            c=model.V@(model.root*model.E**2*np.exp(1j*t*model.E/sqrt(eps)))
            terms=2*np.real(np.conj(a)*model.D*c)+2*abs(b)**2*model.D
            values.append((float(np.sum(terms)),float(np.sum(abs(terms))/max(abs(np.sum(terms)),1e-300))))
        low,H=values[0][0],values[1][0];q=models[-1].q
        assert abs(low-H)<3e-11
        if t==0:assert abs(H-2*(1+q))<2e-12
        mean=-sqrt(q*eps/4)*H
        points.append(dict(u=t,H=H,selected_mean=mean,scaled_mean=mean/sqrt(eps),leading_scaled_mean=-2/cosh(t)**2,
            first_scaled_remainder=(mean/sqrt(eps)+2/cosh(t)**2)/eps,cutoff_difference=abs(low-H),phase_condition=values[-1][1]))
    rows.append(dict(epsilon=eps,dimensions=row['dimensions'],values=points))
inputs=['code/check_neutral_mixed_source.py','code/neutral_selected_response.py','code/slater_composite_transfer.py',parent]
result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),exact=exact,rows=rows,
    environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
    input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},
    scope='Exact selected transfer for the anticommutator with W; numerical mean for A={h,U}. The compact-time result follows from the weighted pair-density proof, not this finite grid.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',exact_identities=sum(x['identities'] for x in exact),times=20,max_cutoff_difference=max(x['cutoff_difference'] for r in rows for x in r['values']),sha256=hashlib.sha256(OUT.read_bytes()).hexdigest())))
