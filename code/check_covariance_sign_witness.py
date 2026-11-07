"""Exact witness that a centered preparation correction need not be positive."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,platform
from selected_bilocal_covariance import JointModel,centered_prefix,field_pair,tilted_pair
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/susskind-selected-bilocal-covariance/sign-witness.json'
if OUT.exists():raise FileExistsError(OUT)
gamma=F(2,3);q=gamma**4;lam=F(1,5);m=JointModel(q,lam)
X,_=centered_prefix(2,2,m)
mean=field_pair(m.H(0),X,0,0,gamma,lam,'selected')
second=field_pair(X,X,0,0,gamma,lam,'selected')
base=tilted_pair(X,X,0,0,gamma,lam)
delta=second-base-mean*mean
formula=gamma**2*lam**2*(1+q)**2*(1+q-4*gamma**2)
assert mean==-2*gamma**2*lam*(1+q)
assert delta==formula and delta<0
assert second-mean*mean>0
inputs=['code/check_covariance_sign_witness.py','code/selected_bilocal_covariance.py',
        'code/slater_higher_wick_sources.py','code/neutral_tower_projection.py']
result=dict(status='passed',q=str(q),lambda_=str(lam),mean=str(mean),selected_second=str(second),
 energy_tilted_variance=str(base),selected_variance=str(second-mean*mean),centered_excess=str(delta),formula=str(formula),
 analytic_sign='For the physical curve, lambda>0 and q>(2-sqrt(3))^2 imply zeta*lambda^2*(1+q)^2*(1+q-4*zeta)<0. The exact rational example checks algebra, not the exponential parameter relation.',
 python=platform.python_version(),input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({'status':'passed','centered_excess':str(delta),'selected_variance':str(second-mean*mean),'sha256':hashlib.sha256(OUT.read_bytes()).hexdigest()}))
