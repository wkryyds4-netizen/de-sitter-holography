"""New neutral two-time algebra and dynamics diagnostics; write once."""
from pathlib import Path
from datetime import datetime,timezone
from fractions import Fraction as F
from itertools import product
from math import ceil,cosh,sqrt
import hashlib,json,platform
import numpy as np
import scipy
import neutral_selected_response as ns
from sector_bilocal import moment as wick
from slater_composite_transfer import moment as selected_moment

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/susskind-neutral-selected-response/response-checks.json'
if OUT.exists():raise FileExistsError(OUT)
exact=[];selected=[]
for q in (F(1,4),F(2,3)):
    count=0
    for a,b,c in product(range(9),repeat=3):
        if a+b+c>8:continue
        word=(1,1)+(0,)*a+(1,)+(0,)*b+(1,)+(0,)*c
        connected=wick(word,q)-wick(word[2:],q)
        expected=(1+q)*ns.chain_moment(a,b,c,q)
        assert connected==expected,(a,b,c,q,connected,expected)
        count+=1
    exact.append(dict(q=str(q),three_block_identities=count))
for gamma,lam in [(F(1,2),F(1,3)),(F(2,3),F(1,5))]:
    q=gamma**4;s=gamma**2;count=0
    for a,b,c in product(range(7),repeat=3):
        if a+b+c>6:continue
        word='H'*a+'W'+'H'*b+'W'+'H'*c
        actual=selected_moment('-'+word+'+',gamma,lam)/lam
        expected=((1-s)*ns.pair_moment(a+c,b,q)+s*ns.pair_moment(a+c+2,b,q)
                  +s*lam*(1+q)*ns.chain_moment(a,b,c,q))
        assert actual==expected,(word,actual,expected)
        count+=1
    selected.append(dict(gamma=str(gamma),lambda_value=str(lam),identities=count))
print('Exact neutral block and selected transfer comparisons passed',flush=True)

pairs=[(0.,0.),(.25,0.),(.5,.5),(1.,1.),(.5,-.5),(.75,.25),(1.,-.5)]
independent=[]
for q in (.25,.64):
    eps=-np.log(q);sp=ns.SpectralKernel(eps,192)
    for u,v in pairs:
        b,_=sp.chain(u,v)
        first,mass1,norm1=ns.angular_chain(q,160,128,u,v)
        second,mass2,norm2=ns.angular_chain(q,224,192,u,v)
        assert max(abs(b-second),abs(first-second),abs(mass2-1),norm2)<5e-11
        independent.append(dict(q=q,u=u,v=v,chain=float(b.real),spectral_angular_difference=float(abs(b-second)),angular_refinement=float(abs(first-second)),mass=mass2,transition_normalization_error=norm2))
print('Independent angular-product comparisons passed',flush=True)

rows=[]
for eps in (.16,.08,.04,.02):
    L=max(128,ceil(12/eps)+64);low=ns.SpectralKernel(eps,L);high=ns.SpectralKernel(eps,L+64)
    points=[]
    for u,v in pairs:
        val,b,condition=high.selected(u,v);small,_,_=low.selected(u,v)
        trace=high.pair(u-v);pred=ns.predicted(u,v);diagerr=None
        if u==v:
            identity=1+sqrt(high.q)*eps/4*(1+high.q)*high.pair(u,d=2)
            diagerr=float(abs(val-identity));assert diagerr<2e-12
        if v==0:assert abs(b-high.pair(u))<2e-12
        assert max(abs(val.imag),abs(b.imag),abs(val-small))<2e-11
        points.append(dict(u=u,v=v,selected=float(val.real),trace=float(trace.real),chain=float(b.real),cutoff_difference=float(abs(val-small)),relative_cutoff_difference=float(abs(val-small)/max(abs(val),1e-300)),phase_sum_condition=condition,diagonal_identity_error=diagerr,
            response_coefficient=float((val.real-trace.real)/eps),predicted_response=pred['response'],
            selected_second_remainder=float((val.real-pred['trace_leading']-eps*pred['selected_first'])/eps**2),
            mixture_second_remainder=float(((val.real+trace.real)/2-pred['trace_leading']-eps*pred['mixture_first'])/eps**2),
            chain_first_scaled_error=float((b.real-cosh(u)**-2*cosh(v)**-2)/eps)))
    rows.append(dict(epsilon=eps,dimensions=[L,L+64],values=points))
    print('epsilon',eps,'passed',flush=True)

inputs=['code/neutral_selected_response.py','code/check_neutral_selected_response.py','code/sector_bilocal.py','code/slater_composite_transfer.py',
        'results/susskind-weighted-mehler-limit/local-evidence.json','results/susskind-sector-bilocal/local-evidence.json','results/susskind-slater-selected-state/local-evidence.json']
result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),exact=exact,selected_transfer=selected,independent=independent,rows=rows,
    parameters=dict(physical_epsilon=[.16,.08,.04,.02],times=pairs,angular_points=[160,224],product_terms=[128,192]),
    environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,dtype='float64/complex128'),
    input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},
    scope='Exact finite three-block moments and independent selected transfer checks; real-time spectral/angular comparisons and cutoff diagnostics, not interval certification. Uniform asymptotic remainder is proved separately in the manuscript.')
OUT.parent.mkdir(parents=True,exist_ok=True)
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',block_identities=sum(x['three_block_identities'] for x in exact),selected_identities=sum(x['identities'] for x in selected),sha256=hashlib.sha256(OUT.read_bytes()).hexdigest())))
