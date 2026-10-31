"""Two-endpoint asymptotics for the balanced KL 4F3; not a 6j limit proof.

All powers use the principal log. Existing finite-mass scientific kernels
are imported read-only. Diagnostics retain complex values and failures.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import platform
import sys
import mpmath
from observer_crossing import CrossingEvaluator

ROOT = Path(__file__).resolve().parents[1]

def b33_parameters(c, mass, q1):
    q = list(map(c.mpf, [q1, '.7', '1.1', '.4']))
    s1,s2,s3,s4 = [c.mpf(mass)+v for v in q]
    nu=c.mpf('.6'); h=c.mpf('.5')
    return (h-c.j*(s3+s4+nu), h-c.j*(s2+s3-nu),
            h-c.j*(s3-s4+nu), h+c.j*(s2-s3+nu),
            1-2*c.j*s3, 1-c.j*(s1+s3), 1+c.j*(s1-s3))

def leading(c, pars, W):
    A,C,b,d,E,F,g=pars
    delta=g-b-d
    pre=c.gammaprod([g],[b,d])
    if abs(delta)<c.mpf(10)**(-c.dps+8):
        return pre*(c.log(W)+c.digamma(1)-c.digamma(b)-c.digamma(d))
    lr=c.loggamma(1+delta)+c.loggamma(b)+c.loggamma(d)-c.loggamma(b+delta)-c.loggamma(d+delta)
    return pre*(c.expm1(lr)-c.expm1(-delta*c.log(W)))/delta

def exact(c, pars):
    A,C,b,d,E,F,g=pars
    return c.hyper([A,C,b,d],[E,F,g],1,strict=True,sum_method='r+s+e')

class EndpointEvaluator(CrossingEvaluator):
    def __init__(self, mass, mode, dps):
        super().__init__(dps,15)
        self.mass=mass; self.mode=mode; self.parameters=[]

    def hyper_f(self,h1,h2,h3,h4,h,k):
        c=self.ctx
        us=[1+k-h1-h4,k+h1-h4,2-h3-h4-h,1+h3-h4-h]
        ls=[2-2*h4,2+k-h2-h4-h,1+k+h2-h4-h]
        # This sorting is only used at the recorded bounded-offset fixtures.
        ui=sorted(range(4),key=lambda n:abs(c.im(us[n])),reverse=True)[:2]
        li=sorted(range(3),key=lambda n:abs(c.im(ls[n])),reverse=True)[:2]
        A,C=[us[n] for n in ui]; b,d=[v for n,v in enumerate(us) if n not in ui]
        E,F=[ls[n] for n in li]; g=[v for n,v in enumerate(ls) if n not in li][0]
        # Both E-A and F-C have real part 1/2; ordering within pairs is free.
        sign=c.sign(c.im(A)); W=2*c.j*self.mass*sign
        assert all(abs(c.im(v-W))<self.mass for v in [A,C,E,F])
        assert all(abs(c.im(v))<self.mass for v in [b,d,g])
        pars=(A,C,b,d,E,F,g); self.parameters.append(pars)
        if self.mode=='one':
            return c.gammaprod([g,g-b-d],[g-b,g-d])
        return leading(c,pars,W)

def main():
    rows=[]; checks=[]; precisions=[30,45]
    for q1 in ['.2','0','.00000001']:
        for mass in [4,16,64]:
            runs=[]
            for dps in precisions:
                c=mpmath.mp.clone(); c.dps=dps+15
                pars=b33_parameters(c,mass,q1)
                val=exact(c,pars); pred=leading(c,pars,-2*c.j*mass)
                delta=pars[-1]-pars[2]-pars[3]
                first=None if q1=='0' else c.gammaprod([pars[-1],delta],[pars[-1]-pars[2],pars[-1]-pars[3]])
                runs.append({'working_dps':c.dps,'exact':[str(c.re(val)),str(c.im(val))],
                    'two_branch':[str(c.re(pred)),str(c.im(pred))],
                    'two_branch_absolute_error':str(abs(val-pred)),
                    'one_branch_absolute_error':None if first is None else str(abs(val-first)),
                    'excess_error':str(abs(pars[4]+pars[5]+pars[6]-sum(pars[:4])-1))})
            c.dps=65
            vals=[c.mpc(*v['exact']) for v in runs]
            err=abs(vals[0]-vals[1])/max(1,abs(vals[1]))
            checks.append({'name':f'precision-q1={q1}-M={mass}', 'scaled_difference':str(err),'tolerance':'1e-25','passed':bool(err<c.mpf('1e-25'))})
            rows.append({'q':[''+q1,'.7','1.1','.4'],'Lambda':mass,'runs':runs})
            print(f'B33 q1={q1}, Lambda={mass}: saved',flush=True)
    # Connection formula at interior points, including close to the coalescence.
    for chi in ['.2','.00000001']:
        c=mpmath.mp.clone(); c.dps=60
        b=c.mpf('.5')-c.j*c.mpf('1.3'); d=c.mpf('.5')+c.j*c.mpf('.2')
        delta=c.j*c.mpf(chi); g=b+d+delta
        for z in ['.2','.9','.999']:
            z=c.mpf(z)
            left=c.hyp2f1(b,d,g,z)
            right=(c.gammaprod([g,delta],[g-b,g-d])*c.hyp2f1(b,d,1-delta,1-z)
                   +c.gammaprod([g,-delta],[b,d])*(1-z)**delta*c.hyp2f1(g-b,g-d,1+delta,1-z))
            err=abs(left-right)/max(1,abs(left))
            checks.append({'name':f'connection-chi={chi}-z={z}','scaled_difference':str(err),'tolerance':'1e-45','passed':bool(err<c.mpf('1e-45'))})
    four=[]
    for mass in [4,16,64]:
        c=mpmath.mp.clone(); c.dps=60
        labels=[c.mpf(mass)+c.mpf(v) for v in ['.2','.7','1.1','.4']]
        ev=CrossingEvaluator(40,15); val=ev.evaluate('.6',*labels)
        one=EndpointEvaluator(mass,'one',40).evaluate('.6',*labels)
        two=EndpointEvaluator(mass,'two',40).evaluate('.6',*labels)
        row={'Lambda':mass,'labels':[str(v) for v in labels], 'working_dps':55}
        for name,result in [('exact',val),('one',one),('two',two)]:
            row[name]={'raw_S':[str(c.re(result.raw_s)),str(c.im(result.raw_s))],
                       'four_terms':[[str(c.re(t)),str(c.im(t))] for t in result.g_terms],
                       'cancellation_ratio':str(result.cancellation_ratio)}
        row['one_relative_error']=str(abs(one.raw_s-val.raw_s)/abs(val.raw_s))
        row['two_relative_error']=str(abs(two.raw_s-val.raw_s)/abs(val.raw_s))
        row['added_branch_relative_size']=str(abs(two.raw_s-one.raw_s)/abs(val.raw_s))
        four.append(row); print(f'Four-branch Lambda={mass}: saved',flush=True)
    bound=['code/hyper_endpoint.py','code/observer_crossing.py','manuscript/endpoint.tex','sources/kolchmeyer-liu-2411.08090v1.html','sources/dlmf-15.8.html']
    result={'created_at_utc':datetime.now(timezone.utc).isoformat(),
        'kind':'finite-point endpoint asymptotic diagnostics; not exact 6j distribution convergence',
        'status':'passed' if all(x['passed'] for x in checks) else 'failed',
        'versions':{'python':sys.version,'mpmath':mpmath.__version__,'platform':platform.platform()},
        'parameters':{'nu':'.6','L':1,'Lambda':[4,16,64],'phase_branch':'principal Log of +/- 2i Lambda'},
        'single_hypergeometric_rows':rows,'connection_and_precision_checks':checks,'checks':len(checks),
        'four_branch_diagnostics':four,
        'limitations':['shared mpmath special-function implementation','empirical refinement, not interval error','four-branch diagnostic away from opposite-label coincidences','no mass integration or delta-test pairing','no numerical proof of uniform remainder'],
        'hashes':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in bound}}
    out=ROOT/'results/stage3/endpoint'; out.mkdir(parents=True,exist_ok=True)
    (out/'diagnostics.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(result['status'],len(checks),flush=True)
    if result['status']!='passed': raise SystemExit(1)

if __name__=='__main__': main()
