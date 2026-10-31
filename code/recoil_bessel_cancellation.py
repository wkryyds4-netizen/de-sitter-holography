"""Second integrated coefficient, auxiliary zeros, and exact-kernel lifting.

The joint growing-K remainder remains open; coefficient asymptotics are separate.
"""
from pathlib import Path
from datetime import datetime, timezone
import sys,platform,json,hashlib,re
import mpmath as mp
import sympy as sp
from recoil_bessel_limit import J,log_q,auxiliary,complete,recurrence
ROOT=Path(__file__).resolve().parents[1]

def parse_complex(s):
    vals=re.findall(r'[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?',s.replace(' ',''))
    if len(vals)!=2:raise ValueError(s)
    return mp.mpc(mp.mpf(vals[0]),mp.mpf(vals[1]))

def second(K,alpha):
    S=K*K+alpha*alpha-1
    dj=K*mp.diff(lambda k:J(k,alpha),K)/2
    return (2*S*dj+(K*K-alpha*alpha+1)*J(K,alpha))/3

def root(alpha):
    return mp.findroot(lambda k:mp.log(mp.besselk(alpha,k))-mp.log(mp.pi)
                       -mp.log(mp.besseli(alpha,k)),(alpha/2,alpha))

def run():
    mp.mp.dps=75
    checks=[];zeros=[];integrals=[];compact=[];comparisons=[];fixed_rows=[];joint_tails=[]
    def check(name,error,tol):
        checks.append(dict(name=name,error=str(abs(error)),tolerance=str(tol),
                           passed=bool(abs(error)<=tol)))
    k,a,t=sp.symbols('k a t',real=True)
    j,dj=sp.symbols('J DJ');S=k*k+a*a-1;q=S/4
    def D(expr):
        return sp.expand(k*sp.diff(expr,k)/2+dj*sp.diff(expr,j)+(dj+q*j)*sp.diff(expr,dj))
    powers=[j]
    for n in range(6):powers.append(D(powers[-1]))
    log2=sp.Rational(8,3)*sp.bernoulli(3,t)+4*k*k*t
    poly=sp.expand(log2-8*(t*t-t)**2)
    raw=sum(poly.coeff(t,n)*powers[n] for n in range(5))
    phase_removed=sp.expand(raw+S*S*j/2)
    target=(2*S*dj+(k*k-a*a+1)*j)/3
    check('symbolic-second-operator-reduction',int(sp.simplify(phase_removed-target)!=0),0)
    check('symbolic-first-operator',int(sp.simplify(4*(powers[2]-powers[1])-S*j)!=0),0)
    L1=4*sp.I*(t*t-t)
    L3=4*sp.I*k*k*t*t-sp.Rational(8,3)*sp.I*(t*t-t)**2
    poly3=sp.expand(L3+L1*log2+L1**3/6)
    raw3=sum(poly3.coeff(t,n)*powers[n] for n in range(7))
    P=(9*k**4-(a*a-1)**2)/6
    third_removed=sp.simplify(raw3-sp.I*S*raw-sp.I*S**3*j/3)
    check('symbolic-third-pure-phase',int(sp.simplify(third_removed-sp.I*P*j)!=0),0)
    Y=mp.mpf(2);D=max(mp.mpf(2),Y*Y/2);cY=8+2*Y
    for K in map(mp.mpf,[8,16,32]):
        alpha=K*Y;A=max(mp.mpf(2),2*K,2*max(alpha*alpha/4-1,mp.mpf(0)));R=cY*K
        for nu in [mp.sqrt(32*(D+1)*K**3),4*K**3]:
            for x in [R,(R+A+1)/2,A+1]:
                value=abs(recurrence(x+1j,K,nu,alpha))
                check(f'joint-tail-contraction-{K}-{nu}-{x}',max(mp.mpf(0),8*value-1),0)
                joint_tails.append(dict(K=str(K),alpha=str(alpha),nu=str(nu),R=str(R),
                                       A=str(A),x=str(x),modulus=str(value),bound='0.125'))
    for kk,aa in [(2,'.7'),(4,8)]:
        K,alpha=mp.mpf(kk),mp.mpf(aa)
        last=None
        for nu in map(mp.mpf,[128,256,512]):
            tt=mp.mpc('2.3','1')
            Q=mp.exp(log_q(tt,K,nu,alpha));U=auxiliary(tt,K,alpha)
            L1=4j*(tt*tt-tt)
            L2=mp.mpf(8)/3*(tt**3-mp.mpf('1.5')*tt**2+tt/2)+4*K*K*tt
            pred=mp.exp(L1/nu+L2/nu**2)
            err=abs(Q/U/pred-1)
            third=4j*K*K*tt*tt-mp.mpf(8)/3*1j*(tt*tt-tt)**2
            err3=abs(Q/U/mp.exp(L1/nu+L2/nu**2+third/nu**3)-1)
            compact.append(dict(K=str(K),alpha=str(alpha),nu=str(nu),
                                remainder=str(err),scaled_nu3=str(nu**3*err),
                                third_log_remainder=str(err3),scaled_nu4=str(nu**4*err3)))
            if last is not None:check(f'compact-third-order-decrease-{K}-{alpha}-{nu}',max(mp.mpf(0),err/last-mp.mpf('.2')),0)
            last=err
    yc=mp.findroot(lambda y:mp.sqrt(1+y*y)-y*mp.asinh(y),(1,2))
    for alpha in map(mp.mpf,['1.5','3.5','7.5','15.5','31.5','63.5']):
        K=root(alpha);I=mp.besseli(alpha,K);T=second(K,alpha)
        scale=K*(mp.besselk(alpha,K)+mp.pi*I)
        check('auxiliary-zero-'+str(alpha),abs(J(K,alpha))/scale,mp.mpf('1e-65'))
        derivative=mp.diff(lambda k:J(k,alpha),K)
        check('simple-zero-Wronskian-'+str(alpha),derivative/(-1j/I)-1,mp.mpf('1e-60'))
        explicit=-1j*K*(K*K+alpha*alpha-1)/(3*I)
        check('zero-lifting-coefficient-'+str(alpha),T/explicit-1,mp.mpf('1e-60'))
        zeros.append(dict(alpha=str(alpha),K=str(K),y=str(alpha/K),
                          scaled_y_gap=str(alpha*alpha*(alpha/K-yc)),
                          derivative=str(derivative),second_coefficient=str(T),branch_envelope=str(scale)))
    for zr in zeros[:2]:
        alpha=mp.mpf(zr['alpha']);K=mp.mpf(zr['K']);T=second(K,alpha)
        last=None
        for nu in [64,128,256]:
            low,_=complete(K,nu,alpha,20,60)
            value,row=complete(K,nu,alpha,32,80)
            check(f'zero-quadrature-{alpha}-{nu}',abs(low/value-1),mp.mpf('1e-14'))
            row['quadrature_change']=str(abs(low/value-1))
            if nu==256:
                high,_=complete(K,nu,alpha,32,100)
                ext,_=complete(K,nu,alpha,32,80,X=44,T=24,eta=mp.mpf('.75'))
                check(f'zero-precision-{alpha}',abs(high/value-1),mp.mpf('1e-50'))
                check(f'zero-contour-extent-{alpha}',abs(ext/value-1),mp.mpf('1e-16'))
                row['precision_change']=str(abs(high/value-1));row['contour_extent_change']=str(abs(ext/value-1))
            phase=mp.exp(1j*(K*K+alpha*alpha-1)/nu)
            prediction=phase*T/nu**2
            err=abs(value/prediction-1)
            row.update(auxiliary_zero_residual=str(J(K,alpha)),second_coefficient=str(T),
                       zero_lifting_prediction=str(prediction),relative_prediction_error=str(err),
                       scaled_nu_prediction_error=str(nu*err))
            if last is not None:check(f'zero-lifting-decrease-{alpha}-{nu}',max(mp.mpf(0),err/last-mp.mpf('.65')),0)
            last=err;integrals.append(row)
            print(f'completed zero alpha={alpha}, nu={nu}',flush=True)
    # Reuse frozen complete integrals; do not modify or rerun them.
    old=json.loads((ROOT/'results/stage4/recoil-bessel-limit/diagnostics.json').read_text())
    for row in old['integral_rows']:
        K,alpha,nu=map(mp.mpf,(row['K'],row['alpha'],row['nu']))
        actual=parse_complex(row['value']);base=J(K,alpha)
        phase=mp.exp(1j*(K*K+alpha*alpha-1)/nu)
        pred=phase*(base+second(K,alpha)/nu**2)
        fixed_rows.append(dict(K=str(K),alpha=str(alpha),nu=str(nu),
                               second_absolute_error=str(abs(actual-pred)),
                               scaled_nu3_absolute_error=str(nu**3*abs(actual-pred))))
    previous=json.loads((ROOT/'results/stage4/recoil-global/diagnostics.json').read_text())
    for row in previous['bessel_comparisons']:
        K,nu,y=map(mp.mpf,(row['K'],row['nu'],row['y']));alpha=K*y
        actual=parse_complex(row['exact_graph_chiral_ratio'])
        base=J(K,alpha)/J(K,0)*mp.exp(1j*alpha*alpha/nu)
        coeff=second(K,alpha)/J(K,alpha)-second(K,0)/J(K,0)
        pred=base*(1+coeff/nu**2)
        err=abs(actual/pred-1);before=abs(actual/base-1)
        check(f'joint-sample-second-improves-{K}-{y}',max(mp.mpf(0),err/before-1),0)
        corrected=(mp.exp(1j*(alpha*alpha/nu-(alpha**4-2*alpha*alpha)/(6*nu**3)))
                   *(J(K,alpha)+second(K,alpha)/nu**2)/(J(K,0)+second(K,0)/nu**2))
        higher_err=abs(actual/corrected-1)
        check(f'joint-sample-cubic-phase-improves-{K}-{y}',max(mp.mpf(0),higher_err/err-1),0)
        comparisons.append(dict(K=str(K),nu=str(nu),y=str(y),phase_only_error=str(before),
                                 phase_and_second_error=str(err),coefficient=str(coeff),
                                 cubic_phase_rational_error=str(higher_err),
                                 coefficient_over_K3=str(coeff/K**3),
                                 predicted_large_K_coefficient=str(((1+y*y)**mp.mpf('1.5')-1)/3
                                     if y<yc else -((1+y*y)**mp.mpf('1.5')+1)/3)))
    inputs=['code/recoil_bessel_cancellation.py','code/recoil_bessel_limit.py',
            'results/stage4/recoil-bessel-limit/diagnostics.json','results/stage4/recoil-global/diagnostics.json',
            'sources/dlmf-5.11.html']
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
              status='passed' if all(c['passed'] for c in checks) else 'failed',
              scope='second fixed-parameter coefficient and auxiliary-zero lifting; joint remainder not proved',
              environment=dict(python=sys.version,platform=platform.platform(),
                               mpmath=mp.__version__,sympy=sp.__version__,dps=[60,75,80,100],orders=[20,32]),
              symbolic=dict(log_second=str(log2),kernel_second=str(poly),phase_removed=str(sp.factor(phase_removed)),
                            log_third=str(L3),third_pure_phase=str(sp.factor(third_removed))),
              checks=checks,zero_rows=zeros,zero_integral_rows=integrals,compact_rows=compact,
              fixed_second_rows=fixed_rows,joint_coefficient_comparisons=comparisons,
              joint_tail_rows=joint_tails,
              critical_y=str(yc),hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    out=ROOT/'results/stage4/recoil-bessel-cancellation';out.mkdir(parents=True,exist_ok=True)
    (out/'diagnostics.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=data['status'],checks=len(checks),failed=[c['name'] for c in checks if not c['passed']])))
    if data['status']!='passed':raise SystemExit(1)

if __name__=='__main__':run()
