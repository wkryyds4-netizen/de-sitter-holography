"""Analytic four-branch assembly diagnostics; no interval certification.

Keeps all existing exact kernels read-only. Retains complex values and
distinguishes finite-point identities from empirical asymptotic errors.
"""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,platform,sys
import mpmath
from observer_crossing import CrossingEvaluator
from observer_crossing_regular import RegularCrossingEvaluator
ROOT=Path(__file__).resolve().parents[1]

def cfun(c,z): return c.sqrt(c.pi)*c.gammaprod([(1-z)/2],[z/2])
def lc(c,z): return -(c.digamma((1-z)/2)+c.digamma(z/2))/2
def parameters(c,offsets,nu):
    a,b,z,d=[c.mpf(x) for x in offsets];nu=c.mpf(nu);h=c.mpf('.5')
    return h+c.j*(d-b-2*nu)/2,h+c.j*(d-b+2*nu)/2,c.j*(a-z),c.j*(a-b+z-d)
def u_closed(c,p,q,e,delta):
    t=(delta-e)/2;v=(delta+e)/2;cf=lambda z:cfun(c,z)
    up=2**(1-delta)*c.pi*c.cos(c.pi*delta/2)*cf(p+v)*cf(q+v)/cf(p-t)/cf(q-t)
    um=2**(1+delta)*c.pi*c.cos(c.pi*delta/2)*cf(p-v)*cf(q-v)/cf(p+t)/cf(q+t)
    return up,um
def coefficients(c,p,q,e,delta):
    t=(delta-e)/2;v=(delta+e)/2;cf=lambda z:cfun(c,z)
    D=cf(1-p-t)*cf(q+v)*cf(1-q+t)*cf(p-v)
    ps=[D*cf(p+v)/cf(p-v)/cf(e),D*cf(1-q-t)/cf(1-q+t)/cf(e),D*cf(q-v)/cf(q+v)/cf(-e),D*cf(1-p+t)/cf(1-p-t)/cf(-e)]
    hs=[(1-p-v,p-v,1-e,delta,-1),(1-q+t,q+t,1-e,-delta,1),(1-q+v,q+v,1+e,-delta,1),(1-p-t,p-t,1+e,delta,-1)]
    regular=[P*c.gammaprod([g],[g-b,g-d]) for P,(b,d,g,_,_) in zip(ps,hs)]
    extra=[-P*c.gammaprod([g],[b,d])*(2*sign*c.j)**(-dd)/dd for P,(b,d,g,dd,sign) in zip(ps,hs)]
    return ps,hs,regular,extra
def leading_raw(c,L,offsets,nu):
    p,q,e,de=parameters(c,offsets,nu);L=c.mpf(L)
    if abs(de)<c.mpf(10)**(-c.dps+8):
        return 4/(c.pi**2*L)*(c.log(L/2)-c.euler+c.fsum(lc(c,z+sgn*e/2) for z in [p,q] for sgn in [-1,1])/2)
    up,um=u_closed(c,p,q,e,de)
    return (L**de*c.gamma(de)*up+L**(-de)*c.gamma(-de)*um)/(c.pi**3*L)
def complex_json(c,z):return [str(c.re(z)),str(c.im(z))]

def main():
    c=mpmath.mp.clone();c.dps=75;checks=[];identities=[]
    fixtures=[['.2','.7','1.1','.4'],['-.3','.2','.5','1.2'],['1.2','-.1','.4','.3']]
    def check(name,err,tol='1e-65'):
        checks.append({'name':name,'error':str(err),'tolerance':tol,'passed':bool(err<c.mpf(tol))})
    for n,off in enumerate(fixtures):
        p,q,e,de=parameters(c,off,'.6');ps,hs,rr,xx=coefficients(c,p,q,e,de);up,um=u_closed(c,p,q,e,de)
        check(f'extra-cancellation-{n}',abs(c.fsum(xx))/max(1,c.fsum(abs(z) for z in xx)))
        check(f'Uplus-{n}',abs(up-rr[0]-rr[3])/max(1,abs(up)))
        check(f'Uminus-{n}',abs(um-rr[1]-rr[2])/max(1,abs(um)))
        zp,zm=u_closed(c,p,q,e,c.zero)
        check(f'residue-plus-{n}',abs(zp-2*c.pi));check(f'residue-minus-{n}',abs(zm-2*c.pi))
        check(f'c-complement-{n}',abs(cfun(c,p)*cfun(c,1-p)-c.pi))
        identities.append({'offsets':off,'P0':[complex_json(c,z) for z in ps],'extra':[complex_json(c,z) for z in xx],'Uplus':complex_json(c,up),'Uminus':complex_json(c,um)})
    # Small complex e,delta displacement tests the algebra used in pole removal.
    p,q,e,de=parameters(c,fixtures[0],'.6');e+=c.mpf('.03');de+=c.mpf('.02')
    _,_,rr,xx=coefficients(c,p,q,e,de);up,um=u_closed(c,p,q,e,de)
    check('complex-tube-extra',abs(c.fsum(xx))/c.fsum(abs(z) for z in xx))
    check('complex-tube-Uplus',abs(up-rr[0]-rr[3])/max(1,abs(up)))
    check('complex-tube-Uminus',abs(um-rr[1]-rr[2])/max(1,abs(um)))
    exact_rows=[]
    cases={'generic':['.2','.7','1.1','.4'], 'delta_zero':['0','.7','1.1','.4'], 'opposite_equal':['.3','.7','.3','.4'], 'all_equal':['.2','.2','.2','.2']}
    for name,off in cases.items():
        for L in [4,16]:
            ev=RegularCrossingEvaluator(dps=25,guard_digits=15,max_working_dps=120)
            cc=ev.ctx;labs=[cc.mpf(L)+cc.mpf(z) for z in off]
            result=ev.evaluate('.6',*labs)
            # raw S is chart-dependent: reconstruct in the ORIGINAL chart.
            raw=result.symmetric_s*(4096*cc.pi**11)/ev.gamma_product(labs,cc.mpf('.6'))
            pred=leading_raw(cc,L,off,'.6')
            exact_rows.append({'case':name,'offsets':off,'Lambda':L,'working_dps':result.working_dps,'method':result.method,'chart_rotation':result.chart_rotation,
                'raw_S':complex_json(cc,raw),'assembled_leading':complex_json(cc,pred),
                'absolute_error':str(abs(raw-pred)),'Lambda_squared_error':str(L**2*abs(raw-pred)),
                'exact_cancellation_ratio':str(result.cancellation_ratio),'exact_reality_relative_error':str(result.reality_relative_error)})
            print(name,L,'error',cc.nstr(abs(raw-pred),6),flush=True)
    # Tests of the explicit leading distribution only, not exact-kernel quadrature.
    tests=[]
    c.dps=35
    for logL in [2,4,8,16]:
        L=c.exp(logL)
        def integrand(chi):
            if abs(chi)==1:return c.zero
            bump=c.exp(1-1/(1-chi**2))
            off=[chi+c.mpf('.2'),c.mpf('.2'),c.mpf('.2'),c.mpf('.2')]
            return L*leading_raw(c,L,off,'.6')*bump
        value=c.quad(integrand,[-1,-.5,0,.5,1])
        tests.append({'log_Lambda':logL,'value':complex_json(c,value),'delta_target':str(4/c.pi),'absolute_difference':str(abs(value-4/c.pi))})
        print('leading compact test',logL,c.nstr(value,9),flush=True)
    bound=['code/crossing_assembly.py','code/observer_crossing.py','code/observer_crossing_regular.py','manuscript/assembly.tex','manuscript/endpoint.tex','sources/kolchmeyer-liu-2411.08090v1.html']
    out={'created_at_utc':datetime.now(timezone.utc).isoformat(),'kind':'finite-point algebra and assembled asymptotic diagnostics; no interval certification',
         'status':'passed' if all(x['passed'] for x in checks) else 'failed','checks':checks,'identity_rows':identities,'exact_kernel_rows':exact_rows,'leading_distribution_tests':tests,
         'parameters':{'nu':'.6','identity_working_dps':75,'exact_requested_dps':25,'exact_guard_digits':15,'leading_integral_dps':35,'bump':'exp(1-1/(1-chi^2)) on (-1,1), zero outside'},
         'versions':{'python':sys.version,'mpmath':mpmath.__version__,'platform':platform.platform()},
         'limitations':['mpmath is shared by algebra/precision computations','finite-point errors do not certify the uniform theorem','compact test uses assembled leading expression, not exact 6j integral','no shrinking layer or growing-time numerical certification'],
         'hashes':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in bound}}
    dest=ROOT/'results/stage3/assembly';dest.mkdir(parents=True,exist_ok=True)
    (dest/'diagnostics.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(out['status'],len(checks),'algebra checks',flush=True)
    if out['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
