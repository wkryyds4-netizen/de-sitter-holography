"""Log-time coefficient/PV diagnostics, not exact 6j quadrature or certification."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, platform, sys
import mpmath
from crossing_assembly import parameters, u_closed, leading_raw

ROOT=Path(__file__).resolve().parents[1]

def pair_coefficients(c,offsets,nu):
    p,q,e,de=parameters(c,offsets,nu)
    up,um=u_closed(c,p,q,e,de)
    return 2**(-de)*up,2**de*um

def kl_b50_coefficient(c,off,nu):
    a,b,z,d=map(c.mpf,off);nu=c.mpf(nu);chi=a-b+z-d
    ap=a-b+nu;am=a-b-nu;half=c.mpf('.5')
    value=2*c.exp(-2*c.pi*nu)*(1+c.j*c.exp(c.pi*ap))*(1-c.j*c.exp(c.pi*(ap-chi)))
    value*=c.gamma(half-c.j*ap)*c.gamma(half+c.j*(ap-chi))
    value/=c.pi**2*c.gamma(half+c.j*am)*(1-c.j*c.exp(c.pi*am))*(1+c.j*c.exp(c.pi*(am-chi)))*c.gamma(half-c.j*(am-chi))
    return value*c.cosh(c.pi*chi/2)

def cj(c,z):return [str(c.re(z)),str(c.im(z))]

def main():
    c=mpmath.mp.clone();c.dps=60
    checks=[];identities=[]
    def check(name,error,tolerance):
        checks.append({'name':name,'absolute_error':str(error),'tolerance':str(tolerance),'passed':bool(error<c.mpf(tolerance))})
    for n,off in enumerate([['.2','.7','1.1','.4'],['-.3','.2','.5','1.2'],['1.2','-.1','.4','.3']]):
        p,q,e,de=parameters(c,off,'.6');vp,vm=pair_coefficients(c,off,'.6')
        ref=kl_b50_coefficient(c,off,'.6')
        check(f'B50-coefficient-{n}',abs(vp/c.pi**3-ref),'1e-50')
        L=c.mpf(16);T=c.mpf('.37');chi=c.im(de)
        rewritten=(c.exp(c.j*(2*c.log(2*L)+T)*chi)*c.gamma(de)*vp+c.exp(c.j*T*chi)*c.gamma(-de)*vm)/c.pi**3
        direct=L*c.exp(c.j*(c.log(2*L)+T)*chi)*leading_raw(c,L,off,'.6')
        check(f'positive-time-branch-{n}',abs(rewritten-direct),'1e-50')
        identities.append({'offsets':off,'chi':str(chi),'Vplus_over_pi3':cj(c,vp/c.pi**3),'B50':cj(c,ref),'dephased_leading':cj(c,direct)})

    # Compact transverse slice: a=.2, c=.4, d=.1, b=.5-chi.
    # Pair +h and -h before summing the PV; no pole is evaluated.
    c.dps=40;data={};orders=[40,80];Ts=[-2,0,2];logs=[2,4,8,16]
    def coeff(chi):
        off=[c.mpf('.2'),c.mpf('.5')-chi,c.mpf('.4'),c.mpf('.1')]
        vp,vm=pair_coefficients(c,off,'.6')
        bump=c.exp(1-1/(1-chi**2))*(1+c.mpf('.2')*chi+c.j*c.mpf('.1')*chi)
        return bump*c.gamma(c.j*chi)*vp/c.pi**3,bump*c.gamma(-c.j*chi)*vm/c.pi**3
    for order in orders:
        nodes,weights=c.gauss_quadrature(order,'legendre');cache=[]
        for left,right in [(c.zero,c.mpf('.5')),(c.mpf('.5'),c.mpf('.9')),(c.mpf('.9'),c.one)]:
            for node,w in zip(nodes,weights):
                h=(left+right)/2+(right-left)*node/2
                cache.append((h,w*(right-left)/2,coeff(h),coeff(-h)))
        limits={};finite=[]
        for T in Ts+[-16,16,-32,32]:
            pv=c.fsum(w*(cp[1]*c.exp(c.j*T*h)+cm[1]*c.exp(-c.j*T*h)) for h,w,cp,cm in cache)
            limits[T]=2/c.pi+pv
        for T in Ts:
            for logL in logs:
                k=2*(logL+c.log(2))+T
                value=c.fsum(w*(cp[0]*c.exp(c.j*k*h)+cm[0]*c.exp(-c.j*k*h)+cp[1]*c.exp(c.j*T*h)+cm[1]*c.exp(-c.j*T*h)) for h,w,cp,cm in cache)
                finite.append((T,logL,value))
        data[order]=(limits,finite)
        print('completed compact leading-expression quadrature',order,flush=True)
    lo,hi=orders;lim_lo,fin_lo=data[lo];lim_hi,fin_hi=data[hi]
    for T in sorted(lim_hi):check(f'PV-order-stability-T{T}',abs(lim_hi[T]-lim_lo[T]),'1e-12')
    rows=[]
    for (T,logL,val),(_,_,lower) in zip(fin_hi,fin_lo):
        check(f'leading-order-stability-T{T}-logL{logL}',abs(val-lower),'1e-12')
        rows.append({'T':T,'log_Lambda':logL,'leading_pairing':cj(c,val),'PV_delta_target':cj(c,lim_hi[T]),'absolute_difference':str(abs(val-lim_hi[T])),'quadrature_order_change':str(abs(val-lower))})
        print('leading',T,logL,'error',c.nstr(abs(val-lim_hi[T]),8),flush=True)
    plateaux=[{'T':T,'PV_delta_value':cj(c,lim_hi[T]),'ordered_plateau_target':str(4/c.pi if T<0 else 0),'absolute_difference':str(abs(lim_hi[T]-(4/c.pi if T<0 else 0)))} for T in [-16,16,-32,32]]
    hashes=['code/logtime_limit.py','code/crossing_assembly.py','manuscript/logtime.tex','manuscript/assembly.tex','manuscript/crosslimit.tex','sources/kolchmeyer-liu-2411.08090v1.html']
    result={'created_at_utc':datetime.now(timezone.utc).isoformat(),'status':'passed' if all(r['passed'] for r in checks) else 'failed',
        'kind':'finite-point coefficient and compact leading-expression PV diagnostics; no exact 6j distribution integration',
        'checks':checks,'coefficient_rows':identities,'leading_distribution_rows':rows,'plateau_rows':plateaux,
        'parameters':{'nu':'.6','coefficient_dps':60,'quadrature_dps':40,'Gauss_Legendre_orders_per_subinterval':orders,'positive_subintervals':['0','.5','.9','1'],'slice':'a=.2,c=.4,d=.1,b=.5-chi','test':'exp(1-1/(1-chi^2))*(1+.2*chi+.1*i*chi) on (-1,1)','physical_time':'log(2 Lambda)+T'},
        'versions':{'python':sys.version,'mpmath':mpmath.__version__,'platform':platform.platform()},
        'limitations':['Gauss-order agreement is empirical, not an interval error bound','Only compact leading-expression slice, not full physical packet or exact 6j quadrature','Finite-point identities supplement analytic derivation, not independent proof certification','No finite-G result or exponent extracted'],
        'hashes':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in hashes}}
    dest=ROOT/'results/stage3/logtime';dest.mkdir(parents=True,exist_ok=True)
    (dest/'diagnostics.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(result['status'],len(checks),'checks',flush=True)
    if result['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
