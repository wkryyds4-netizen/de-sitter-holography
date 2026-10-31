"""Full deformed Barnes quadrature and a two-branch Bessel comparison.

Finite numerical evidence is not a uniform global asymptotic theorem.
"""
from pathlib import Path
from datetime import datetime, timezone
import sys, platform, json, hashlib
import mpmath as mp
from recoil_saddle import labels,saddle,phase,phase_second
from recoil_local import log_prefactor
from recoil_internal import barnes
ROOT=Path(__file__).resolve().parents[1]
NODES={}

def left_log(w,a,r):
    return (mp.log(8*mp.pi**3)+4j*mp.pi*w+mp.pi*sum(r)
            +mp.loggamma(w+mp.mpf('.5'))
            +sum(mp.loggamma(w+mp.mpf('.5')-1j*t) for t in r)
            -sum(mp.loggamma(1+w-1j*t) for t in a))

def left_factor(w,a):
    return (1+mp.exp(2j*mp.pi*w))/mp.fprod(1-mp.exp(2j*mp.pi*(w-1j*t)) for t in a)

def right_log(z,a,r):
    return (-1j*mp.pi*z+sum(mp.loggamma(z+1j*t) for t in a)
            +sum(mp.loggamma(mp.mpf('.5')-z-1j*t) for t in r)
            -mp.loggamma(z+mp.mpf('.5')))

def complete(K,nu,y,order=32,dps=60,line_gap='1',Xfactor=8,T=18):
    with mp.workdps(dps):
        K,nu,y,gap=map(mp.mpf,(K,nu,y,line_gap))
        M=nu*nu/(2*K);a,r=labels(M,nu,1j*K*y)
        z,D=saddle(a,r);L=mp.re(z);c=M+nu-gap
        if not (max(mp.re(t) for t in a)<c<min(mp.re(t) for t in r) and L>K*abs(y)/2):
            raise ValueError('Deformed contour geometry not valid.')
        base=phase(z,a,r)+log_prefactor(z,a)
        def fn(zz):return mp.exp(right_log(zz,a,r)-base)
        def vl(t):
            w=L+1j*(c+t)
            return mp.exp(left_log(w,a,r)-base)*left_factor(w,a)
        key=(order,dps)
        if key not in NODES:NODES[key]=mp.gauss_quadrature(order,'legendre')
        nodes,weights=NODES[key]
        def integ(f,cuts):
            ans=0
            for lo,hi in zip(cuts[:-1],cuts[1:]):
                for nd,wt in zip(nodes,weights):
                    ans+=wt*(hi-lo)/2*f((lo+hi)/2+nd*(hi-lo)/2)
            return ans
        X=max(mp.mpf(32),Xfactor*K)
        negcuts=sorted(set([-L,mp.mpf(0)]+[mp.mpf(i) for i in range(int(mp.ceil(-L)),0)]))
        poscuts=sorted(set([mp.mpf(0),X,L]+[mp.mpf(i) for i in range(1,int(X))]))
        neg=integ(lambda x:fn(x-1j*c),negcuts)
        pos=integ(lambda x:fn(x-1j*c),poscuts)
        vc=[mp.mpf(i) for i in [0,1,2,4,8,12,T]]
        vertical=integ(vl,vc)
        total=pos+neg+1j*vertical
        gauss=mp.sqrt(2*mp.pi/(-phase_second(z,a,r)))
        return dict(K=str(K),nu=str(nu),M=str(M),y=str(y),order=order,dps=dps,line_gap=str(gap),
                    X=str(X),T=T,saddle=str(z),base=str(base),scaled_total=str(total),
                    scaled_positive=str(pos),scaled_negative=str(neg+1j*vertical),
                    scaled_vertical=str(vertical),local_gaussian=str(gauss),
                    full_over_local=str(total/gauss),
                    segment_cancellation=str((abs(pos)+abs(neg)+abs(vertical))/abs(total)),
                    _total=total,_base=base)

def bessel_branch(K,y):
    alpha=K*y
    return mp.exp(1j*mp.pi*alpha)*mp.besselk(alpha,K)+1j*mp.pi*mp.besseli(alpha,K)

def run():
    mp.mp.dps=65
    checks=[];identities=[];rows=[];ratios=[]
    def check(name,error,tol):
        checks.append(dict(name=name,error=str(abs(error)),tolerance=str(tol),passed=bool(abs(error)<=tol)))
    for yy in ['0','1','2']:
        K,nu=mp.mpf(4),mp.mpf(80);y=mp.mpf(yy)
        a,r=labels(nu*nu/(2*K),nu,1j*K*y);z,_=saddle(a,r)
        for fac in [1,2]:
            w=fac*z.real-1j*z.imag
            ex=mp.exp(left_log(w,a,r))*left_factor(w,a)
            check('negative-reflection-'+yy+'-'+str(fac),ex/barnes(-w,a,r)-1,mp.mpf('1e-50'))
            f=lambda t:t*(mp.log(t)-1)
            psi=4j*mp.pi*w+mp.pi*sum(r)+f(w)+sum(f(w-1j*t) for t in r)-sum(f(w-1j*t) for t in a)
            check('opposite-exponent-'+yy+'-'+str(fac),psi-(2j*mp.pi*w-mp.conj(phase(mp.conj(w),a,r))),mp.mpf('1e-50'))
        gap=2*mp.pi*(-z.imag)+2*mp.re(phase(z,a,r))
        identities.append(dict(y=yy,K=str(K),nu=str(nu),exponent_gap=str(gap),
                               limiting_gap_over_K=str(2*(mp.sqrt(1+y*y)-y*mp.asinh(y)))))
    critical=mp.findroot(lambda y:mp.sqrt(1+y*y)-y*mp.asinh(y),(1,2))
    # Exact auxiliary residue sum equals the continued K branch; noninteger orders.
    for alpha in map(mp.mpf,['.3','1.2','2.7']):
        K=mp.mpf(3)
        J=-(mp.exp(-1j*mp.pi*alpha)*mp.gamma(-alpha)*(K/2)**(1+alpha)*mp.hyper([], [1+alpha],K*K/4)
            +mp.exp(1j*mp.pi*alpha)*mp.gamma(alpha)*(K/2)**(1-alpha)*mp.hyper([], [1-alpha],K*K/4))
        check('auxiliary-Bessel-connection-'+str(alpha),J/(-K*bessel_branch(K,alpha/K))-1,mp.mpf('1e-50'))
        jf=lambda k:-k*bessel_branch(k,alpha/k)
        differential=(K*K*mp.diff(jf,K,2)-K*mp.diff(jf,K))/4
        check('auxiliary-first-correction-'+str(alpha),
              differential/((K*K+alpha*alpha-1)*jf(K)/4)-1,mp.mpf('1e-50'))
    for KK in ['4','8','16']:
        K=mp.mpf(KK);nu=4*K**3;family={}
        for yy in (['0','1','2'] if KK=='8' else ['0','2']):
            coarse=complete(K,nu,yy,20,50)
            row=complete(K,nu,yy,36,65)
            change=abs(coarse['_total']/row['_total']-1)
            check('whole-quadrature-'+KK+'-'+yy,change,mp.mpf('1e-15'))
            row['quadrature_change']=str(change)
            if yy=='2':
                contour=complete(K,nu,yy,36,65,line_gap='.75')
                extended=complete(K,nu,yy,36,65,Xfactor=10,T=24)
                precision=complete(K,nu,yy,36,90)
                for name,other,tol in [('contour',contour,'1e-15'),('extent',extended,'1e-18'),('precision',precision,'1e-40')]:
                    err=abs(other['_total']/row['_total']-1)
                    check(name+'-'+KK,err,mp.mpf(tol));row[name+'_change']=str(err)
            family[yy]=row
            print('completed full contour K='+KK+', y='+yy,flush=True)
        for yy,row in family.items():
            if yy=='0':continue
            ref=family['0'];value=row['_total']/ref['_total']*mp.exp(row['_base']-ref['_base'])
            target=bessel_branch(K,mp.mpf(yy))/bessel_branch(K,mp.mpf(0))
            ratios.append(dict(K=KK,nu=str(nu),y=yy,exact_graph_chiral_ratio=str(value),
                               auxiliary_Bessel_ratio=str(target),relative_difference=str(abs(value/target-1)),
                               compact_kernel_phase_correction=str(mp.exp(1j*(K*mp.mpf(yy))**2/nu)),
                               phase_corrected_relative_difference=str(abs(value/(target*mp.exp(1j*(K*mp.mpf(yy))**2/nu))-1)),
                               K_branch_over_I_branch=str(abs(mp.besselk(K*mp.mpf(yy),K)/(mp.pi*mp.besseli(K*mp.mpf(yy),K))))))
        for row in family.values():
            row.pop('_total');row.pop('_base');rows.append(row)
    inputs=['code/recoil_global_probe.py','code/recoil_saddle.py','code/recoil_local.py','code/recoil_internal.py',
            'sources/dlmf-5.11.html']
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='passed' if all(c['passed'] for c in checks) else 'failed',
              scope='full finite-parameter deformed Barnes quadrature, exact reflection and auxiliary Bessel identity; no proved uniform Barnes-to-Bessel limit',
              environment=dict(python=sys.version,platform=platform.platform(),mpmath=mp.__version__,
                               dps=[50,65,90],orders=[20,36]),
              checks=checks,negative_identities=identities,critical_y=str(critical),full_rows=rows,bessel_comparisons=ratios,
              hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    out=ROOT/'results/stage4/recoil-global';out.mkdir(parents=True,exist_ok=True)
    (out/'diagnostics.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=data['status'],checks=len(checks),failed=[c['name'] for c in checks if not c['passed']])))
    if data['status']!='passed':raise SystemExit(1)

if __name__=='__main__':run()
