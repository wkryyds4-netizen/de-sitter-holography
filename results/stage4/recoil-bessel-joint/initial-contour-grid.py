"""Phase-preserving Barnes split; joint-window diagnostics, not interval certificates."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, platform, re, sys
import mpmath as mp
from recoil_bessel_limit import log_q, J

ROOT=Path(__file__).resolve().parents[1]

def log_r(t,K,nu):
    M=nu*nu/(2*K); v=nu/2
    b=mp.mpf('.5')-1j*(M-nu); c=mp.mpf('.5')-1j*(M+nu)
    return (1j*mp.pi*t-2*t*mp.log(K/2)
            +4*(mp.loggamma(t-1j*v)-mp.loggamma(-1j*v))
            +mp.loggamma(b-t)-mp.loggamma(b)+mp.loggamma(c)-mp.loggamma(c+t))

def log_g(t,K,alpha):
    return (2*t*mp.log(K/2)+mp.loggamma(mp.mpf('.5')-t+alpha/2)
            +mp.loggamma(mp.mpf('.5')-t-alpha/2))

def log_v(t,K,alpha):
    return (-1j*mp.pi*t+2*t*mp.log(K/2)
            +mp.loggamma(mp.mpf('.5')-t+alpha/2)
            -mp.loggamma(mp.mpf('.5')+t+alpha/2))

def parse_complex(s):
    v=re.findall(r'[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?',s.replace(' ',''))
    if len(v)!=2: raise ValueError(s)
    return mp.mpc(*map(mp.mpf,v))

def gauss_segment(fn,a,b,nodes,weights,step=1):
    pieces=max(1,int(mp.ceil(abs(b-a)/step))); ans=0; absolute=0
    for j in range(pieces):
        lo=a+(b-a)*j/pieces; hi=a+(b-a)*(j+1)/pieces
        for z,w in zip(nodes,weights):
            value=fn((lo+hi)/2+(hi-lo)*z/2)
            ans+=(hi-lo)*w*value/2
            absolute+=abs(hi-lo)*w*abs(value)/2
    return ans,absolute

def safe_right(K,alpha):
    # At least one of these eight offsets avoids both real pole families.
    best=max([12*K+mp.mpf(j)/8 for j in range(8)],
             key=lambda x:min(abs(mp.cos(mp.pi*(x-alpha/2))),
                              abs(mp.cos(mp.pi*(x+alpha/2)))))
    return best

def branch_value(K,nu,alpha,order=24,dps=65,height_factor=8):
    with mp.workdps(dps):
        K,nu,alpha=map(mp.mpf,(K,nu,alpha)); y=alpha/K
        L=mp.sqrt(K*K+alpha*alpha)/2; d=mp.mpf('.5')-L
        T=min(height_factor*K,nu/4); right=safe_right(K,alpha)
        nodes,weights=mp.gauss_quadrature(order,'legendre')
        ek=K*mp.besselk(alpha,K); ei=K*mp.besseli(alpha,K); E=ek+mp.pi*ei
        logek=mp.log(ek)
        # Scale each analytic branch by its own positive Bessel magnitude.
        f=lambda s:mp.exp(log_r(d+1j*s,K,nu)+log_g(d+1j*s,K,alpha)-logek)/(2*mp.pi)
        vk,abs_vk=gauss_segment(f,-T,T,nodes,weights,step=max(1,mp.sqrt(K)/2))
        a=(1+alpha)/2; nmax=int(mp.floor(right-a))
        terms=[]; weights_only=[]
        for n in range(nmax+1):
            t=a+n
            lw=2*t*mp.log(K/2)-mp.loggamma(n+1)-mp.loggamma(alpha+n+1)
            wn=2*mp.exp(lw)/ei
            weights_only.append(wn); terms.append(wn*mp.exp(log_r(t,K,nu)))
        vi=mp.fsum(terms); abs_vi=mp.fsum(abs(x) for x in terms)
        value=-mp.exp(1j*mp.pi*alpha)*ek*vk-1j*mp.pi*ei*vi
        phase=mp.exp(1j*(K*K+alpha*alpha-1)/nu); lead=phase*J(K,alpha)
        err=abs(value-lead)/E; rho=K**3/nu**2
        return value,dict(K=str(K),nu=str(nu),alpha=str(alpha),y=str(y),rho=str(rho),
                         d=str(d),T=str(T),right=str(right),nmax=nmax,order=order,dps=dps,
                         value=str(value),leading=str(lead),envelope=str(E),
                         K_branch_factor=str(vk),I_branch_factor=str(vi),
                         K_branch_condition=str(abs_vk/abs(vk)),I_branch_condition=str(abs_vi/abs(vi)),
                         sum_condition=str((ek*abs(vk)+mp.pi*ei*abs(vi))/abs(value)),
                         envelope_error=str(err),error_over_rho=str(err/rho),
                         I_weight_tail=str(1-mp.fsum(weights_only)),
                         scope='two finite analytic branches; connector error controlled analytically, sampled by extent changes')

def finite_contour(K,nu,alpha,T,right,order,dps):
    with mp.workdps(dps):
        K,nu,alpha,T,right=map(mp.mpf,(K,nu,alpha,T,right))
        L=mp.sqrt(K*K+alpha*alpha)/2; N=int(mp.ceil(L)); d=mp.mpf('.5')-L
        A=-N-1j*T; Bb=-N+1j; C=right+1j
        D=d-1j*T; E=d+1j*T; F=right+1j*T; G=right-1j*T
        nodes,weights=mp.gauss_quadrature(order,'legendre')
        def integ(fn,a,b):return gauss_segment(fn,a,b,nodes,weights,mp.mpf('.5'))[0]
        fq=lambda t:mp.exp(log_q(t,K,nu,alpha))
        fg=lambda t:mp.exp(log_r(t,K,nu)+log_g(t,K,alpha))
        fv=lambda t:mp.exp(log_r(t,K,nu)+log_v(t,K,alpha))
        original=(integ(fq,A,Bb)+integ(fq,Bb,C))/(2j*mp.pi)
        vertical=integ(fg,D,E)/(2j*mp.pi)
        kb=(integ(fg,A,D)+integ(fg,E,F)+integ(fg,F,C))/(2j*mp.pi)
        vb=(integ(fv,C,G)+integ(fv,G,A))/(2j*mp.pi)
        terms=[]
        a=(1+alpha)/2
        for n in range(int(mp.floor(right-a))+1):
            t=a+n
            terms.append(mp.exp(log_r(t,K,nu)+2*t*mp.log(K/2)
                                -mp.loggamma(n+1)-mp.loggamma(alpha+n+1)))
        main=-mp.exp(1j*mp.pi*alpha)*vertical-2j*mp.pi*mp.fsum(terms)
        correction=-mp.exp(1j*mp.pi*alpha)*kb-2*mp.pi*mp.exp(1j*mp.pi*alpha/2)*vb
        scale=K*(mp.besselk(alpha,K)+mp.pi*mp.besseli(alpha,K))
        return dict(K=str(K),nu=str(nu),alpha=str(alpha),T=str(T),right=str(right),
                    order=order,dps=dps,original=str(original),branch_main=str(main),
                    boundary_correction=str(correction),envelope=str(scale),
                    normalized_identity_error=str(abs(original-main-correction)/scale),
                    boundary_size_over_E=str(abs(correction)/scale))

def run():
    mp.mp.dps=75; checks=[]; finite=[]; moments=[]; joint=[]; legacy=[]
    def check(name,err,tol):
        checks.append(dict(name=name,error=str(abs(err)),tolerance=str(tol),passed=bool(abs(err)<=tol)))
    for K,alpha in [(2,mp.mpf('.7')),(4,8),(4,mp.mpf('3.5'))]:
        K=mp.mpf(K); alpha=mp.mpf(alpha); nu=32*K**mp.mpf('1.5')
        for t in [mp.mpc('-1.7','.9'),mp.mpc('3.2','-2.1')]:
            g=mp.exp(log_g(t,K,alpha)); v=mp.exp(log_v(t,K,alpha))
            lhs=mp.exp(-2j*mp.pi*t)*g
            rhs=-mp.exp(1j*mp.pi*alpha)*g+2*mp.pi*mp.exp(1j*mp.pi*alpha/2)*v
            check(f'split-{K}-{alpha}-{t}',abs(lhs-rhs)/(abs(lhs)+abs(rhs)),mp.mpf('1e-65'))
            check(f'exact-multiplier-{K}-{alpha}-{t}',
                  mp.exp(log_r(t,K,nu)+log_g(t,K,alpha)-2j*mp.pi*t-log_q(t,K,nu,alpha))-1,
                  mp.mpf('1e-65'))
        # Positive series moment identities, summed beyond the numerical right edge.
        terms=[]; m1=[]; m2=[]; tmean=[]; nm=[]
        S=K*K+alpha*alpha-1
        for n in range(int(30*K+50)):
            t=(1+alpha)/2+n; w=2*mp.exp(2*t*mp.log(K/2)-mp.loggamma(n+1)-mp.loggamma(alpha+n+1))
            q=4*(t*t-t)-S; terms.append(w);m1.append(w*q);m2.append(w*q*q);tmean.append(w*t);nm.append(w*n*(n+alpha))
        norm=mp.fsum(terms)
        check(f'positive-series-{K}-{alpha}',norm/(K*mp.besseli(alpha,K))-1,mp.mpf('1e-65'))
        check(f'centered-first-moment-{K}-{alpha}',mp.fsum(m1)/norm,mp.mpf('1e-65'))
        check(f'centered-second-moment-{K}-{alpha}',mp.fsum(m2)/(8*K*K*mp.fsum(tmean))-1,mp.mpf('1e-65'))
        check(f'n-moment-{K}-{alpha}',mp.fsum(nm)/norm-K*K/4,mp.mpf('1e-65'))
        moments.append(dict(K=str(K),alpha=str(alpha),first=str(mp.fsum(m1)/norm),
                            second=str(mp.fsum(m2)/norm),second_upper=str(4*K*K*(1+alpha+K))))
    for K,alpha in [(2,'.7'),(2,4),(3,'3.5')]:
        lo=finite_contour(K,128,alpha,'1.2',mp.mpf(4*K)+mp.mpf('.37'),20,60)
        hi=finite_contour(K,128,alpha,'1.2',mp.mpf(4*K)+mp.mpf('.37'),32,80)
        check(f'finite-contour-identity-{K}-{alpha}',mp.mpf(hi['normalized_identity_error']),mp.mpf('1e-35'))
        check(f'finite-contour-refine-{K}-{alpha}',
              abs(parse_complex(lo['original'])-parse_complex(hi['original']))/mp.mpf(hi['envelope']),mp.mpf('1e-25'))
        finite.append(hi)
    previous=json.loads((ROOT/'results/stage4/recoil-bessel-limit/diagnostics.json').read_text())
    for old in previous['integral_rows']:
        K,alpha,nu=map(mp.mpf,(old['K'],old['alpha'],old['nu']))
        value,row=branch_value(K,nu,alpha,32,75)
        err=abs(value-parse_complex(old['value']))/mp.mpf(row['envelope'])
        check(f'legacy-full-{K}-{alpha}-{nu}',err,mp.mpf('1e-13'))
        row['legacy_envelope_change']=str(err);legacy.append(row)
    # A genuine rho -> 0 path, with K^2/nu growing: nu=8*K^(7/4).
    for yi in [0,1,2]:
        for ki in [8,16,32,64]:
            K=mp.mpf(ki);alpha=K*yi;nu=8*K**mp.mpf('1.75')
            low,_=branch_value(K,nu,alpha,20,55)
            value,row=branch_value(K,nu,alpha,32,75)
            change=abs(low-value)/mp.mpf(row['envelope'])
            check(f'joint-quadrature-{ki}-{yi}',change,mp.mpf('1e-18'))
            row['quadrature_envelope_change']=str(change)
            if ki==64:
                high,_=branch_value(K,nu,alpha,32,95)
                ext,_=branch_value(K,nu,alpha,32,75,height_factor=10)
                pe=abs(high-value)/mp.mpf(row['envelope']);ee=abs(ext-value)/mp.mpf(row['envelope'])
                check(f'joint-precision-{yi}',pe,mp.mpf('1e-55'))
                check(f'joint-extent-{yi}',ee,mp.mpf('1e-35'))
                row['precision_envelope_change']=str(pe);row['extent_envelope_change']=str(ee)
            joint.append(row)
            print(f'joint K={ki}, y={yi}, error/rho={row["error_over_rho"][:16]}',flush=True)
    # Approach an auxiliary cancellation zero without dividing by J.
    from recoil_bessel_cancellation import root
    alpha=mp.mpf('15.5'); K=root(alpha); nu=8*K**mp.mpf('1.75')
    zero,zero_row=branch_value(K,nu,alpha,32,80)
    zhi,_=branch_value(K,nu,alpha,40,100)
    check('zero-envelope-refine',abs(zero-zhi)/mp.mpf(zero_row['envelope']),mp.mpf('1e-25'))
    inputs=['code/recoil_bessel_joint.py','code/recoil_bessel_limit.py','code/recoil_bessel_cancellation.py',
            'results/stage4/recoil-bessel-limit/diagnostics.json',
            'notes/stage4-recoil-bessel-endpoints.md','sources/dlmf-5.11.html']
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
              status='passed' if all(c['passed'] for c in checks) else 'failed',
              scope='phase-preserving branch split and joint rho-to-zero diagnostics; no interval or finite-G certificate',
              environment=dict(python=sys.version,platform=platform.platform(),mpmath=mp.__version__,
                               orders=[20,32,40],dps=[55,60,75,80,95,100]),
              joint_path='nu=8*K^(7/4); rho=1/(64*sqrt(K)); K^2/nu=K^(1/4)/8',
              checks=checks,finite_contours=finite,positive_moments=moments,legacy_rows=legacy,
              joint_rows=joint,zero_row=zero_row,
              hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    out=ROOT/'results/stage4/recoil-bessel-joint';out.mkdir(parents=True,exist_ok=True)
    (out/'diagnostics.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=data['status'],checks=len(checks),failed=[c['name'] for c in checks if not c['passed']])))
    if data['status']!='passed':raise SystemExit(1)

if __name__=='__main__':run()
