"""Transverse Barnes kernel and two-branch diagnostics; not interval certificates."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, platform, sys
import mpmath as mp
from recoil_bessel_joint import gauss_segment, log_g, log_v, parse_complex
from recoil_internal import parameters, barnes

ROOT=Path(__file__).resolve().parents[1]

def effective(K,nu,u,v):
    M=nu**2/(2*K)
    ws=[(nu-e*u-f*v)/2 for e in [-1,1] for f in [-1,1]]
    if min(ws)<=0: raise ValueError('Requires |u|+|v|<nu.')
    k=2*mp.sqrt(mp.fprod(ws))/M
    lam=mp.fsum(1/w for w in ws)/8
    return M,ws,k,lam

def log_r(t,K,nu,u,v):
    M,ws,k,lam=effective(K,nu,u,v)
    b=mp.mpf('.5')-1j*(M-nu); c=mp.mpf('.5')-1j*(M+nu)
    return (1j*mp.pi*t-2*t*mp.log(k/2)
        +mp.fsum(mp.loggamma(t-1j*w)-mp.loggamma(-1j*w) for w in ws)
        +mp.loggamma(b-t)-mp.loggamma(b)+mp.loggamma(c)-mp.loggamma(c+t))

def log_q(t,K,nu,alpha,u,v):
    M,ws,k,lam=effective(K,nu,u,v)
    b=mp.mpf('.5')-1j*(M-nu); c=mp.mpf('.5')-1j*(M+nu)
    return (-1j*mp.pi*t+mp.fsum(mp.loggamma(t-1j*w)-mp.loggamma(-1j*w) for w in ws)
        +mp.loggamma(b-t)-mp.loggamma(b)+mp.loggamma(c)-mp.loggamma(c+t)
        +mp.loggamma((1+alpha)/2-t)+mp.loggamma((1-alpha)/2-t))

def log_c(K,nu,u,v):
    M,ws,k,lam=effective(K,nu,u,v)
    return (-mp.pi*(M+nu)+mp.fsum(mp.loggamma(-1j*w) for w in ws)
        +mp.loggamma(mp.mpf('.5')-1j*(M-nu))
        -mp.loggamma(mp.mpf('.5')-1j*(M+nu)))

def recurrence(t,K,nu,alpha,u,v):
    M,ws,k,lam=effective(K,nu,u,v)
    return mp.fprod(t-1j*w for w in ws)/(((t+mp.mpf('.5')-1j*nu)**2+M*M)
        *((t+mp.mpf('.5'))**2-alpha*alpha/4))

def finite_contour(K,nu,alpha,u,v,order,dps):
    with mp.workdps(dps):
        K,nu,alpha,u,v=map(mp.mpf,(K,nu,alpha,u,v))
        M,ws,k,lam=effective(K,nu,u,v)
        T=mp.mpf('1.2'); right=4*K+mp.mpf('.37')
        L=mp.sqrt(k*k+alpha*alpha)/2; N=int(mp.ceil(L)); d=mp.mpf('.5')-L
        A=-N-1j*T; B=-N+1j; C=right+1j
        D=d-1j*T; E=d+1j*T; F=right+1j*T; H=right-1j*T
        nodes,weights=mp.gauss_quadrature(order,'legendre')
        def integ(fn,a,b):return gauss_segment(fn,a,b,nodes,weights,mp.mpf('.2'))[0]
        fq=lambda t:mp.exp(log_q(t,K,nu,alpha,u,v))
        fg=lambda t:mp.exp(log_r(t,K,nu,u,v)+log_g(t,k,alpha))
        fv=lambda t:mp.exp(log_r(t,K,nu,u,v)+log_v(t,k,alpha))
        original=(integ(fq,A,B)+integ(fq,B,C))/(2j*mp.pi)
        vertical=integ(fg,D,E)/(2j*mp.pi)
        kb=(integ(fg,A,D)+integ(fg,E,F)+integ(fg,F,C))/(2j*mp.pi)
        vb=(integ(fv,C,H)+integ(fv,H,A))/(2j*mp.pi)
        a=(1+alpha)/2
        terms=[mp.exp(log_r(a+n,K,nu,u,v)+2*(a+n)*mp.log(k/2)
            -mp.loggamma(n+1)-mp.loggamma(alpha+n+1)) for n in range(int(mp.floor(right-a))+1)]
        main=-mp.exp(1j*mp.pi*alpha)*vertical-2j*mp.pi*mp.fsum(terms)
        correction=-mp.exp(1j*mp.pi*alpha)*kb-2*mp.pi*mp.exp(1j*mp.pi*alpha/2)*vb
        envelope=k*(mp.besselk(alpha,k)+mp.pi*mp.besseli(alpha,k))
        return dict(K=str(K),nu=str(nu),alpha=str(alpha),u=str(u),v=str(v),T=str(T),right=str(right),
            order=order,dps=dps,original=str(original),main=str(main),correction=str(correction),
            envelope=str(envelope),identity_error=str(abs(original-main-correction)/envelope),
            boundary_over_envelope=str(abs(correction)/envelope))

def branches(K,nu,alpha,u,v,order=28,dps=70,height=8):
    with mp.workdps(dps):
        K,nu,alpha,u,v=map(mp.mpf,(K,nu,alpha,u,v))
        M,ws,k,lam=effective(K,nu,u,v)
        L=mp.sqrt(k*k+alpha*alpha)/2; d=mp.mpf('.5')-L
        T=min(height*K,min(ws)/2,(M-nu)/2)
        base=max(36*K,8*K+2*alpha)
        right=max([base+mp.mpf(j)/8 for j in range(8)],key=lambda x:
            min(abs(mp.cos(mp.pi*(x-alpha/2))),abs(mp.cos(mp.pi*(x+alpha/2)))))
        ek=k*mp.besselk(alpha,k); ei=k*mp.besseli(alpha,k); envelope=ek+mp.pi*ei
        nodes,weights=mp.gauss_quadrature(order,'legendre')
        f=lambda s:mp.exp(log_r(d+1j*s,K,nu,u,v)+log_g(d+1j*s,k,alpha)-mp.log(ek))/(2*mp.pi)
        vk,avk=gauss_segment(f,-T,T,nodes,weights,step=max(1,mp.sqrt(k)/2))
        a=(1+alpha)/2; terms=[]; prob=[]
        for n in range(int(mp.floor(right-a))+1):
            t=a+n
            wn=2*mp.exp(2*t*mp.log(k/2)-mp.loggamma(n+1)-mp.loggamma(alpha+n+1))/ei
            prob.append(wn);terms.append(wn*mp.exp(log_r(t,K,nu,u,v)))
        vi=mp.fsum(terms); avi=mp.fsum(abs(x) for x in terms)
        phase=mp.exp(1j*lam*(k*k+alpha*alpha-1))
        val=-mp.exp(1j*mp.pi*alpha)*ek*vk-1j*mp.pi*ei*vi
        lead=phase*(-mp.exp(1j*mp.pi*alpha)*ek-1j*mp.pi*ei)
        return val,dict(K=str(K),nu=str(nu),alpha=str(alpha),u=str(u),v=str(v),k=str(k),lam=str(lam),
            theta=str((abs(u)+abs(v))/nu),rho=str(K**3/nu**2),order=order,dps=dps,T=str(T),right=str(right),
            value=str(val),leading=str(lead),envelope=str(envelope),
            K_factor=str(vk),I_factor=str(vi),K_condition=str(avk/abs(vk)),I_condition=str(avi/abs(vi)),
            sum_condition=str((ek*abs(vk)+mp.pi*ei*abs(vi))/abs(val)),
            error_over_rho=str(abs(val-lead)/(envelope*K**3/nu**2)),
            envelope_error=str(abs(val-lead)/envelope),auxiliary_I_tail=str(1-mp.fsum(prob)))

def run():
    mp.mp.dps=85; checks=[]; algebra=[]; finite=[]; joint=[]; transverse=[]; tails=[]
    def check(name,error,tol):
        checks.append(dict(name=name,error=str(abs(error)),tolerance=str(tol),passed=bool(abs(error)<=tol)))
    for ui,vi in [('0','0'),('.3','.2'),('-.25','.4')]:
        K=mp.mpf(4);nu=mp.mpf(256);alpha=mp.mpf('3.7');u=nu*mp.mpf(ui);v=nu*mp.mpf(vi)
        M,ws,k,lam=effective(K,nu,u,v)
        check('effective-product-'+ui+vi,k/(mp.sqrt((nu**2-(u+v)**2)*(nu**2-(u-v)**2))/(2*M))-1,mp.mpf('1e-70'))
        cyc=[M+1j*alpha/2+u,M-1j*alpha/2+v,M+1j*alpha/2-u,M-1j*alpha/2-v]
        aa,rr=parameters([cyc[0],nu,cyc[1],cyc[2],nu,cyc[3]])
        for t in [mp.mpc('-1.7','.9'),mp.mpc('3.2','-2.1')]:
            q=mp.exp(log_q(t,K,nu,alpha,u,v))
            original=barnes(t-1j*(M+nu),aa,rr)/mp.exp(log_c(K,nu,u,v))
            check('raw-normalization-'+ui+vi+str(t),q/original-1,mp.mpf('1e-65'))
            check('exact-recurrence-'+ui+vi+str(t),mp.exp(log_q(t+1,K,nu,alpha,u,v)-log_q(t,K,nu,alpha,u,v))/recurrence(t,K,nu,alpha,u,v)-1,mp.mpf('1e-65'))
            split=mp.exp(log_r(t,K,nu,u,v))*(-mp.exp(1j*mp.pi*alpha+log_g(t,k,alpha))+2*mp.pi*mp.exp(1j*mp.pi*alpha/2+log_v(t,k,alpha)))
            check('exact-split-'+ui+vi+str(t),q/split-1,mp.mpf('1e-65'))
        exact=abs(mp.exp(log_c(K,nu,u,v)-log_c(K,nu,0,0)))
        predicted=K/k*(1-mp.exp(-mp.pi*nu))**2/mp.sqrt(mp.fprod(1-mp.exp(-2*mp.pi*w) for w in ws))
        check('normalization-modulus-'+ui+vi,exact/predicted-1,mp.mpf('1e-65'))
        algebra.append(dict(u=str(u),v=str(v),K=str(K),nu=str(nu),k=str(k),lam=str(lam),C_modulus_ratio=str(exact)))
    for yi in [0,2]:
        K=mp.mpf(8);nu=32*K**mp.mpf('1.5');alpha=yi*K;u=mp.mpf('.4')*nu;v=mp.mpf('.3')*nu
        M,ws,k,lam=effective(K,nu,u,v);dmax=max(alpha*alpha/4-1,0);A=max(2,8*K,2*dmax);R=(32+2*yi)*K
        N=int(mp.ceil(mp.sqrt(k*k+alpha*alpha)/2));rate=2*mp.pi-N*(1/(M-nu)+1/(M+nu))
        for x in [R,A,2*A,10*A]:
            actual=abs(recurrence(x+1j,K,nu,alpha,u,v));bound=(x/(x+1))**mp.mpf('1.5')
            check(f'positive-tail-{yi}-{x}',max(mp.mpf(0),actual/bound-1),mp.mpf('1e-65'))
            tails.append(dict(branch='positive',K=str(K),nu=str(nu),alpha=str(alpha),u=str(u),v=str(v),x=str(x),ratio=str(actual),bound=str(bound)))
        for s in [1,K,M-nu,2*M]:
            logratio=mp.re(log_q(-N-1j*s,K,nu,alpha,u,v)-log_q(-N,K,nu,alpha,u,v))+rate*s
            check(f'negative-tail-{yi}-{s}',max(mp.mpf(0),logratio),mp.mpf('1e-60'))
            tails.append(dict(branch='negative',s=str(s),K=str(K),nu=str(nu),alpha=str(alpha),u=str(u),v=str(v),N=N,log_ratio_to_bound=str(logratio)))
    for alpha,u,v in [('.7','12','7'),('4','25','-16')]:
        low=finite_contour(2,128,alpha,u,v,20,60);hi=finite_contour(2,128,alpha,u,v,32,80)
        check('finite-identity-'+alpha,mp.mpf(hi['identity_error']),mp.mpf('1e-32'))
        check('finite-refine-'+alpha,abs(parse_complex(low['original'])-parse_complex(hi['original']))/mp.mpf(hi['envelope']),mp.mpf('1e-22'))
        finite.append(hi)
    for kind in ['natural','fractional']:
        for yi in [0,2]:
            for ki in [8,16,32]:
                K=mp.mpf(ki);nu=8*K**mp.mpf('1.75');M=nu*nu/(2*K);alpha=yi*K
                u,v=(mp.sqrt(M),mp.sqrt(M)/2) if kind=='natural' else (nu*mp.mpf('.3'),nu*mp.mpf('.2'))
                low,_=branches(K,nu,alpha,u,v,20,55)
                val,row=branches(K,nu,alpha,u,v,32,75)
                change=abs(low-val)/mp.mpf(row['envelope']);row['quadrature_change']=str(change);row['family']=kind
                check(f'joint-refine-{kind}-{yi}-{ki}',change,mp.mpf('1e-17'))
                if ki==32:
                    hi,_=branches(K,nu,alpha,u,v,32,95);ext,_=branches(K,nu,alpha,u,v,32,75,height=10)
                    pe=abs(hi-val)/mp.mpf(row['envelope']);ee=abs(ext-val)/mp.mpf(row['envelope'])
                    check(f'precision-{kind}-{yi}',pe,mp.mpf('1e-55'));check(f'extent-{kind}-{yi}',ee,mp.mpf('1e-25'))
                    row['precision_change']=str(pe);row['extent_change']=str(ee)
                joint.append(row)
                print(f'{kind} K={ki} y={yi} error/rho={row["error_over_rho"][:16]}',flush=True)
    # Independent Bessel functions, fixed alpha during changes in their argument.
    for yi in [0,1,2]:
        errors=[]
        for ki in [16,64,256]:
            K=mp.mpf(ki);nu=8*K**mp.mpf('1.75');M=nu*nu/(2*K);alpha=yi*K
            u=mp.sqrt(M);v=mp.sqrt(M)/2;_,_,k,lam=effective(K,nu,u,v)
            h=mp.sqrt(1+yi*yi)*mp.mpf('1.25')/2
            ri=k*mp.besseli(alpha,k)/(K*mp.besseli(alpha,K));rk=k*mp.besselk(alpha,k)/(K*mp.besselk(alpha,K))
            error=max(abs(mp.log(ri)+h),abs(mp.log(rk)-h));errors.append(error)
            transverse.append(dict(K=str(K),y=yi,k=str(k),I_ratio=str(ri),K_ratio=str(rk),I_limit=str(mp.exp(-h)),K_limit=str(mp.exp(h)),log_error=str(error),K_log_error=str(K*error)))
        check(f'transverse-rate-{yi}',errors[-1]/errors[0],mp.mpf('.08'))
    inputs=['code/recoil_transverse.py','code/recoil_bessel_joint.py','code/recoil_bessel_limit.py','code/recoil_internal.py',
        'sources/chandra-liska-post-2610.30247v1.html','sources/dlmf-5.11.html']
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='passed' if all(c['passed'] for c in checks) else 'failed',
        scope='transverse internal chiral Barnes graph; empirical diagnostics, no interval or physical outer-integral certificate',
        environment=dict(python=sys.version,platform=platform.platform(),mpmath=mp.__version__,dps=[55,60,75,80,85,95],orders=[20,32]),
        checks=checks,algebra_rows=algebra,tail_rows=tails,finite_contours=finite,joint_rows=joint,transverse_rows=transverse,
        hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    out=ROOT/'results/stage4/recoil-transverse';out.mkdir(parents=True,exist_ok=True)
    (out/'diagnostics.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=data['status'],checks=len(checks),failed=[c for c in checks if not c['passed']])))
    if data['status']!='passed':raise SystemExit(1)

if __name__=='__main__':run()
