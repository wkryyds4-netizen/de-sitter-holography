"""Fixed-K full Barnes limit: recurrence/tail ingredients and direct quadrature.

Numerical checks accompany the analytic estimates; they are not interval proofs.
"""
from pathlib import Path
from datetime import datetime, timezone
import json, hashlib, sys, platform, re
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1]
NODES={}

def log_q(t,K,nu,alpha):
    """Log of B(t-i(M+nu))/C, using the exact unreduced Gamma kernel."""
    M=nu*nu/(2*K);v=nu/2
    b=mp.mpf('.5')-1j*(M-nu);c=mp.mpf('.5')-1j*(M+nu)
    return (-1j*mp.pi*t+4*(mp.loggamma(t-1j*v)-mp.loggamma(-1j*v))
            +mp.loggamma(b-t)-mp.loggamma(b)+mp.loggamma(c)-mp.loggamma(c+t)
            +mp.loggamma(mp.mpf('.5')-t+alpha/2)
            +mp.loggamma(mp.mpf('.5')-t-alpha/2))

def auxiliary(t,K,alpha):
    return mp.exp(-2j*mp.pi*t+2*t*mp.log(K/2)
                  +mp.loggamma(mp.mpf('.5')-t+alpha/2)
                  +mp.loggamma(mp.mpf('.5')-t-alpha/2))

def recurrence(t,K,nu,alpha):
    M=nu*nu/(2*K)
    return (t-1j*nu/2)**4/(((t+mp.mpf('.5')-1j*nu)**2+M*M)
                            *((t+mp.mpf('.5'))**2-alpha*alpha/4))

def J(K,alpha):
    return -K*(mp.exp(1j*mp.pi*alpha)*mp.besselk(alpha,K)
               +1j*mp.pi*mp.besseli(alpha,K))

def complete(K,nu,alpha,order=24,dps=60,X=36,T=18,eta=1):
    with mp.workdps(dps):
        K,nu,alpha,eta=map(mp.mpf,(K,nu,alpha,eta))
        N=int(mp.ceil(abs(alpha)/2))+1
        key=(order,dps)
        if key not in NODES: NODES[key]=mp.gauss_quadrature(order,'legendre')
        nodes,weights=NODES[key]
        def integ(fn,cuts):
            total=0
            for lo,hi in zip(cuts[:-1],cuts[1:]):
                for nd,wt in zip(nodes,weights):
                    total+=wt*(hi-lo)/2*fn((hi+lo)/2+(hi-lo)*nd/2)
            return total
        hc=list(map(mp.mpf,range(-N,X+1)))
        vc=sorted(set(map(mp.mpf,[-eta,0,1,2,4,8,12,T])))
        horizontal=integ(lambda x:mp.exp(log_q(x+1j*eta,K,nu,alpha)),hc)
        vertical=integ(lambda s:mp.exp(log_q(-N-1j*s,K,nu,alpha)),vc)
        total=(horizontal+1j*vertical)/(2j*mp.pi)
        return total,dict(K=str(K),nu=str(nu),alpha=str(alpha),N=N,eta=str(eta),
                          order=order,dps=dps,X=X,T=T,value=str(total),
                          horizontal=str(horizontal),vertical=str(vertical),
                          cancellation=str((abs(horizontal)+abs(vertical))/abs(horizontal+1j*vertical)))

def run():
    mp.mp.dps=65
    checks=[];bounds=[];points=[];rows=[]
    def check(name,error,tol):
        checks.append(dict(name=name,error=str(abs(error)),tolerance=str(tol),
                           passed=bool(abs(error)<=tol)))
    for K,alpha in [(mp.mpf(2),mp.mpf('.7')),(mp.mpf(4),mp.mpf(8)),
                    (mp.mpf('1.5'),mp.mpf('-3.2'))]:
        eta=mp.mpf(1);d=max(mp.mpf(0),alpha*alpha/4-eta*eta)
        A=max(mp.mpf(2),2*K,2*d)
        for nu in [mp.mpf(64),mp.mpf(256)]:
            for t in [mp.mpc('-2.3','-1.1'),mp.mpc('3.7','1')]:
                err=mp.exp(log_q(t+1,K,nu,alpha)-log_q(t,K,nu,alpha))/recurrence(t,K,nu,alpha)-1
                check(f'recurrence-{K}-{alpha}-{nu}-{t}',err,mp.mpf('1e-50'))
            M=nu*nu/(2*K);a=nu/2-eta;Nmass=M-nu+eta
            rhs=(1+1/(2*A))**2*(2*a*a+a**4/A**2)
            check(f'mass-inequality-{K}-{alpha}-{nu}',max(mp.mpf(0),rhs/Nmass**2-1),0)
            for x in [A,A+mp.mpf('.3'),2*A,max(A,nu),max(A,M),max(A,10000*M)]:
                t=x+1j*eta;ratio=abs(recurrence(t,K,nu,alpha))
                major=(x/(x+1))**mp.mpf('1.5')
                check(f'tail-step-{K}-{alpha}-{nu}-{x}',max(mp.mpf(0),ratio/major-1),0)
                bounds.append(dict(K=str(K),alpha=str(alpha),nu=str(nu),x=str(x),
                                   exact_modulus=str(ratio),majorant=str(major)))
    for K,alpha in [(mp.mpf(2),mp.mpf('.7')),(mp.mpf(4),mp.mpf(8))]:
        for nu in map(mp.mpf,[64,128,256]):
            # Fixed points only; the note proves the logarithmically growing box.
            t=mp.mpc('2.3','1')
            value=mp.exp(log_q(t,K,nu,alpha))/auxiliary(t,K,alpha)
            first=1+4j*(t*t-t)/nu
            points.append(dict(K=str(K),alpha=str(alpha),nu=str(nu),
                               compact_residual=str(abs(value-first)),
                               scaled_nu2_residual=str(nu*nu*abs(value-first))))
    for K,alpha,ns in [(2,'.7',[64,128,256]),(4,8,[128,256,512])]:
        family=[]
        for nu in ns:
            coarse,_=complete(K,nu,alpha,18,50)
            val,row=complete(K,nu,alpha,30,65)
            err=abs(coarse/val-1)
            check(f'quadrature-{K}-{alpha}-{nu}',err,mp.mpf('1e-14'))
            row['quadrature_change']=str(err)
            if nu==ns[-1]:
                high,_=complete(K,nu,alpha,30,85)
                extended,_=complete(K,nu,alpha,30,65,X=44,T=24,eta=mp.mpf('.75'))
                check(f'precision-{K}-{alpha}',abs(high/val-1),mp.mpf('1e-40'))
                check(f'contour-extent-{K}-{alpha}',abs(extended/val-1),mp.mpf('1e-16'))
                row['precision_change']=str(abs(high/val-1))
                row['contour_extent_change']=str(abs(extended/val-1))
            KK,aa=mp.mpf(K),mp.mpf(alpha);target=J(KK,aa)
            first=target*(1+1j*(KK*KK+aa*aa-1)/nu)
            row.update(auxiliary=str(target),relative_leading_error=str(abs(val/target-1)),
                       first_absolute_error=str(abs(val-first)),
                       first_relative_error=str(abs(val/first-1)),
                       scaled_nu2_absolute_error=str(nu*nu*abs(val-first)))
            family.append((nu,abs(val-first)));rows.append(row)
            print(f'completed K={K}, alpha={alpha}, nu={nu}',flush=True)
        for (n0,e0),(n1,e1) in zip(family[:-1],family[1:]):
            # Finite empirical convergence, not the analytic O(nu^-2) proof.
            check(f'fixed-K-first-error-decreases-{K}-{alpha}-{n1}',max(mp.mpf(0),e1/e0-mp.mpf('.4')),0)
    # Match the original z-contour computation, including C and 2*pi*i.
    previous=json.loads((ROOT/'results/stage4/recoil-global/diagnostics.json').read_text(encoding='utf-8'))
    old=next(r for r in previous['full_rows'] if r['K']=='4.0' and r['y']=='2.0')
    def parse_complex(s):
        vals=re.findall(r'[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?',s.replace(' ',''))
        if len(vals)!=2:raise ValueError(s)
        return mp.mpc(mp.mpf(vals[0]),mp.mpf(vals[1]))
    K,nu=mp.mpf(4),mp.mpf(256);M=nu*nu/(2*K)
    logC=(-mp.pi*(M+nu)+4*mp.loggamma(-1j*nu/2)
          +mp.loggamma(mp.mpf('.5')-1j*(M-nu))
          -mp.loggamma(mp.mpf('.5')-1j*(M+nu)))
    legacy=mp.exp(parse_complex(old['base'])-logC)*parse_complex(old['scaled_total'])/(2j*mp.pi)
    now=next(r for r in rows if r['K']=='4.0' and r['nu']=='256.0')
    legacy_error=abs(legacy/parse_complex(now['value'])-1)
    check('original-z-contour-normalization',legacy_error,mp.mpf('1e-18'))
    inputs=['code/recoil_bessel_limit.py','code/recoil_global_probe.py',
            'results/stage4/recoil-global/diagnostics.json','sources/dlmf-5.11.html']
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
              status='passed' if all(c['passed'] for c in checks) else 'failed',
              scope='finite-point diagnostics supporting a fixed-K uniform-tail proof; no interval certificate or growing-K assertion',
              environment=dict(python=sys.version,platform=platform.platform(),mpmath=mp.__version__,
                               orders=[18,30],dps=[50,65,85]),
              checks=checks,tail_bound_rows=bounds,compact_rows=points,integral_rows=rows,
              original_contour_comparison=dict(K='4',nu='256',alpha='8',relative_change=str(legacy_error)),
              hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    out=ROOT/'results/stage4/recoil-bessel-limit';out.mkdir(parents=True,exist_ok=True)
    (out/'diagnostics.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=data['status'],checks=len(checks),
                         failed=[c['name'] for c in checks if not c['passed']])))
    if data['status']!='passed':raise SystemExit(1)

if __name__=='__main__':run()
