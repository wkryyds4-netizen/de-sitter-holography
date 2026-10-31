"""Independent source-coordinate and high-precision controls; retained errors explicit."""
import hashlib,json,platform
from pathlib import Path
import numpy as np
import mpmath as mp
from scipy.special import expi,exp1
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss
from gaussian_eikonal import j,smeared,enhancement

def direct_four(m,gamma,t,sigma,n,na):
    # Original four independent shifts; no midpoint/relative Gaussian reduction.
    z,w=hermgauss(n);d=np.sqrt(2)*sigma*z;w=w/np.sqrt(np.pi)
    shifts=np.meshgrid(d,d,d,d,indexing='ij');weights=np.meshgrid(w,w,w,w,indexing='ij')
    d1,d2,d3,d4=[v.ravel() for v in shifts];weight=np.prod(np.array(weights),axis=0).ravel()
    th1=-np.pi/2+1j*(t+d1);th2=1j*d2;th3=np.pi/2+1j*(t+d3);th4=np.pi+1j*d4
    A=np.exp(-1j*th4)-np.exp(-1j*th2);B=np.exp(1j*th1)-np.exp(1j*th3)
    pref=np.exp(.5j*(th1+th3-th2-th4));alpha=2*np.pi*(1-gamma)
    disc=(-pref/(4*alpha*alpha*A*B)).real
    aa=-(np.exp(-d4)+np.exp(-d2));bb=-1j*(np.exp(t+d1)+np.exp(t+d3))
    xp,wp=leggauss(na);phi=np.pi*(xp+1)/2;wp=wp/2
    h=-np.cos((1-gamma)*(np.pi-phi))/(2*m*np.sinc(gamma))
    def kernel(c,side):
        c=c.real;out=np.empty(c.shape,complex);small=np.abs(c)<.002
        val=np.ones(np.count_nonzero(small));term=val.copy()
        for k in range(1,16):term*=(-c[small])*k;val+=term
        out[small]=val
        pos=(c>=.002);a=1/c[pos];out[pos]=a*np.exp(a)*exp1(a)
        neg=(c<=-.002);a=-1/c[neg];out[neg]=a*np.exp(-a)*expi(a)
        neg=c<0;a=-1/c[neg];out[neg]-=side*1j*np.pi*a*np.exp(-a)
        return out
    total=0j
    for q,v in zip(h,wp):total+=v*np.sum(weight*disc*(kernel(1j*q/(A*B),1)+kernel(1j*q/(aa*bb),-1)-1))
    norm=float(np.sum(weight*disc))
    return dict(value=[float(total.real/norm),float(total.imag/norm)],baseline=norm,n=n,na=na)

def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/gaussian-eikonal';dest=base/'diagnostics.json'
    if dest.exists():raise RuntimeError('Preserve diagnostics')
    study=json.loads((base/'study.json').read_text());mp.mp.dps=70;checks=[];refs=[]
    def check(label,error,tol):checks.append(dict(label=label,error=float(error),tolerance=float(tol),passed=bool(error<=tol)))
    for c in [-8,-1,-.1,-1/128,-.005,-1e-6,0,1e-6,.005,1/128,.1,1,8]:
        x=mp.mpf(c)
        ref=1 if x==0 else (mp.exp(1/x)*mp.e1(1/x)/x if x>0 else (-1/x)*mp.exp(1/x)*(mp.ei(-1/x)-1j*mp.pi))
        got=complex(j(c));err=abs(got-complex(ref));check('J mpmath '+str(c),err,3e-14)
        refs.append(dict(c=c,value=[str(mp.re(ref)),str(mp.im(ref))],error=err))
    # Direct unintegrated G.2 rational momentum integral for positive h.
    source=[]
    for shifts in [(0,0,0,0),(.2,-.3,.6,-.1),(-.8,.5,-.2,.9)]:
        t=mp.mpf('.7');h=mp.mpf('.2');a=mp.mpf('4.5');d1,d2,d3,d4=map(mp.mpf,shifts)
        th=[-mp.pi/2+1j*(t+d1),1j*d2,mp.pi/2+1j*(t+d3),mp.pi+1j*d4]
        A=mp.exp(-1j*th[3])-mp.exp(-1j*th[1]);B=mp.exp(1j*th[0])-mp.exp(1j*th[2])
        pref=mp.exp(.5j*(th[0]+th[2]-th[1]-th[3]))/(4*a*a)
        raw=pref*mp.quad(lambda v:mp.exp(-1j*B*v)/(1j*(A-h*v)),[0,1,mp.inf])
        T=t+(d1+d3-d2-d4)/2;C=mp.cosh((d1-d3)/2)*mp.cosh((d2-d4)/2);c=h*mp.exp(T)/(4*C)
        closed=mp.exp(1/c)*mp.e1(1/c)/(16*a*a*C*c)
        check('original G.2 momentum integral',abs(raw-closed),mp.mpf('1e-60'))
        source.append(dict(shifts=list(shifts),raw=[str(mp.re(raw)),str(mp.im(raw))],closed=str(closed)))
    moments=[]
    for sigma in [.1,.2,.4,.8]:
        A=[mp.quad(lambda x:mp.exp(-x*x)/mp.sqrt(mp.pi)/mp.cosh(mp.mpf(sigma)*x)**r,[-mp.inf,0,mp.inf]) for r in [1,2,3]]
        for n in [1,2]:
            ref=mp.exp(n*n*mp.mpf(sigma)**2/2)*(A[n]/A[0])**2
            check('moment coefficient',abs(enhancement(n,sigma)-float(ref)),3e-11)
            moments.append(dict(sigma=sigma,n=n,value=str(ref)))
    direct=[]
    for sigma in [.4,.8]:
        r=direct_four(1000,.25,np.log(1000),sigma,16,48)
        f=direct_four(1000,.25,np.log(1000),sigma,24,72)
        target=next(x['fine'] for x in study['rows'] if x['gamma']==.25 and x['sigma']==sigma and x['u']==0)
        err=abs(complex(*f['value'])-complex(*target['value']))
        check('direct four Gaussian/source coordinates',err,2e-5)
        direct.append(dict(sigma=sigma,coarse=r,fine=f,reduced=target,difference=err))
        print('direct',sigma,'difference',err,flush=True)
    for row in study['rows']:
        check('coarse/fine full integral',row['difference'],2e-4)
        check('residue imaginary integral',abs(row['fine']['value'][1]-row['fine']['residue_imag']),1e-13)
        check('negative imaginary part',max(row['fine']['value'][1],0),0)
    check('zero-time imaginary',abs(study['controls'][0]['result']['value'][1]),1e-14)
    pos=next(x['fine'] for x in study['rows'] if x['gamma']==.25 and x['sigma']==.8 and x['u']==0)
    neg=study['controls'][1]['result'];check('time reversal',abs(complex(*neg['value'])-complex(*pos['value']).conjugate()),1e-13)
    out=dict(status='passed' if all(c['passed'] for c in checks) else 'failed',checks=checks,references=refs,
             original_momentum_integrals=source,moments=moments,direct_four=direct,
             scope='Independent original four-coordinate quadrature, momentum integrals, and high-precision special-function/moment checks. Empirical quadrature tolerances, not certified total error.',
             dependencies=dict(python=platform.python_version(),mpmath=mp.__version__),
             hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in ['code/gaussian_eikonal.py','code/check_gaussian_eikonal.py','results/stage4/gaussian-eikonal/study.json']})
    dest.write_text(json.dumps(out,indent=2),encoding='utf-8');print(out['status'],len(checks),flush=True)
    if out['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
