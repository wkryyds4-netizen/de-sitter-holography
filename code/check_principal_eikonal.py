"""Independent formulas, deterministic integrals and recorded numerical controls."""
import hashlib,json,math
from pathlib import Path
import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.hermite import hermgauss
from principal_eikonal import kernel,spectral_moment
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/stage4/principal-eikonal'

def pair_rapidity(n,nu,sigma,order=80,hn=48,R=14):
    z,w=leggauss(order);r=R*(z+1)/2;w=w*R/2
    x,hw=hermgauss(hn);x*=sigma;hw/=np.sqrt(np.pi)
    rr=np.cosh(r);ww=w[:,None]*w[None,:]*np.cos(nu*r)[:,None]*np.cos(nu*r)[None,:]
    return float(sum(v*math.factorial(n+1)*np.sum(ww/(np.exp(y)*rr[:,None]+np.exp(-y)*rr[None,:])**(n+2)) for y,v in zip(x,hw)))

def fixed_rapidity(lam,nu,nup,X,Y,order,R=14,conformal=False):
    z,w=leggauss(order);r=R*(z+1)/2;w*=R/2
    cr=np.cosh(r)
    A=(np.exp(Y)*cr[:,None]+np.exp(-Y)*cr[None,:]).ravel()
    B=(np.exp(X)*cr[:,None]+np.exp(-X)*cr[None,:]).ravel()
    f=np.cosh(r/2) if conformal else np.cos(nu*r)
    g=np.cosh(r/2) if conformal else np.cos(nup*r)
    wa=((w*f)[:,None]*(w*f)[None,:]).ravel()/A**2
    wb=((w*g)[:,None]*(w*g)[None,:]).ravel()/B**2
    result=0j
    for i in range(0,len(A),128):
        result+=np.sum(wa[i:i+128,None]*wb[None,:]*kernel(lam/(A[i:i+128,None]*B[None,:])))
    return result

def fixed_bessel(lam,nu,nup,X,Y,order):
    # Independently integrates x dx y dy in logarithmic momenta.
    z,w=leggauss(order);r=-18+22*(z+1)/2;w*=11;x=np.exp(r)
    with mp.workdps(40):
        f=np.array([float((mp.besselk(1j*nu,v*mp.exp(Y))*mp.besselk(1j*nu,v*mp.exp(-Y))).real) for v in x])
        g=np.array([float((mp.besselk(1j*nup,v*mp.exp(X))*mp.besselk(1j*nup,v*mp.exp(-X))).real) for v in x])
    return float(np.sum((w*x*x*f)[:,None]*(w*x*x*g)[None,:]*np.exp(-lam*x[:,None]*x[None,:])))

def main():
    if (OUT/'diagnostics.json').exists():raise RuntimeError('Do not overwrite diagnostics')
    checks=[];details=[];mp.mp.dps=60
    def check(name,error,tol,extra=None):
        checks.append(dict(name=name,error=float(error),tolerance=float(tol),passed=bool(error<tol)))
        if extra is not None:details.append(dict(name=name,data=extra))
    for c in [-100,-2,-1,-.3,-.01,-1/128,-.007,.007,1/128,.01,.1,1,100]:
        c0=mp.mpf(str(c));a=-1/c0
        if c<0:j=a*mp.exp(-a)*(mp.ei(a)-1j*mp.pi)
        else:j=mp.exp(1/c0)*mp.e1(1/c0)/c0
        exact=((1+c0)*j-1)/c0**2
        check('lateral kernel '+str(c),abs(complex(kernel(c))-complex(exact)),2e-11)
        check('conjugate '+str(c),abs(complex(kernel(c,-1))-complex(kernel(c)).conjugate()),1e-14)
    moments=json.loads((OUT/'moments.json').read_text())
    for row in moments:
        n,nu,sigma=row['n'],row['nu'],row['sigma'];b=float(row['value'])
        independent=pair_rapidity(n,nu,sigma)
        check('smeared moment '+str((n,nu,sigma)),abs(independent-b)/b,2e-8,dict(spectral=b,rapidity=independent))
        if n==0:
            nu0=mp.mpf(str(nu));s=mp.mpf(str(sigma))
            def g(x):
                val=mp.pi*nu0/(2*mp.sinh(mp.pi*nu0)) if not x else mp.pi*mp.sin(2*nu0*x)/(2*mp.sinh(mp.pi*nu0)*mp.sinh(2*x))
                return 2*mp.exp(-x*x/(s*s))*val/(s*mp.sqrt(mp.pi))
            alt=mp.quad(g,[0,1,4,mp.inf])
            check('two point '+str((nu,sigma)),abs(alt-mp.mpf(row['value'])),mp.mpf('1e-45'))
    for lam,X,Y in [(.2,.2,-.3),(2.,0.,0.)]:
        vals=[fixed_bessel(lam,.6,.9,X,Y,n) for n in [100,160]]
        rapid=[fixed_rapidity(lam,.6,.9,X,Y,n) for n in [32,48]]
        check('positive full Bessel '+str(lam),abs(vals[-1]-rapid[-1]),3e-10,dict(bessel=vals,rapidity=[[v.real,v.imag] for v in rapid]))
        check('positive Bessel refinement '+str(lam),abs(vals[1]-vals[0]),3e-10)
        check('rapidity refinement '+str(lam),abs(rapid[1]-rapid[0]),3e-8)
    for lam in [.5,-2.]:
        X=.2;Y=-.3;C=np.cosh(X)*np.cosh(Y);c=mp.mpf(str(lam/(4*C)))
        j=mp.exp(1/c)*mp.e1(1/c)/c
        if lam<0:j=mp.conj(j)
        exact=np.pi**2/(16*C)*complex(j)
        val=fixed_rapidity(lam,0,0,X,Y,64,R=20,conformal=True)
        check('conformal normalization '+str(lam),abs(val-exact),3e-7,dict(rapidity=[val.real,val.imag],exact=[exact.real,exact.imag]))
    # Source theta phases and Jacobians: independent direct complex substitution.
    for d in [[.1,-.2,.4,.3],[-.4,.6,.1,-.2]]:
        t=.7;x=.8;y=1.2;s13=(d[0]+d[2])/2;s24=(d[1]+d[3])/2
        theta=np.array([-np.pi/2+1j*(t+d[0]),1j*d[1],np.pi/2+1j*(t+d[2]),np.pi+1j*d[3]])
        pp=-1j*np.exp(-s24)*x;pm=-np.exp(t+s13)*y
        args=np.array([-1j*np.exp(-1j*theta[3])*pp,1j*np.exp(-1j*theta[1])*pp,1j*np.exp(1j*theta[2])*pm,-1j*np.exp(1j*theta[0])*pm])
        X=(d[0]-d[2])/2;Y=(d[1]-d[3])/2
        target=np.array([x*np.exp(-Y),x*np.exp(Y),y*np.exp(X),y*np.exp(-X)])
        check('source Bessel arguments '+str(d),max(abs(args-target)),1e-14)
        pref=np.exp(-1j*(theta[1]+theta[3])+1j*(theta[0]+theta[2]))*(-np.exp(-2*s24))*np.exp(2*t+2*s13)
        check('source prefactor '+str(d),abs(pref-1),1e-14)
        check('source phase '+str(d),abs(1j*pp*pm+np.exp(t+s13-s24)*x*y),1e-14)
    inputs=['code/principal_eikonal.py','code/check_principal_eikonal.py','results/stage4/principal-eikonal/moments.json']
    out=dict(status='passed' if all(r['passed'] for r in checks) else 'failed',checks=checks,details=details,
             scope='Formula, positive-coupling full Bessel integrals and conformal controls; numerical tolerances are empirical, not interval enclosures.',
             hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    (OUT/'diagnostics.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(dict(status=out['status'],count=len(checks),failures=[r for r in checks if not r['passed']]),indent=2))

if __name__=='__main__':main()
