"""Exact covariance algebra, finite matrix loop identity, and Mellin slopes."""
import hashlib,itertools,json
from pathlib import Path
import mpmath as mp
import numpy as np
import sympy as sp
from numpy.polynomial.hermite import hermgauss
from mass_filter import covariance,slope
from principal_eikonal import spectral_moment
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/stage4/mass-filter'

def main():
    if (OUT/'diagnostics.json').exists():raise RuntimeError('Do not overwrite')
    checks=[];details=[]
    def ck(name,error,tol):checks.append(dict(name=name,error=float(error),tolerance=float(tol),passed=bool(error<tol)))
    C=sp.Matrix([[5,-1,-3,-1],[-1,5,-1,-3],[-3,-1,5,-1],[-1,-3,-1,5]])/8
    D=sp.eye(4)-sp.Matrix([[0,1,0,0],[0,0,1,0],[0,0,0,1],[1,0,0,0]])
    P=sp.eye(4)-sp.ones(4)/4
    H=sp.Matrix([[1,0,-1,0],[0,1,0,-1],[1,-1,1,-1]])/2
    for name,res in [('characteristic',D.T*C*D-2*P),('reduced covariance',H*C*H.T-sp.eye(3)/2),('zero center',C*sp.ones(4,1)),('pseudoinverse',(D*D.T)*(C/2)-P)]:
        ck(name,0 if res==sp.zeros(*res.shape) else 1,1e-30)
    rng=np.random.default_rng(20261031);m=np.array([.3,.7,1.2]);M=.8;t=.4
    v=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3));v=(v+v.conj().T)/2
    w=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3));w=(w+w.conj().T)/2
    for delta in [.2,.6,1.]:
        y=np.diag(np.exp(np.pi*m/2-delta**2*(m-M)**2));vt=v*np.exp(1j*t*(m[:,None]-m[None,:]))
        direct=np.trace(y@vt@y@w@y@vt@y@w);vals=[];meanonly=0j
        for order in [8,16]:
            z,hw=hermgauss(order);hw/=np.sqrt(np.pi);total=0j
            for i,j,k,l in itertools.product(range(3),repeat=4):
                mm=m[[i,j,k,l]];mu=mm.mean();omega=mm-np.roll(mm,-1)
                amplitude=v[i,j]*w[j,k]*v[k,l]*w[l,i]*np.exp(1j*t*(omega[0]+omega[2]))
                weight=np.exp(2*np.pi*mu-4*delta**2*(mu-M)**2)
                coefs=[omega[0]-omega[2],omega[1]-omega[3],(omega[0]-omega[1]+omega[2]-omega[3])/2]
                char=np.prod([np.sum(hw*np.exp(1j*delta*z*c)) for c in coefs])
                total+=weight*amplitude*char
                if order==8:meanonly+=weight*amplitude
            vals.append(total)
        scale=max(abs(direct),1.)
        ck('finite matrix Fourier '+str(delta),abs(vals[-1]-direct)/scale,2e-13)
        ck('finite matrix quadrature '+str(delta),abs(vals[0]-vals[1])/scale,3e-6)
        ck('finite matrix trace real '+str(delta),abs(direct.imag)/scale,2e-14)
        details.append(dict(delta=delta,direct=[direct.real,direct.imag],quadrature=[[r.real,r.imag] for r in vals],mean_only=[meanonly.real,meanonly.imag]))
    mp.mp.dps=50;slopes=[]
    for nu in [mp.mpf('.6'),mp.mpf('.9')]:
        for n in [0,1,2,3]:
            a=mp.mpf(n+2)/2
            exact=(a*a+nu*nu)/(2*a+1)
            In=2**(n-1)*mp.gamma(a)**2*abs(mp.gamma(a+1j*nu))**2/mp.gamma(n+2)
            def density(o):
                b=mp.mpf(n+2)/4
                return o*o*abs(2**(mp.mpf(n)/2-1-1j*o)*mp.gamma(b+1j*(nu-o)/2)*mp.gamma(b-1j*(nu+o)/2))**2/mp.pi
            numerical=mp.quad(density,[0,1,4,12,mp.inf])/In
            ck('spectral derivative '+str((n,nu)),abs(numerical-exact),mp.mpf('1e-43'))
            slopes.append(dict(n=n,nu=str(nu),exact=str(exact),spectral=str(numerical)))
    for n in [1,2]:
        values=[]
        for delta in [.02,.01]:
            factors=[]
            for nu in [.6,.9]:
                bd=spectral_moment(n,nu,delta,50)/spectral_moment(0,nu,delta,50)
                b0=spectral_moment(n,nu,0,50)/spectral_moment(0,nu,0,50)
                factors.append(bd/b0)
            enhancement=mp.exp(n*n*mp.mpf(str(delta))**2/4)*factors[0]*factors[1]
            values.append(float((enhancement-1)/delta**2))
        prediction=slope(n,.6,.9)
        ck('small filter slope '+str(n),abs(values[-1]-prediction),.001)
        details.append(dict(n=n,small_delta_slopes=values,prediction=prediction,independent_prediction=slope(n,.6,.9,1.)))
    inputs=['code/check_mass_filter.py','code/mass_filter.py','code/principal_eikonal.py']
    out=dict(status='passed' if all(r['passed'] for r in checks) else 'failed',checks=checks,details=details,spectral_slopes=slopes,
             scope='Exact cycle algebra and finite matrix identity; conditional moment derivatives. No gravitational Hilbert identification.',
             hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    (OUT/'diagnostics.json').write_text(json.dumps(out,indent=2));print(json.dumps(dict(status=out['status'],count=len(checks),failures=[c for c in checks if not c['passed']])))
if __name__=='__main__':main()
