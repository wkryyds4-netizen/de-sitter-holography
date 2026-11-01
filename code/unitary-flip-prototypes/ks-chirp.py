import numpy as np,json,time
from scipy.fft import fft2,ifft2,fftfreq
from importlib import import_module
m=import_module('ks-prototype')
def run(N,L,h,flat=False):
    t=np.linspace(-L,L,N,endpoint=False);p=2*np.pi**2*fftfreq(N,d=2*L/N)
    x=t[:,None];y=t[None,:]
    fun=(lambda x,y,h:np.exp(-1j*y*np.logaddexp(0,x)/np.pi)) if flat else m.phase
    a=fun(x,y,h);b=fun(p[None,:],p[:,None],h);chirp=np.exp(1j*x*y/np.pi)
    def B(f):return ifft2(b*fft2(f,norm='ortho'),norm='ortho')
    f=np.exp(-((x+.3)**2+(y-.4)**2)/2+.2j*x-.15j*y);f/=np.linalg.norm(f)
    ab=a*B(f);ba=B(a*f);rhs=B(np.conj(chirp)*B(chirp*a*f))
    d=ab-rhs;ov=np.vdot(ab,ba)
    return dict(N=N,L=L,h=h,flat=flat,pentagon=float(np.linalg.norm(d)),central=float(np.linalg.norm(d[(abs(x)+abs(y))<L/2])),overlap=[ov.real,ov.imag],D=float(np.linalg.norm(ab-ba)**2))
if __name__=='__main__':
    for N,L,h,flat in [(256,16,.5,False),(512,24,.5,False),(256,16,1,False),(512,24,1,False),(256,16,1,True)]:print(json.dumps(run(N,L,h,flat)),flush=True)
