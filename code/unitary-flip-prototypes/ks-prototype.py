import numpy as np
from scipy.fft import fft,ifft,fftfreq
import time,json

def logone(a,b):
    # Principal log(1+exp(a+ib)); avoid overflow without changing its argument.
    e=np.exp(-np.abs(a))
    z=np.where(a>0,np.exp(1j*b)+e,1+e*np.exp(1j*b))
    return np.maximum(a,0)+np.log(z)

def phase(x,y,h):
    x,y=np.broadcast_arrays(x,y);s=np.zeros(x.shape)
    for a,b,rate,sign in [(x,h*y,np.pi*h,-1),(y,-x/h,np.pi/h,1)]:
        n=max(1,int(np.ceil((float(np.max(a))+38)/ (2*rate))))
        for j in range(n):s+=sign*logone(a-(2*j+1)*rate,b).imag
    return np.exp(2j*s)

def run(N,L,h,sh=(0,0,0,0)):
    t=np.linspace(-L,L,N,endpoint=False);p=2*np.pi**2*fftfreq(N,d=2*L/N)
    a,b,c,d=sh
    f1=phase(p[:,None]+a,t[None,:]+b,h)
    f2=phase(p[None,:]+c,-t[:,None]+d,h)
    fs=phase(p[:,None]+a+c,t[None,:]+b+d,h)
    i=np.arange(N)[:,None];j=np.arange(N)[None,:]
    ix=(i+j-N//2)%N; inv=(j-i+N//2)%N
    def A(f):return ifft(f1*fft(f,axis=0,norm='ortho'),axis=0,norm='ortho')
    def B(f):return ifft(f2*fft(f,axis=1,norm='ortho'),axis=1,norm='ortho')
    def C(f):
        g=f[i,ix];g=ifft(fs*fft(g,axis=0,norm='ortho'),axis=0,norm='ortho')
        return g[i,inv]
    X=t[:,None];Y=t[None,:]
    f=np.exp(-((X+.3)**2+(Y-.4)**2)/2+.2j*X-.15j*Y);f/=np.linalg.norm(f)
    ab=A(B(f));ba=B(A(f));rhs=B(C(A(f)))
    z=np.vdot(ab,ba)
    return dict(N=N,L=L,h=h,sh=sh,pentagon=float(np.linalg.norm(ab-rhs)),overlap=[z.real,z.imag],D=float(np.linalg.norm(ab-ba)**2),norms=[float(np.linalg.norm(v)) for v in [ab,ba,rhs]])

if __name__=='__main__':
    for N,L in [(128,8),(256,8),(256,12),(512,12)]:
        start=time.time();r=run(N,L,.5);r['seconds']=time.time()-start;print(json.dumps(r),flush=True)
