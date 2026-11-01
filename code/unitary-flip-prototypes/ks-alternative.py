import sys,json,time,numpy as np
from importlib import import_module
from scipy.fft import fft,ifft,fftfreq
m=import_module('ks-prototype')
def run(N,L,h):
    t=np.linspace(-L,L,N,endpoint=False);p=2*np.pi**2*fftfreq(N,d=2*L/N)
    i=np.arange(N)[:,None];j=np.arange(N)[None,:]
    ix=(i+j-N//2)%N; inv=(j-i+N//2)%N
    f1=m.phase(-t[:,None],-p[None,:],h)
    f2=m.phase(t[None,:],-p[:,None],h)
    fs=m.phase(t[None,:],-p[:,None],h)
    def A(f):return ifft(f1*fft(f,axis=1,norm='ortho'),axis=1,norm='ortho')
    def B(f):return ifft(f2*fft(f,axis=0,norm='ortho'),axis=0,norm='ortho')
    def C(f):
        g=f[i,ix];g=ifft(fs*fft(g,axis=0,norm='ortho'),axis=0,norm='ortho')
        return g[i,inv]
    X=t[:,None];Y=t[None,:]
    f=np.exp(-((X+.3)**2+(Y-.4)**2)/2+.2j*X-.15j*Y);f/=np.linalg.norm(f)
    lhs=A(B(f));rhs=B(C(A(f)));d=lhs-rhs
    return dict(N=N,L=L,h=h,pentagon=float(np.linalg.norm(d)),residual_central_norm=float(np.linalg.norm(d[(abs(X)+abs(Y))<L/2])),inner=[np.vdot(lhs,rhs).real,np.vdot(lhs,rhs).imag],lhs_tail=float(np.linalg.norm(lhs[(abs(X)+abs(Y))>L/2])),rhs_tail=float(np.linalg.norm(rhs[(abs(X)+abs(Y))>L/2])))
if __name__=='__main__':
    for N,L,h in [(256,16,.5),(512,24,.5),(256,16,1),(512,24,1)]:print(json.dumps(run(N,L,h)),flush=True)
