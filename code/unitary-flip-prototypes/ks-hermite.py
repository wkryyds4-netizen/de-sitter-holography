import numpy as np,json,time
from scipy.linalg import eigh_tridiagonal
from importlib import import_module
m=import_module('ks-prototype')
def run(N,h,flat=False):
    z,Q=eigh_tridiagonal(np.zeros(N),np.sqrt(np.arange(1,N)/2))
    # Q[n,j] is the nth normalized oscillator polynomial at node j, times sqrt(weight).
    T=(Q.T*((-1j)**np.arange(N)))@Q
    t=np.sqrt(np.pi)*z;x=t[:,None];y=t[None,:]
    fun=(lambda x,y,h:np.exp(-1j*y*np.logaddexp(0,x)/np.pi)) if flat else m.phase
    a=fun(x,y,h);b=fun(y,x,h);chirp=np.exp(1j*x*y/np.pi)
    def B(f):return T.conj().T@(b*(T@f@T.T))@T.conj()
    f=Q[0,:,None]*Q[0,None,:]
    ab=a*B(f);ba=B(a*f);rhs=B(np.conj(chirp)*B(chirp*a*f));d=ab-rhs
    ov=np.vdot(ab,ba)
    return dict(N=N,h=h,flat=flat,pentagon=float(np.linalg.norm(d)),projected=[np.vdot(f,d).real,np.vdot(f,d).imag],overlap=[ov.real,ov.imag],D=float(np.linalg.norm(ab-ba)**2),unitarity=float(np.max(abs(T.conj().T@T-np.eye(N)))))
if __name__=='__main__':
    for N in [64,128,256,512]:
        for h,flat in [(1,False),(.5,False),(1,True)]: print(json.dumps(run(N,h,flat)),flush=True)
