"""Exact charge-sector bilocal compression and finite-time memory checks.

No overwrite. Matrix identities use Gaussian integers; propagation uses
independent full-H and charge-sector diagonalizations. Js=1 for DS tables.
"""
from pathlib import Path
from math import comb
from fractions import Fraction
import json, hashlib, platform
import numpy as np
import scipy
from scipy.linalg import expm
from numpy.polynomial.legendre import leggauss
import mpmath as mp
from complex_syk_dictionary import fermions, creators, comm

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-bilocal-channel'

def rel(a,b):return float(np.linalg.norm(a-b)/max(1.,np.linalg.norm(a),np.linalg.norm(b)))
def norm(a):return float(np.vdot(a,a).real/a.shape[1])
def integer_zero(a):
    assert np.max(abs(a))<2**48
    assert np.array_equal(a,np.zeros_like(a))

def blocks(N,p,H,cs,m):
    i=np.array([s for s in range(2**N) if s.bit_count()==m])
    j=np.array([s for s in range(2**N) if s.bit_count()==m-1])
    h=H[np.ix_(i,i)];hm=H[np.ix_(j,j)]
    C=np.vstack([c[np.ix_(j,i)] for c in cs])
    K=np.kron(np.eye(N),hm)
    return i,h,C,K

def deterministic():
    rng=np.random.default_rng(110631);rows=[]
    for N,p,m in [(3,1,2),(3,2,2),(4,2,3),(5,2,3)]:
        cs=fermions(N);bs=creators(cs,p);L=len(bs)
        raw=rng.integers(-2,3,(L,L))+1j*rng.integers(-2,3,(L,L));J=raw+raw.conj().T
        H=sum(J[a,b]*(A@B.conj().T) for a,A in enumerate(bs) for b,B in enumerate(bs))
        A=sum(comm(H,c).conj().T@comm(H,c) for c in cs)
        inds,h,C,K=blocks(N,p,H,cs,m);aa=A[np.ix_(inds,inds)]
        integer_zero(C.conj().T@C-m*np.eye(len(h)))
        integer_zero(C.conj().T@K@C-(m-p)*h)
        Z=m*K@C-(m-p)*C@h
        integer_zero(Z.conj().T@Z-m*m*aa+m*p*p*h@h)
        # Symmetric-bilocal second Taylor coefficient, without finite differences.
        Th2=C.conj().T@K@K@C
        integer_zero(Th2-(m-2*p)*h@h-aa)
        second=-Th2/2+(m-p)*h@h-m*h@h/2
        integer_zero(2*second+aa)
        # Scale this deterministic Hamiltonian for stable finite-time tests.
        scale=max(1.,np.linalg.norm(H,2));H=H/scale;h=h/scale;K=K/scale;aa=aa/scale**2
        V=C/np.sqrt(m);P=V@V.conj().T;Qc=np.eye(len(K))-P
        he=(1-p/m)*h;R=Qc@K@V;G=R.conj().T@R
        Kq=Qc@K@Qc;F=lambda t:V.conj().T@expm(1j*t*K)@V
        W=lambda t:expm(1j*t*K)@V@expm(-1j*t*h)
        checks={}
        checks['G_formula']=rel(G,aa/m-p*p*h@h/m**2)
        checks['G_commutator']=rel(comm(h,G),comm(h,aa)/m)
        minG=float(np.linalg.eigvalsh(G).min());assert minG>-1e-13
        z=.7+.9j
        full=V.conj().T@np.linalg.inv(z*np.eye(len(K))-K)@V
        sig=R.conj().T@np.linalg.inv(z*np.eye(len(K))-Kq)@R
        reduced=np.linalg.inv(z*np.eye(len(h))-he-sig)
        checks['resolvent']=rel(full,reduced)
        # Full Fock-space definition, independent of the dilation formula.
        def S(a,b):
            ua=expm(1j*a*H);ub=expm(1j*b*H)
            out=sum((ua@c@ua.conj().T).conj().T@(ub@c@ub.conj().T) for c in cs)
            return out[np.ix_(inds,inds)]/m
        ts=[-.31,.17,.59]
        ws=[W(t) for t in ts]
        gram=np.block([[a.conj().T@b for b in ws] for a in ws])
        assert np.linalg.eigvalsh(gram).min()>-2e-12
        checks['bilocal']=max(rel(S(a,b),W(a).conj().T@W(b)) for a,b in [(-.31,.17),(.17,.59)])
        a,b,d=ws
        defect=a.conj().T@(np.eye(len(K))-b@b.conj().T)@d
        checks['composition']=rel((a.conj().T@b)@(b.conj().T@d),a.conj().T@d-defect)
        quad=[]
        t=.73;lhs=1j*V.conj().T@K@expm(1j*t*K)@V
        for count in [24,48]:
            x,w=leggauss(count);x=(x+1)*t/2;w=w*t/2
            integ=sum(wi*(R.conj().T@expm(1j*(t-si)*Kq)@R@F(si)) for si,wi in zip(x,w))
            quad.append(rel(lhs,1j*he@F(t)-integ))
        checks['memory_quadrature']=quad
        finite=[]
        for t in [.2,.7,1.4]:
            diff=norm(expm(1j*t*K)@V-V@expm(1j*t*he))
            bound=t*t*np.trace(G).real/len(h)
            assert diff<=bound+2e-13
            collapsed=norm(F(t)@expm(-1j*t*h)-np.eye(len(h)))
            collapse_bound=2*t*t*(np.trace(G).real/len(h)+(p/m)**2*np.trace(h@h).real/len(h))
            assert collapsed<=collapse_bound+2e-13
            finite.append(dict(time=t,dilation_error_squared=diff,bound=bound,bilocal_error_squared=collapsed,identity_bound=collapse_bound))
        assert max(v for k,v in checks.items() if k!='memory_quadrature')<2e-12
        assert max(quad)<2e-12
        rows.append(dict(N=N,p=p,m=m,scale=scale,exact_identities=5,checks=checks,min_G_eigenvalue=minG,
           memory_initial_norm=norm(G),finite_time=finite))
    return rows

def covariance_checks():
    rows=[]
    for N,p,m in [(3,1,2),(3,2,2),(4,2,3),(4,2,4),(5,2,3)]:
        cs=fermions(N);bs=creators(cs,p);value=Fraction(0)
        for A in bs:
            for B in bs:
                H=A@B.conj().T
                _,h,C,K=blocks(N,p,H,cs,m)
                Z=m*K@C-(m-p)*C@h
                zz=np.vdot(Z,Z).real
                assert zz==round(zz)
                value+=Fraction(round(zz),m**3*len(h))
        C2=comb(m,p)*comb(N-m+p,p)
        target=Fraction(C2*p*(m-p)*(N+1),m*m*(N-m+1))
        assert value==target,(N,p,m,value,target)
        rows.append(dict(N=N,p=p,m=m,monomials=len(bs)**2,E_trace_G_over_w=str(value),formula=str(target)))
    return rows

def double_scale():
    mp.mp.dps=80;rows=[]
    for lam in ['0.25','0.5','1']:
        l=mp.mpf(lam)
        for p in [8,16,32,64,128,256]:
            N=int(p*p/l);m=N//2
            w=mp.mpf(4)**p/l*(mp.factorial(p)/mp.mpf(N)**p)**2
            C=w*mp.binomial(m,p)*mp.binomial(N-m+p,p)
            G=C*p*(m-p)*(N+1)/(mp.mpf(m)**2*(N-m+1))
            identity_bound=2*(G+(mp.mpf(p)/m)**2*C)
            rows.append(dict(N=N,p=p,m=m,lambda_value=lam,C=str(C),p_E_trace_G=str(p*G),
               limit='4',identity_bound_at_t1=str(identity_bound)))
    return rows

def main():
    OUT.mkdir(parents=True,exist_ok=True);target=OUT/'checks.json'
    if target.exists():raise FileExistsError(target)
    r=dict(status='passed',seed=110631,deterministic=deterministic(),covariance=covariance_checks(),double_scaling=double_scale(),
       environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__),
       error_scope='Exact Gaussian-integer identities and monomial covariance; finite-dimensional float propagation and quadrature, not interval certification. DS table evaluates exact sector formulas at 80 digits. Uniform finite-time conclusion follows from the analytic dilation bound.',
       input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['code/bilocal_channel.py','code/complex_syk_dictionary.py']})
    target.write_text(json.dumps(r,indent=2)+'\n')
    print('passed:',len(r['deterministic']),'matrix systems,',len(r['covariance']),'exact covariances,',len(r['double_scaling']),'DS rows')
    for row in r['double_scaling']:
        if row['p']==256:print(row['lambda_value'],row['p_E_trace_G'],row['identity_bound_at_t1'])

if __name__=='__main__':main()
