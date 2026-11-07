"""Renewal coefficients, centered memory, and an exact force-comparison bound."""
from pathlib import Path
import json,hashlib,platform
import numpy as np
import sympy as sp
import scipy
from scipy.linalg import expm
from complex_syk_dictionary import fermions,creators,comm,sample_hamiltonian
from bilocal_channel import blocks
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-force'

def exact(a):return sp.Matrix([[int(z.real)+sp.I*int(z.imag) for z in row] for row in a])
def rel(a,b):return float(np.linalg.norm(a-b)/max(1.,np.linalg.norm(a),np.linalg.norm(b)))

def integer_renewal():
    rng=np.random.default_rng(110633);rows=[]
    for N,p,m in [(3,1,2),(3,2,2),(4,2,3)]:
        cs=fermions(N);bs=creators(cs,p);L=len(bs)
        x=rng.integers(-1,2,(L,L))+1j*rng.integers(-1,2,(L,L));J=x+x.conj().T
        H=sum(J[a,b]*(A@B.conj().T) for a,A in enumerate(bs) for b,B in enumerate(bs))
        inds,h,C,K=blocks(N,p,H,cs,m)
        A=sum(comm(H,c).conj().T@comm(H,c) for c in cs)
        M3=sum(c.conj().T@comm(H,comm(H,comm(H,c))) for c in cs)
        T=(M3+M3.conj().T)/2
        h,C,K,aa,tt=exact(h),exact(C),exact(K),exact(A[np.ix_(inds,inds)]),exact(T[np.ix_(inds,inds)])
        moments={n:C.conjugate().T*K**n*C for n in range(1,7)}
        beta={};qbar=m*sp.eye(K.rows)-C*C.conjugate().T;checks=[]
        for n in range(1,7):
            beta[n]=m**(n-1)*moments[n]-sum((beta[k]*m**(n-k-1)*moments[n-k] for k in range(1,n)),sp.zeros(h.rows))
            direct=C.conjugate().T*K*(qbar*K)**(n-1)*C
            assert beta[n].expand()==direct.expand()
            checks.append(n)
        JHA=(h*aa+aa*h)/2
        b3=beta[3]/m**3
        r=sp.Rational(p,m)
        formula=tt/m+(1+2*r)*JHA/m-r*r*(1+r)*h**3
        assert (b3-formula).expand()==sp.zeros(h.rows)
        G=beta[2]/m**2
        centered=b3-(h*G+G*h)/2
        assert (m*centered-tt-2*r*JHA+m*r**3*h**3).expand()==sp.zeros(h.rows)
        rows.append(dict(N=N,p=p,m=m,renewal_orders=checks,exact_first_memory_moment=True,exact_centered_first_moment=True))
    return rows

def bounds():
    rows=[]
    for N,p,m in [(3,1,2),(4,2,3),(5,2,3)]:
        cs,J,H=sample_hamiltonian(N,p,110640+N)
        _,h,C,K=blocks(N,p,H,cs,m);V=C/np.sqrt(m)
        Q=np.eye(len(K))-V@V.conj().T;R=Q@K@V;G=R.conj().T@R
        D=K@V-V@h;Kq=Q@K@Q;dim=len(h)
        g2=float(np.trace(G@G).real/dim);gn=float(np.trace(G).real/dim);hn=float(np.trace(h@h).real/dim)
        first_bound=2*(p/m)*np.sqrt(gn*hn)+(p/m)**2*hn
        values=[]
        for t in [0,.2,.7,1.3]:
            U=expm(1j*t*K);ext=expm(-1j*t*h)
            mem=np.trace(R.conj().T@expm(1j*t*Kq)@R@ext)/dim
            full=np.trace(R.conj().T@U@R@ext)/dim
            force=np.trace(D.conj().T@U@D@ext)/dim
            orth_bound=t*t*g2/2
            assert abs(mem-full)<=orth_bound+2e-12
            assert abs(full-force)<=first_bound+2e-12
            values.append(dict(time=t,memory=[mem.real,mem.imag],force=[force.real,force.imag],
              orthogonal_difference=float(abs(mem-full)),orthogonal_bound=orth_bound,
              parallel_difference=float(abs(full-force)),parallel_bound=float(first_bound)))
        rows.append(dict(N=N,p=p,m=m,mean_G=gn,mean_G_squared=g2,values=values))
    return rows

def main():
    OUT.mkdir(parents=True,exist_ok=True);target=OUT/'bridge.json'
    if target.exists():raise FileExistsError(target)
    record=dict(status='passed',integer_renewal=integer_renewal(),comparison_bounds=bounds(),
      environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=sp.__version__,scipy=scipy.__version__),
      error_scope='Exact noncommuting matrix renewal through order six, and finite-system float tests of analytic scalar memory-force bounds. No fit of finite systems to the DS kernel.',
      input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['code/memory_force_bridge.py','code/bilocal_channel.py','code/complex_syk_dictionary.py']})
    target.write_text(json.dumps(record,indent=2)+'\n')
    print('Exact renewal cases',len(record['integer_renewal']),'scalar bound systems',len(record['comparison_bounds']))

if __name__=='__main__':main()
