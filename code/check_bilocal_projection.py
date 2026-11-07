"""Independent unequal-time four-term algebra and charged spectral reconstruction."""
from pathlib import Path
import json,hashlib,platform
import numpy as np
from scipy.linalg import expm
import scipy
from complex_syk_dictionary import sample_hamiltonian, comm

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-bilocal-channel'
def hs(a):return float(np.vdot(a,a).real/len(a))
def rel(a,b):return float(np.linalg.norm(a-b)/max(1.,np.linalg.norm(a),np.linalg.norm(b)))

def main():
    target=OUT/'projection-checks.json'
    if target.exists():raise FileExistsError(target)
    rows=[]
    for N,p,m in [(4,2,3),(5,2,3)]:
        seed=110632+N;cs,J,H=sample_hamiltonian(N,p,seed)
        inds=np.array([s for s in range(2**N) if s.bit_count()==m])
        prev=np.array([s for s in range(2**N) if s.bit_count()==m-1])
        take=lambda a:a[np.ix_(inds,inds)]
        times=[-.47,.13,.71,-.22]
        ct=[]
        for t in times:
            u=expm(1j*t*H);ct.append([u@c@u.conj().T for c in cs])
        s12=sum(a.conj().T@b for a,b in zip(ct[0],ct[1]))
        s34=sum(a.conj().T@b for a,b in zip(ct[2],ct[3]))
        parts=[np.zeros_like(H) for _ in range(4)]
        anti=lambda a,b:a@b+b@a
        for i in range(N):
            for j in range(N):
                a,b,c,d=ct[0][i].conj().T,ct[1][i],ct[2][j].conj().T,ct[3][j]
                terms=[a@anti(b,c)@d,-c@anti(d,a)@b,c@a@anti(d,b),-anti(a,c)@b@d]
                for k in range(4):parts[k]+=terms[k]
        lhs=take(comm(s12,s34));full=take(sum(parts))
        algebra_error=rel(lhs,full);assert algebra_error<2e-12
        # Spectral bridge: keep coherent intermediate-state amplitudes.
        h=take(H);hm=H[np.ix_(prev,prev)]
        e,U=np.linalg.eigh(h);em,Um=np.linalg.eigh(hm)
        maps=[Um.conj().T@c[np.ix_(prev,inds)]@U for c in cs]
        deltas=[-.6,.2,.8];ops=[];spectral=[]
        for delta in deltas:
            up=expm(1j*delta*H/2);un=up.conj().T
            center=sum((un@c@up).conj().T@(up@c@un) for c in cs)
            ops.append(take(center))
            mid=sum(c.conj().T@np.diag(np.exp(1j*delta*em))@c for c in maps)
            phase=np.exp(-.5j*delta*e)
            spectral.append(U@(phase[:,None]*mid*phase[None,:])@U.conj().T)
        spectral_error=max(rel(a,b) for a,b in zip(ops,spectral));assert spectral_error<2e-12
        gram=np.array([[np.trace(a.conj().T@b)/len(h) for b in ops] for a in ops])
        min_gram=float(np.linalg.eigvalsh(gram).min());assert min_gram>-2e-12
        # Even-relative-time jet: an independent finite-difference convergence check.
        A=sum(comm(H,c).conj().T@comm(H,c) for c in cs);aa=take(A);errs=[]
        for delta in [.04,.02,.01]:
            u=expm(1j*delta*H/2)
            S=sum((u.conj().T@c@u).conj().T@(u@c@u.conj().T) for c in cs)
            reconstructed=2*(m*np.eye(len(h))-take((S+S.conj().T)/2))/delta**2
            errs.append(rel(reconstructed,aa))
        assert errs[2]<errs[1]/3 and errs[1]<errs[0]/3
        rows.append(dict(N=N,p=p,m=m,seed=seed,times=times,
           commutator_HS_squared=hs(lhs),four_term_HS_squared=[hs(take(x)) for x in parts],
           four_term_error=algebra_error,coherent_spectral_error=spectral_error,
           minimum_covariance_gram_eigenvalue=min_gram,jet_steps=[.04,.02,.01],jet_relative_errors=errs))
    record=dict(status='passed',records=rows,
       error_scope='Finite-system floating-point check of full projected four-term algebra, coherent spectral reconstruction and quadratic jet convergence. Does not extrapolate projected commutators to DS.',
       environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
       input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['code/check_bilocal_projection.py','code/complex_syk_dictionary.py']})
    target.write_text(json.dumps(record,indent=2)+'\n')
    for r in rows:print(r['N'],r['commutator_HS_squared'],r['four_term_HS_squared'],r['jet_relative_errors'])

if __name__=='__main__':main()
