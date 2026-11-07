"""Independent finite coefficient-tensor and charged-probe checks."""
from pathlib import Path
from itertools import permutations
from math import factorial,comb
import hashlib,json,platform
import numpy as np
import scipy
from scipy.linalg import expm
from complex_syk_dictionary import fermions,creators
from bilocal_channel import blocks
from memory_wick_rank import distribution

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-wick-rank'


def tensor_checks():
    rows=[]
    for N,p,m in [(3,1,2),(4,1,2),(4,2,2)]:
        cs=fermions(N);bs=creators(cs,p);ids=[x for x in range(2**N) if x.bit_count()==m]
        d=len(ids);L=len(bs);C=comb(m,p)*comb(N-m+p,p)
        ts=[[ (a@b.conj().T)[np.ix_(ids,ids)]/np.sqrt(C) for b in bs] for a in bs]
        basis=[]
        for i in range(L):
            basis.append(ts[i][i])
            for j in range(i):
                basis.extend([(ts[i][j]+ts[j][i])/np.sqrt(2),1j*(ts[i][j]-ts[j][i])/np.sqrt(2)])
        basis=np.array(basis);I=np.eye(d);cas=np.zeros((d*d,d*d),complex)
        E=[[ (a.conj().T@b)[np.ix_(ids,ids)] for b in cs] for a in cs]
        for i in range(N):
            for j in range(N):
                ca=np.kron(E[i][j],I)-np.kron(I,E[i][j].T)
                cb=np.kron(E[j][i],I)-np.kron(I,E[j][i].T)
                cas+=ca@cb
        ev,U=np.linalg.eigh(cas);rank=(N+1-np.sqrt((N+1)**2-2*ev))/2
        assert np.max(abs(rank-np.rint(rank)))<1e-10
        rank=np.rint(rank)
        for k in [1,2,3]:
            tensor=basis if k==1 else (np.einsum('aij,bjk->abik',basis,basis) if k==2 else np.einsum('aij,bjk,ckl->abcil',basis,basis,basis))
            sym=sum(tensor.transpose(tuple(perm)+(k,k+1)) for perm in permutations(range(k)))/factorial(k)
            independent=(tensor.reshape(-1,d*d)@U.conj())*(k*p-rank)
            coherent=(sym.reshape(-1,d*d)@U.conj())*(k*p-rank)
            raw=float(np.vdot(independent,independent).real/d)
            actual=float(factorial(k)*np.vdot(coherent,coherent).real/d)
            pr=distribution(N,p,m,k);expected=float(sum(pr[r]*(k*p-r)**2 for r in range(len(pr))))
            assert abs(raw-expected)<2e-10*max(1,expected)
            assert actual<=factorial(k)*expected+2e-9
            rows.append(dict(N=N,p=p,m=m,k=k,real_Gaussian_coordinates=len(basis),
                             coefficient_tensor_entries=int(tensor.size),
                             normalized_highest_chaos_deficit_norm_squared=actual,
                             independent_bound=factorial(k)*expected,
                             independent_rank_probability_error=abs(raw-expected)))
    return rows


def probe_checks():
    rng=np.random.default_rng(1107001);rows=[]
    for N,p,m in [(4,1,2),(5,2,3),(6,2,3)]:
        cs=fermions(N);bs=creators(cs,p);L=len(bs)
        z=rng.integers(-2,3,(L,L))+1j*rng.integers(-2,3,(L,L));J=z+z.conj().T
        H=sum(J[a,b]*(A@B.conj().T) for a,A in enumerate(bs) for b,B in enumerate(bs))
        scale=max(1,np.linalg.norm(H,2));H=H/scale
        ids,h,C,K=blocks(N,p,H,cs,m);V=C/np.sqrt(m);a=(1-p/m)*h
        R=K@V-V@a;lower=[x for x in range(2**N) if x.bit_count()==m-1]
        hm=H[np.ix_(lower,lower)];cc=np.array([c[np.ix_(lower,ids)] for c in cs])
        phases=np.exp(2j*np.pi*np.outer(np.arange(N),np.arange(N))/N)
        probes=np.einsum('ki,iab->kab',phases,cc)
        s,u,T=.61,-.37,.29
        def colored(t):
            e=expm(1j*t*hm)
            return [b.conj().T@e@b/m for b in probes]
        Fs=V.conj().T@expm(1j*s*K)@V;Fu=V.conj().T@expm(1j*u*K)@V
        Us=expm(1j*T*h);colors_s=colored(s);colors_u=colored(u)
        four=sum(np.trace(x.conj().T@Us@y@Us.conj().T)/len(h) for x in colors_s for y in colors_u)/(N*N)
        direct=np.trace(Fs.conj().T@Us@Fu@Us.conj().T)/len(h)
        # Analytic endpoint derivatives of e^(iat1) F(t2-t1)e^(-iat2).
        eK=expm(1j*s*K);Fp=1j*V.conj().T@K@eK@V;Fpp=-V.conj().T@K@K@eK@V
        derivative=-Fpp+1j*(a@Fp+Fp@a)+a@Fs@a
        lam=R.conj().T@eK@R
        err=max(abs(four-direct),np.linalg.norm(derivative-lam))
        assert err<2e-12
        rows.append(dict(N=N,p=p,m=m,time_parameters=[s,u,T],fourpoint_absolute_error=float(abs(four-direct)),
                         derivative_matrix_error=float(np.linalg.norm(derivative-lam)),
                         coupling_sha256=hashlib.sha256(J.tobytes()).hexdigest(),Hamiltonian_scale=float(scale)))
    return rows


if __name__=='__main__':
    dest=OUT/'chaos-and-probe-checks.json'
    if dest.exists():raise FileExistsError(dest)
    result=dict(highest_chaos=tensor_checks(),charged_probes=probe_checks(),seed=1107001,
                scope='Complete finite highest-chaos coefficient tensors through degree three and deterministic charged-probe endpoint identities; no asymptotic fit.',
                environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
                input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in
                              ['code/check_memory_wick_interface.py','code/memory_wick_rank.py','code/bilocal_channel.py','code/complex_syk_dictionary.py']})
    dest.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='passed',highest_chaos_tensors=len(result['highest_chaos']),probe_systems=len(result['charged_probes']))))
