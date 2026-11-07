"""Finite CAR checks of the finite-return expansion; no asymptotic fit.

Creates fresh evidence only. All propagators share the same Hamiltonian.
"""
from pathlib import Path
from math import comb, factorial
from itertools import product
import hashlib
import json
import platform
import numpy as np
import scipy
from scipy.linalg import expm, eigh
from complex_syk_dictionary import fermions, creators, comm
from bilocal_channel import blocks

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'results/susskind-memory-return-expansion'

def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0

def hs(a):
    return float(np.linalg.norm(a) / np.sqrt(a.shape[1]))

def discrepancy(a, b):
    return float(np.linalg.norm(a-b) / max(1., np.linalg.norm(a), np.linalg.norm(b)))

def force_counts():
    """Exact integer covariance maps, including overlapping creation/annihilation."""
    rows = []
    for N, p in [(3,1), (3,2), (4,2), (5,2), (5,3)]:
        cs = fermions(N); bs = creators(cs,p); dim = 2**N
        maxerr = 0.; checks = 0
        for i, c in enumerate(cs):
            vin = np.zeros((dim,dim), complex); vout = vin.copy()
            for a in bs:
                for b in bs:
                    f = comm(a @ b.conj().T, c)
                    vin += f.conj().T @ f
                    vout += f @ f.conj().T
            ref_in = []; ref_out = []
            for state in range(dim):
                m = state.bit_count(); eta = (state >> i) & 1
                f0 = choose(m,p)*choose(N-m+p-1,p-1)
                f1 = (choose(m-1,p)*choose(N-m+p,p-1)
                      + choose(m-1,p-1)*choose(N-m+p-1,p-1))
                ref_in.append(f1 if eta else f0)
                ref_out.append(choose(m-eta,p-1)*choose(N-m+p-1,p))
            e = max(np.max(abs(vin-np.diag(ref_in))), np.max(abs(vout-np.diag(ref_out))))
            assert e == 0
            maxerr = max(maxerr,float(e)); checks += 2*dim
        rows.append(dict(N=N,p=p,diagonal_entries_checked=checks,max_integer_error=maxerr))
    # Rational bound before the exponential simplification; no floating arithmetic.
    bound_checks = 0
    for N in range(2,71):
        for p in range(1,min(N,12)+1):
            # count <= (N+p)^ (2p-1) / (4^(p-1) p! (p-1)!).
            rhs = (N+p)**(2*p-1)
            den = 4**(p-1)*factorial(p)*factorial(p-1)
            for m in range(N+1):
                for eta in [0,1]:
                    if eta > m or m-eta > N-1: continue
                    a = (choose(m,p)*choose(N-m+p-1,p-1) if not eta else
                         choose(m-1,p)*choose(N-m+p,p-1)+choose(m-1,p-1)*choose(N-m+p-1,p-1))
                    b = choose(m-eta,p-1)*choose(N-m+p-1,p)
                    assert max(a,b)*den <= rhs
                    bound_checks += 2
    return dict(exact_CAR=rows,rational_bound_checks=bound_checks)

def cascade(K, a, R, s):
    dims = [len(K),len(a),len(K),len(a),len(K)]
    edges = [R,R.conj().T,R,R.conj().T]
    offsets = np.cumsum([0]+dims)
    B = np.zeros((offsets[-1],offsets[-1]),complex)
    for j, d in enumerate(dims):
        sl = slice(offsets[j],offsets[j+1])
        B[sl,sl] = 1j*(K if j%2 == 0 else a)
        if j < 4: B[sl,offsets[j+1]:offsets[j+2]] = edges[j]
    E = expm(s*B)
    return [R.conj().T @ E[:len(K),offsets[2*n]:offsets[2*n+1]] @ R for n in range(3)]

def simplex(s, nodes):
    b = np.diag(1j*np.array(nodes)) + np.diag(np.ones(len(nodes)-1),1)
    return expm(s*b)[0,-1]

def spectral_paths(h, lower, cs, inds, lowerinds, p, m, s):
    e, u = eigh(h); eps, v = eigh(lower); ae = (1-p/m)*e; d=len(e)
    cc = np.array([v.conj().T @ c[np.ix_(lowerinds,inds)] @ u for c in cs])
    W = np.zeros((len(eps),d,d),complex)
    for al in range(len(eps)):
        z = cc[:,al,:]*(eps[al]-ae)[None,:]
        W[al] = z.conj().T @ z / m
    out = [np.einsum('a,aij->ij',np.exp(1j*eps*s),W), np.zeros((d,d),complex), np.zeros((d,d),complex)]
    for al,be,c in product(range(len(eps)),range(len(eps)),range(d)):
        out[1] += np.outer(W[al,:,c],W[be,c,:])*simplex(s,[eps[al],ae[c],eps[be]])
        for ga,f in product(range(len(eps)),range(d)):
            out[2] += np.outer(W[al,:,c],W[ga,f,:])*W[be,c,f]*simplex(s,[eps[al],ae[c],eps[be],ae[f],eps[ga]])
    return [u @ x @ u.conj().T for x in out]

def finite_checks():
    rows=[]; rng=np.random.default_rng(1106701)
    for N,p,m in [(3,1,2),(4,2,3),(5,2,3),(6,2,3)]:
        cs=fermions(N); bs=creators(cs,p); L=len(bs)
        raw=rng.integers(-2,3,(L,L))+1j*rng.integers(-2,3,(L,L)); J=raw+raw.conj().T
        H=sum(J[a,b]*(A@B.conj().T) for a,A in enumerate(bs) for b,B in enumerate(bs))
        scale=max(1.,np.linalg.norm(H,2)); H=H/scale
        inds,h,C,K=blocks(N,p,H,cs,m); V=C/np.sqrt(m); Q=np.eye(len(K))-V@V.conj().T
        a=(1-p/m)*h; R=Q@K@V; D=Q@K@Q; G=R.conj().T@R
        e,u=eigh(h); gaps=abs(e[:,None]-e[None,:]); mask=gaps>1e-10
        assert np.min(gaps[np.triu_indices(len(e),1)])>1e-8
        def pinch_dynamic(x):
            return u @ ((u.conj().T@x@u)*mask) @ u.conj().T
        def field(x,s):
            z=expm(-.5j*s*h)
            return pinch_dynamic(m*z@x@z)
        times=[-.8,.35,1.1,2.0]; data=[]; allkernels={}
        for s in times:
            M=R.conj().T@expm(1j*s*D)@R
            js=cascade(K,a,R,s); sums=np.cumsum(js,axis=0)
            E=expm(1j*s*K); F=V.conj().T@E@V
            Fp=1j*V.conj().T@K@E@V; Fpp=-V.conj().T@K@K@E@V
            derivative=-Fpp+1j*(a@Fp+Fp@a)+a@F@a
            dererr=discrepancy(js[0],derivative); assert dererr<2e-12
            checks=[]
            ev=np.maximum(np.linalg.eigvalsh(G),0)
            for j in range(3):
                err=hs(M-sums[j]); bound=abs(s)**(2*j+2)/factorial(2*j+2)*np.sqrt(np.mean(ev**(2*j+4)))
                assert err <= bound+2e-13
                checks.append(dict(returns=j,error=err,bound=float(bound),ratio=float(err/bound) if bound else 0.))
            adj=max(discrepancy(x.conj().T,y) for x,y in zip(js,cascade(K,a,R,-s))); assert adj<2e-12
            data.append(dict(s=s,remainder=checks,derivative_error=dererr,adjoint_error=adj))
            allkernels[s]=(M,sums)
        s,u_time,T=-.8,1.1,.67
        M1,L1=allkernels[s]; M2,L2=allkernels[u_time]
        F1,F2=field(M1,s),field(M2,u_time); X1,X2=field(L1[2],s),field(L2[2],u_time)
        U=expm(1j*T*h)
        cov=lambda x,y:np.trace(x.conj().T@U@y@U.conj().T)/len(h)
        aa=u.conj().T@L1[2]@u; bb=u.conj().T@L2[2]@u
        phase=np.exp(1j*T*(e[:,None]-e[None,:])+.5j*(s-u_time)*(e[:,None]+e[None,:]))
        spectral=m*m*np.sum(mask*phase*aa.conj()*bb)/len(h)
        phaseerr=float(abs(spectral-cov(X1,X2))); assert phaseerr<2e-12
        coverr=float(abs(cov(F1,F2)-cov(X1,X2)))
        covbound=hs(F1-X1)*hs(F2)+hs(X1)*hs(F2-X2)
        assert coverr<=covbound+2e-13
        paths=None
        if N==3:
            lowerinds=np.array([b for b in range(2**N) if b.bit_count()==m-1])
            lower=H[np.ix_(lowerinds,lowerinds)]
            sp=spectral_paths(h,lower,cs,inds,lowerinds,p,m,.35)
            paths=[discrepancy(x,y) for x,y in zip(sp,cascade(K,a,R,.35))]
            assert max(paths)<2e-12
        rows.append(dict(N=N,p=p,m=m,seed=1106701,hamiltonian_scale=float(scale),
                         coupling_sha256=hashlib.sha256(J.tobytes()).hexdigest(),times=data,
                         covariance=dict(s=s,u=u_time,T=T,error=coverr,bound=covbound,spectral_phase_error=phaseerr),
                         spectral_paths_errors=paths))
    # Exact degeneracy: pinching must remove the whole degenerate energy block.
    energy=np.array([0.,0.,1.]); xx=np.ones((3,3)); mask=energy[:,None]!=energy[None,:]
    assert np.count_nonzero(xx*mask)==4
    return dict(CAR=rows,degenerate_pinching_retained_entries=4)

def scalar_checks():
    rows=[]; K=np.array([[0.,1.],[1.,0.]]); a=np.zeros((1,1)); R=np.array([[0.],[1.]])
    for s in [-2.,-.7,0.,.7,2.]:
        js=[z[0,0] for z in cascade(K,a,R,s)]
        exact=[np.cos(s),s*np.sin(s)/2,(s*np.sin(s)-s*s*np.cos(s))/8]
        err=max(abs(np.array(js)-exact)); assert err<2e-14
        for j in range(3): assert abs(1-sum(js[:j+1]))<=abs(s)**(2*j+2)/factorial(2*j+2)+1e-14
        rows.append(dict(s=s,cascade_error=float(err),remainders=[float(abs(1-sum(js[:j+1]))) for j in range(3)]))
    return rows

if __name__=='__main__':
    OUT.mkdir(parents=True,exist_ok=True)
    dest=OUT/'finite-checks.json'
    if dest.exists(): raise FileExistsError(dest)
    result=dict(force_square_functions=force_counts(),scalar=scalar_checks(),finite=finite_checks(),
                environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
                precision='complex128; exact integer covariance entries; analytical asymptotic bounds, no DS fit')
    dest.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='passed',output=str(dest),CAR_systems=len(result['finite']['CAR']),
                         rational_bound_checks=result['force_square_functions']['rational_bound_checks'])))
