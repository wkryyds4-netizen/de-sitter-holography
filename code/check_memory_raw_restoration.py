"""Independent CAR/renewal check of the restored raw identity and counterterms."""
from pathlib import Path
import json,hashlib,platform
import numpy as np
import sympy as s
from complex_syk_dictionary import fermions,creators,comm
from bilocal_channel import blocks
from memory_force_bridge import exact
from memory_raw_spectrum import coefficients,q,lam

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-raw-spectrum'


def finite(N,p,m):
    cs=fermions(N);bs=creators(cs,p);rng=np.random.default_rng(110689+N*100+p*10+m)
    z=rng.integers(-1,2,(len(bs),len(bs)))+1j*rng.integers(-1,2,(len(bs),len(bs)))
    z=z+z.conj().T;H=sum(z[a,b]*aa@bb.conj().T for a,aa in enumerate(bs) for b,bb in enumerate(bs))
    ix,h,V,K=blocks(N,p,H,cs,m);h,V,K=map(exact,[h,V,K]);rr=s.Rational(p,m)
    ex=lambda z:z.applyfunc(s.expand)
    J=lambda z:ex((h*z+z*h)/2)
    L=lambda z:ex(h*z-z*h)
    forms=[[c] for c in cs]
    for _ in range(4):
        for row in forms:row.append(comm(H,row[-1]))
    Bf=lambda a,b:exact(sum(row[a].conj().T@row[b] for row in forms)[np.ix_(ix,ix)])
    A=Bf(1,1);T3=ex((Bf(0,3)+Bf(3,0))/2);core=ex(Bf(2,2)+L(L(A))/4)
    charge=np.diag([i.bit_count() for i in range(1<<N)])
    D=lambda X:charge@X-sum(c.conj().T@X@c for c in cs)
    tcal3=exact((3*p*(H@H@H)-D(H@H@H))[np.ix_(ix,ix)])
    tcal4=exact((4*p*(H@H@H@H)-D(H@H@H@H))[np.ix_(ix,ix)])
    E=ex((A*h*h+h*h*A)/2);M=ex(h*A*h);S=ex(J(tcal3)-E-M)
    CC=ex(tcal4-2*E-M-2*S)
    co={k:s.factor(v.subs({q:s.Rational(1,2),lam:s.Rational(p*p,N)})) for k,v in coefficients().items()}
    raw=ex(CC-co['alpha']*E-co['beta']*M-co['gamma']*S-co['delta']*h**4)
    B={j:ex(V.conjugate().T*K**j*V/m) for j in range(1,5)};bn={}
    for n in range(1,5):bn[n]=ex(B[n]-sum((bn[j]*B[n-j] for j in range(1,n)),s.zeros(h.rows)))
    Z2=ex(bn[4]-2*J(bn[3])+J(J(bn[2])))
    # Arbitrary commuting polynomials test the algebra without asymptotic assumptions.
    KA=5*s.eye(h.rows)+s.Rational(2,7)*h*h;P=3*h+h**3/11
    Y=ex(A-KA);U=ex(T3-P)
    D2=ex(2*rr*J(P)+3*rr*rr*J(J(KA))-KA**2/m-m*rr**4*h**4)
    eps=ex((rr/2+rr**2/4)*L(L(Y))+2*rr*J(U)+3*rr**2*J(J(Y))-(KA*Y+Y*KA+Y*Y)/m)
    counter=ex(D2+co['gJ']*h*h*KA+co['bR']*h*P+co['delta']*h**4)
    restored=ex(raw+co['gJ']*J(J(Y))+co['bD']*L(L(Y))+co['bR']*J(U)+eps)
    checks=[S-J(T3)-J(J(A)),core-raw-co['gJ']*J(J(A))-co['bD']*L(L(A))-co['bR']*J(T3)-co['delta']*h**4,
            m*Z2-counter-restored,h*counter-counter*h]
    assert all(all(s.expand(v)==0 for v in z) for z in checks)
    norm=lambda z:s.expand(s.trace(z.conjugate().T*z)/h.rows)
    norms=[norm(z) for z in [A,core,m*Z2,restored]]
    assert norms[0]>0 and norms[1]>0 and norms[3]>0
    if m>p:assert norms[2]>0
    print('raw finite',N,p,m,'passed',flush=True)
    return dict(N=N,p=p,m=m,seed=110689+N*100+p*10+m,exact_identities=len(checks),
                norms_A_core_mZ2_restored=list(map(str,norms)),q_for_algebra='1/2',lambda_for_algebra=str(s.Rational(p*p,N)),
                scope='Exact identity for arbitrary rational q and commuting KA/P, not a finite-N spectral equality.')


def main():
    target=OUT/'finite-restoration.json'
    if target.exists():raise FileExistsError(target)
    rows=[finite(*case) for case in [(4,1,2),(4,2,2),(5,2,3),(6,2,3)]]
    names=['code/check_memory_raw_restoration.py','code/memory_raw_spectrum.py','code/complex_syk_dictionary.py',
           'code/bilocal_channel.py','code/memory_force_bridge.py']
    target.write_text(json.dumps(dict(status='passed',records=rows,
        environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')


if __name__=='__main__':main()
