"""Nonzero integer systems and direct finite CAR words for fourth marks."""
from pathlib import Path
from itertools import combinations,permutations
from fractions import Fraction as F
import json,hashlib,platform
import numpy as np
import sympy as s
from complex_syk_dictionary import fermions,creators,comm
from bilocal_channel import blocks
from memory_force_bridge import exact
from check_multimark_car import local_coeff
from memory_quartic import weights

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-quartic'


def finite():
    rows=[]
    for N,p,m in [(3,1,2),(4,2,2),(5,2,3)]:
        cs=fermions(N);bs=creators(cs,p);rng=np.random.default_rng(110670+N)
        j=rng.integers(-2,3,(len(bs),len(bs)))+1j*rng.integers(-2,3,(len(bs),len(bs)));j=j+j.conj().T
        H=sum(j[a,b]*aa@bb.conj().T for a,aa in enumerate(bs) for b,bb in enumerate(bs))
        ix,h,C,K=blocks(N,p,H,cs,m);h,C,K=map(exact,[h,C,K]);r=s.Rational(p,m)
        f=[[c] for c in cs]
        for _ in range(4):
            for row in f:row.append(comm(H,row[-1]))
        BB=lambda a,b:exact(sum(row[a].conj().T@row[b] for row in f)[np.ix_(ix,ix)])
        A=BB(1,1);T3=(BB(0,3)+BB(3,0))/2;T4=(BB(0,4)+BB(4,0))/2;B22=BB(2,2)
        LH=lambda x:h*x-x*h
        JH=lambda x:(h*x+x*h)/2
        B={0:s.eye(h.rows),**{j:(C.conjugate().T*K**j*C/m).applyfunc(s.expand) for j in range(1,5)}};beta={}
        for n in range(1,5):beta[n]=(B[n]-sum((beta[j]*B[n-j] for j in range(1,n)),s.zeros(h.rows))).applyfunc(s.expand)
        G=beta[2];Z2=(beta[4]-2*JH(beta[3])+JH(JH(G))).applyfunc(s.expand)
        direct=B22+(1+r)**2*LH(LH(A))/4+2*r*JH(T3)+3*r*r*JH(JH(A))-A*A/m-m*r**4*h**4
        assert all(s.expand(v)==0 for v in m*Z2-direct)
        assert all(s.expand(v)==0 for v in T4-B22-2*LH(LH(A)))
        assert all(s.expand(v)==0 for v in BB(0,4)-T4-2*LH(T3))
        # Exact separation of the return/drift terms for arbitrary conserved polynomials.
        c=s.Rational(2,7);a0=5;KA=a0*s.eye(h.rows)+c*h*h;PT=3*h+h**3/11
        Y=A-KA;U=T3-PT
        counter=2*r*JH(PT)+3*r*r*JH(JH(KA))-KA**2/m-m*r**4*h**4
        remainder=(r/2+r*r/4)*LH(LH(Y))+2*r*JH(U)+3*r*r*JH(JH(Y))-(KA*Y+Y*KA+Y*Y)/m
        core=B22+LH(LH(A))/4
        assert all(s.expand(v)==0 for v in m*Z2-counter-core-remainder)
        normalized=[x.applyfunc(s.expand) for x in [A,B22,m*Z2]]
        norms=[s.expand(s.trace(x.conjugate().T*x)/h.rows) for x in normalized]
        assert norms[0]>0 and norms[1]>0
        if m>p:assert norms[2]>0
        rows.append(dict(N=N,p=p,m=m,identities=4,norm_squared=list(map(str,norms)),zero_memory_control=m==p))
    return rows


def words():
    rng=np.random.default_rng(110675);records=[];nonzero=negative=0
    for N,p in [(3,1),(4,2)]:
        cs=fermions(N);sets=list(combinations(range(N),p));bs=dict(zip(sets,creators(cs,p)))
        dim=1<<N;Q=np.diag([i.bit_count() for i in range(dim)])
        D=lambda x:Q@x-sum(c.conj().T@x@c for c in cs)
        mark=lambda seq:len(seq)*p*product(seq)-D(product(seq))
        def basis(seq):
            a,b,c,d=seq;k01=mark([a,b])@c@d;k12=a@mark([b,c])@d;k23=a@b@mark([c,d])
            k02=mark([a,b,c])@d-k01-k12;k13=a@mark([b,c,d])-k12-k23
            k03=mark(seq)-k01-k12-k23-k02-k13
            return [(k01+k23)/2,k12,(k02+k13)/2,k03,product(seq)]
        for case in range(32):
            perm=list(permutations(range(4)))[case%24];pa=tuple((i,4+perm[i]) for i in range(4))
            pi={a:b for a,b in pa}|{b:a for a,b in pa}
            chosen=[sets[rng.integers(len(sets))] for _ in range(8)]
            if case<12:chosen=[sets[case%len(sets)]]*8
            ops=[bs[chosen[i]]@bs[chosen[pi[i]]].conj().T for i in range(8)]
            lft,rgt=basis(ops[:4]),basis(ops[4:]);cache={}
            for i,wi in enumerate(weights):
                for j,wj in enumerate(weights):
                    expected=F(0)
                    for ma,ca in (wi or {None:F(1)}).items():
                        for mb,cb in (wj or {None:F(1)}).items():
                            marks=tuple(x for x in [ma,None if mb is None else tuple(v+4 for v in mb)] if x is not None)
                            if marks not in cache:cache[marks]=local_coeff(N,p,chosen,pa,marks)[-1]
                            expected+=ca*cb*cache[marks]
                    value=np.trace(lft[i]@rgt[j])/dim
                    assert value.imag==0 and F(float(value.real))==expected
                    nonzero+=int(expected!=0);negative+=int(expected<0)
            records.append(dict(N=N,p=p,pairing=pa,sets=chosen))
    assert nonzero>=100
    return dict(words=len(records),comparisons=25*len(records),nonzero=nonzero,negative=negative,records=records)


def product(seq):
    out=np.eye(len(seq[0]),dtype=complex)
    for a in seq:out=out@a
    return out


def main():
    path=OUT/'independent-car.json'
    if path.exists():raise FileExistsError(path)
    a=finite();print('Nonzero deterministic systems verified',flush=True);b=words()
    names=['code/check_memory_quartic_car.py','code/memory_quartic.py','code/check_multimark_car.py',
           'code/complex_syk_dictionary.py','code/bilocal_channel.py','code/memory_force_bridge.py']
    path.write_text(json.dumps(dict(status='passed',finite=a,local_words=b,
        environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')
    print('CAR words',b['words'],'nonzero',b['nonzero'],'negative',b['negative'],flush=True)


if __name__=='__main__':main()
