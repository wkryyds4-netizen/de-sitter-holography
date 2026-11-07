"""Exact CAR identity linking the output Casimir deficit to both force orientations."""
from pathlib import Path
from math import comb
import hashlib,json,platform
import numpy as np
import sympy as s
from complex_syk_dictionary import fermions,creators
from memory_rank_blocks import casimir_matrix,projectors,exact_integer

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-oriented-kernel'


def main():
    target=OUT/'force-bridge.json'
    if target.exists():raise FileExistsError(target)
    rows=[];rng=np.random.default_rng(110687)
    for N,p in [(4,1),(4,2),(6,1),(6,2),(6,3)]:
        m=N//2;d1=m-p+1;cs=fermions(N);bs=creators(cs,p)
        Q=np.diag([i.bit_count() for i in range(1<<N)]);I=np.eye(1<<N)
        D=lambda X:Q@X-sum(c.conj().T@X@c for c in cs)
        Db=lambda X:(N*I-Q)@X-sum(c@X@c.conj().T for c in cs)
        ids=[i for i in range(1<<N) if i.bit_count()==m];dim=len(ids)
        C2=casimir_matrix(N,m);Ps=projectors(C2,N,m)
        den=int(s.ilcm(*[dd for _,dd in Ps]))
        knum=sum(a*num*(den//dd) for a,(num,dd) in enumerate(Ps))
        dn=p*den*np.eye(dim*dim,dtype=np.int64)-knum
        for case in range(3):
            J=rng.integers(-2,3,(len(bs),len(bs)))+1j*rng.integers(-2,3,(len(bs),len(bs)))
            J=J+J.conj().T
            H=sum(J[a,b]*aa@bb.conj().T for a,aa in enumerate(bs) for b,bb in enumerate(bs))
            h=H[np.ix_(ids,ids)];sq=h@h
            hdef=(dn@h.ravel()).reshape(dim,dim);hdef2=(dn@hdef.ravel()).reshape(dim,dim)
            hk=(knum@h.ravel()).reshape(dim,dim)
            ch=(C2@h.ravel()).reshape(dim,dim);csq=(C2@sq.ravel()).reshape(dim,dim)
            dg=[H@c-c@H for c in cs]
            Af=sum(x.conj().T@x for x in dg);Acf=sum(x@x.conj().T for x in dg)
            barA=2*p*(H@H)-Db(H@H);double=2*p*barA-D(barA)
            A=Af[np.ix_(ids,ids)];Ac=Acf[np.ix_(ids,ids)];E2=double[np.ix_(ids,ids)]
            anti=lambda x,y:x@y+y@x
            assert np.array_equal(Af,2*p*H@H-D(H@H))
            assert np.array_equal(barA,Acf+anti(H,p*H-Db(H)))
            assert np.array_equal(2*d1*Db(H)[np.ix_(ids,ids)],ch-2*p*m*h)
            assert np.array_equal(d1*den**2*(p*H-Db(H))[np.ix_(ids,ids)],
                                  2*d1*den*hdef+hdef2-den*hdef)
            U2num=anti(h,hdef2)+2*hdef@hdef
            left=d1*(den**2*(anti(h,ch)-csq)-4*hk@hk)
            right=(den**2*d1*N*(A+Ac)+2*den**2*d1*((1-2*p)*A-2*p*Ac+E2)
                   -2*den*d1*anti(h,hdef)-2*d1*U2num
                   +2*(m-2*p)*anti(h,hdef2-den*hdef))
            exact_integer(left);exact_integer(right);assert np.array_equal(left,right),(N,p,case)
            # Direct operator superoperator identity, independent of rank projectors.
            rhs=((m+1)*D(H@H)+m*Db(H@H)-D(Db(H@H)))[np.ix_(ids,ids)]
            assert np.array_equal(2*rhs,csq)
            rows.append(dict(N=N,p=p,m=m,case=case,rank_denominator=den,identities=6,
                             max_cleared_entry=float(np.max(abs(left)))))
        print('bridge',N,p,'passed',flush=True)
    N,p,u=s.symbols('N p u');m=N/2
    xp=p*(N+1-p);xa=(p-u)*(N+1-p+u)
    assert s.factor((xp-xa)/(m-p+1)-2*u-u*(u-1)/(m-p+1))==0
    names=['code/memory_force_orientation_bridge.py','code/memory_rank_blocks.py','code/complex_syk_dictionary.py']
    result=dict(status='passed',records=rows,finite_identities=sum(r['identities'] for r in rows),
                exact_error_identity='N*(S_Casimir-(O+O_creation)/2)=(1-2p)*O-2p*O_creation-W2{h,Dinput}-W2 Uinput2+W2 E_double+(m-2p)*W2{h,Einput}, with Einput=sum u(u-1)h_a/(m-p+1).',
                analytic_scope='Finite raw matrix identities precede W2 projection. L2 control of the double mark and both orientations follows from absolute marked-overlap moments; this code alone does not establish the asymptotic bound.',
                environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__),
                input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
