"""General-charge CAR bridge to the two actual force orientations."""
from pathlib import Path
import hashlib,json,platform
import numpy as np
import sympy as s
from complex_syk_dictionary import fermions,creators
from memory_rank_blocks import casimir_matrix,projectors,exact_integer

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-central-spectrum'

def main():
    target=OUT/'force-bridge.json'
    if target.exists():raise FileExistsError(target)
    rows=[];rng=np.random.default_rng(110691)
    for N,p,m in [(5,1,2),(5,2,3),(6,1,2),(6,1,4),(6,2,2),(6,2,4)]:
        n=N-m+1;d1=m-p+1;cs=fermions(N);bs=creators(cs,p)
        Q=np.diag([i.bit_count() for i in range(1<<N)]);eye=np.eye(1<<N)
        D=lambda X:Q@X-sum(c.conj().T@X@c for c in cs)
        Db=lambda X:(N*eye-Q)@X-sum(c@X@c.conj().T for c in cs)
        ids=[i for i in range(1<<N) if i.bit_count()==m];dim=len(ids)
        C2=casimir_matrix(N,m);Ps=projectors(C2,N,min(m,N-m))
        den=int(s.ilcm(*[dd for _,dd in Ps]));knum=sum(a*num*(den//dd) for a,(num,dd) in enumerate(Ps))
        dn=p*den*np.eye(dim*dim,dtype=np.int64)-knum
        anti=lambda a,b:a@b+b@a
        for case in range(3):
            J=rng.integers(-2,3,(len(bs),len(bs)))+1j*rng.integers(-2,3,(len(bs),len(bs)));J=J+J.conj().T
            H=sum(J[a,b]*aa@bb.conj().T for a,aa in enumerate(bs) for b,bb in enumerate(bs))
            h=H[np.ix_(ids,ids)];sq=h@h
            hd=(dn@h.ravel()).reshape(dim,dim);hd2=(dn@hd.ravel()).reshape(dim,dim)
            hk=(knum@h.ravel()).reshape(dim,dim);ch=(C2@h.ravel()).reshape(dim,dim);csq=(C2@sq.ravel()).reshape(dim,dim)
            dg=[H@c-c@H for c in cs];Af=sum(d.conj().T@d for d in dg);Acf=sum(d@d.conj().T for d in dg)
            barA=2*p*(H@H)-Db(H@H);double=2*p*barA-D(barA)
            A=Af[np.ix_(ids,ids)];Ac=Acf[np.ix_(ids,ids)];E2=double[np.ix_(ids,ids)]
            assert np.array_equal(Af,2*p*H@H-D(H@H))
            assert np.array_equal(barA,Acf+anti(H,p*H-Db(H)))
            assert np.array_equal(2*d1*Db(H)[np.ix_(ids,ids)],ch-2*p*(N-m)*h)
            ein=hd2+(N-2*m-1)*den*hd
            assert np.array_equal(d1*den**2*(p*H-Db(H))[np.ix_(ids,ids)],2*d1*den*hd+ein)
            U2=anti(h,hd2)+2*hd@hd
            left=d1*(den**2*(anti(h,ch)-csq)-4*hk@hk)
            right=(2*d1*den**2*((n-2*p)*A+(m-2*p)*Ac+E2)
                   +2*d1*den*(2*m-N-1)*anti(h,hd)-2*d1*U2+2*(m-2*p)*anti(h,ein))
            exact_integer(left);exact_integer(right);assert np.array_equal(left,right),(N,p,m,case)
            rhs=(n*D(H@H)+m*Db(H@H)-D(Db(H@H)))[np.ix_(ids,ids)]
            assert np.array_equal(2*rhs,csq)
            rows.append(dict(N=N,p=p,m=m,case=case,identities=6,rank_denominator=den,max_cleared_entry=float(np.max(abs(left)))))
        print('general force bridge',N,p,m,'passed',flush=True)
    NN,pp,mm,u=s.symbols('N p m u');xp=pp*(NN+1-pp);xa=(pp-u)*(NN+1-pp+u)
    assert s.factor((xp-xa)/(mm-pp+1)-2*u-u*(u-1+NN-2*mm)/(mm-pp+1))==0
    names=['code/memory_central_force_bridge.py','code/memory_rank_blocks.py','code/complex_syk_dictionary.py']
    data=dict(status='passed',records=rows,finite_identities=sum(r['identities'] for r in rows),
        exact_error_identity='N*(S_proxy-(O_minus+O_plus)/2)=(N/2-m+1-2p)O_minus+(m-N/2-2p)O_plus+(2m-N-1)W2{h,D_in h}-W2 U_in2+W2 E_double+(m-2p)W2{h,E_in}',
        E_in='sum_u u*(u-1+N-2m)/(m-p+1) h_(p-u)',
        scope='Exact full-Fock CAR and cleared-rational rank identities at general charge, before W2. The central L2 bound uses the independent trace/overlap argument; no half-filled parity substitution.',
        environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__),
        input_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names})
    target.write_text(json.dumps(data,indent=2)+'\n')

if __name__=='__main__':main()
