"""Check derived coefficients on coherent rank inputs, preserving complex entries."""
from pathlib import Path
from functools import reduce
from math import gcd
import hashlib,json,platform,time
import numpy as np
import sympy as sp
from singlet_wick_moments import gue_basis
from memory_quadratic_contraction import Contraction,adj,anti
from memory_rank_blocks import casimir_matrix,projectors,exact_integer
from memory_two_charge_filter import filters

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-two-charge'

def words(o,h,k):
    h2=h@h;k2=k@k;r2=o.rem(k2);I=np.eye(o.dm)
    ph=o.phi('hh',h);pk=o.phi('kk',k)
    u=o.gamma(k,I)+o.gamma(I,k)
    v=o.gamma(I,k2)+o.gamma(k,k)+o.gamma(k2,I)
    uh=sum(o.phi('hk',h@adj(c)@k)@c+o.phi('hk',h@adj(c))@k@c for c in o.cs)
    return dict(D4=o.rem(3*o.Cm*k2+anti(k,pk)+o.phi('kk',k2)),
        D13=anti(h,o.rem(2*o.Cm*k+pk))+v+adj(v),
        D22=2*o.C*r2+2*o.m*o.Cm*h2+h@u+uh+adj(u)@h+adj(uh),
        D121=o.phi('hh',r2)+u@h+h@adj(u)+o.m*o.Cm*h2,
        D40=3*o.C*h2+anti(h,ph)+o.phi('hh',h2),Q=o.contract(h,k),
        plus=h2,minus=r2)

def check_combination(actual,plus,minus,g):
    for component in ['real','imag']:
        arrays=[getattr(z,component).ravel() for z in [actual,plus,minus]]
        assert all(np.array_equal(x,np.rint(x)) for x in arrays)
        assert all(np.max(abs(x))<2**50 for x in arrays)
        for x,y,z in zip(*arrays):
            assert int(x)==g[0]*int(y)+g[1]*int(z),(x,y,z,g)

def main():
    OUT.mkdir(exist_ok=True);target=OUT/'filter-checks.json'
    if target.exists():raise FileExistsError(target)
    rows=[];start=time.time();total=0
    for N,p,m in [(4,1,2),(5,2,3),(6,2,3),(6,2,4),(6,3,4)]:
        o=Contraction(N,p,m);pin=projectors(casimir_matrix(N,p),N,min(p,N-p))
        pout=projectors(casimir_matrix(N,m),N,min(m,N-m))
        _,basis,_=gue_basis(N,p);ix=np.array([x for x in range(1<<N) if x.bit_count()==p])
        bp=basis[:,ix[:,None],ix].reshape(len(o.rs),-1)
        hs=[];ks=[]
        for a in range(min(p,m,N-m)+1):
            x=np.array([((7+2*a)*i+3+a)%17-8 for i in range(len(o.rs))])
            J=x@bp;num,den=pin[a];projected=num@J
            exact_integer(projected)
            common=reduce(gcd,[abs(int(v)) for v in np.r_[projected.real,projected.imag]])
            assert common
            projected/=common
            coeff=(bp.conj()@projected)*o.rs/2
            assert not np.any(coeff.imag);exact_integer(coeff)
            hs.append(np.einsum('a,aij->ij',coeff.real,o.hs))
            ks.append(np.einsum('a,aij->ij',coeff.real,o.ks))
        diagonal=[words(o,h,k) for h,k in zip(hs,ks)]
        checks=0;nonzero=0
        for a in range(len(hs)):
            for b in range(a,len(hs)):
                if a==b:ww=diagonal[a]
                else:
                    mixed=words(o,hs[a]+hs[b],ks[a]+ks[b])
                    ww={key:mixed[key]-diagonal[a][key]-diagonal[b][key] for key in mixed}
                for r,(num,den) in enumerate(pout):
                    projected={key:8*(num@x.ravel()) for key,x in ww.items()}
                    f=filters(N,p,m,a,b,r)
                    for key in ['D4','D13','D22','D121','D40','Q']:
                        g=f['g'] if key=='Q' else f['words'][key]
                        check_combination(projected[key],projected['plus'],projected['minus'],g)
                        checks+=1
                    nonzero+=int(np.any(projected['Q']))
        total+=checks;rows.append(dict(N=N,p=p,m=m,identities=checks,nonzero_Q_blocks=nonzero))
        print(rows[-1],flush=True)
    paths=['code/memory_two_charge_filter.py','code/check_memory_two_charge_filter.py',
        'code/memory_quadratic_contraction.py','code/memory_rank_blocks.py','code/singlet_wick_moments.py']
    result=dict(status='passed',rows=rows,identities=total,elapsed_seconds=time.time()-start,
        scope='Exact rank-pure input realization checks of five separate single contractions and their core combination; analytic proof needed for general tensors.',
        environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=sp.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in paths})
    target.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
