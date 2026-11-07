"""Adjacent-charge isometry: symbolic scalar and complete finite-map checks."""
from pathlib import Path
import hashlib,json,platform,time
from fractions import Fraction as F
import numpy as np
import sympy as s
from memory_quadratic_contraction import Contraction
from memory_rank_blocks import casimir_matrix,projectors,exact_integer

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-central-cross'


def main():
    target=OUT/'transport.json'
    if target.exists():raise FileExistsError(target)
    start=time.time()
    e,lam,x,rho=s.symbols('e lam x rho',positive=True)
    p=1/e;N=p*p/lam;m=N/2+x*p;n=N-m+1;r=rho*p
    ratio=(m-p)*(n+p)/(m*n)
    ell2=m*(m-r)*(n-r)/n
    drift=s.factor(s.limit(s.cancel((ratio**2*ell2-(m-r)**2)/(2*m)),e,0))
    expected=2*lam*(x*(2-rho)-1)
    assert s.factor(drift-expected)==0
    ratio_second=s.factor(s.limit((ratio-1)/e**2,e,0))
    assert s.factor(ratio_second-4*lam**2*(2*x-1))==0
    q,u,v=s.symbols('q u v');a=p-u;b=p-v
    xa=a*(N+1-a);xb=b*(N+1-b)
    cf=m*(N-m)*(N+1)/N;g=(m-p)**2/m;eta=p*p/m
    neutral=-2*a*b*(1+q)+xa*xb/cf
    candidate=(eta*neutral+g*(xa*xb/(m*(m-1)*(N+1))+(m-2*a*b)*(-4*lam*q/m)))/N
    cross_limit=s.factor(s.limit(s.cancel(candidate),e,0))
    assert s.factor(cross_limit-(4*lam**2*(1+q)+2*lam*(1-q)))==0
    # Complete operator-map check, not selected vectors. Common rank is >=2p.
    NN,pp,mm=6,1,3;o=Contraction(NN,pp,mm);nn=NN-mm+1
    lower=projectors(casimir_matrix(NN,mm-1),NN,min(mm-1,NN-mm+1))
    upper=projectors(casimir_matrix(NN,mm),NN,min(mm,NN-mm))
    rem=exact_integer(sum(np.kron(c.conj().T,c.T) for c in o.cs)).real.astype(np.int64)
    ranks=[]
    for rr,(num,den) in enumerate(lower[:2*pp+1]):
        ell=F(mm*(mm-rr)*(nn-rr),nn)
        lhs=o.dm*ell.denominator*(num@rem.T@rem@num)
        rhs=o.d*ell.numerator*den*num
        assert np.array_equal(lhs,rhs)
        ranks.append(dict(rank=rr,dimension=int(np.trace(num)//den),ell_squared=str(ell),entries=lhs.size))
    count=len(o.rs)
    hp=np.einsum('aij,bjk->abik',o.hs,o.hs);hp=(hp+hp.swapaxes(0,1))/2
    kp=np.einsum('aij,bjk->abik',o.ks,o.ks);kp=(kp+kp.swapaxes(0,1))/2
    bp=hp.reshape(count,count,-1)/o.C;bm=kp.reshape(count,count,-1)/o.Cm
    normalized_diff=np.zeros_like(bp);identity_error=0.;iso_error=0.
    weights=o.rs[:,None,None]*o.rs[None,:,None]/2
    def norm2(z,dim):return float(np.sum(weights*abs(z)**2)/dim)
    for rr,((ln,ld),(un,ud)) in enumerate(zip(lower,upper)):
        vl=bm@ln.T/ld;vu=bp@un.T/ud
        el=np.sqrt(float(F(mm*(mm-rr)*(nn-rr),nn)))
        transported=vl@rem.T/el
        iso_error=max(iso_error,abs(norm2(transported,o.d)-norm2(vl,o.dm)))
        tr=float(F(o.Cm,o.C))*el
        oblock=(kp.reshape(count,count,-1)@ln.T/ld@rem.T+(2*pp-mm)*hp.reshape(count,count,-1)@un.T/ud)/o.C
        identity_error=max(identity_error,float(abs(tr*transported-(mm-2*pp)*vu-oblock).max()))
        normalized_diff+=transported-vu
    assert identity_error<1e-12 and iso_error<1e-12
    names=['code/check_memory_central_transport.py','code/memory_quadratic_contraction.py','code/memory_rank_blocks.py','code/complex_syk_dictionary.py','code/singlet_wick_moments.py']
    data=dict(status='passed',scalar_drift=str(drift),variance_ratio_second=str(ratio_second),candidate_cross_limit=str(cross_limit),
        finite=dict(N=NN,p=pp,m=mm,rank_map_checks=ranks,ordered_gaussian_pairs=count*count,
            max_isometry_norm_error=iso_error,max_transport_identity_error=identity_error,
            adjacent_product_norm_squared=norm2(normalized_diff,o.d)),
        scope='Exact symbolic scalar asymptotics and complete rational isometry on common output ranks 0..2p. Floating full Gaussian-coefficient transport identity uses square roots. Uniform telescoping/tail proof is in the note. The candidate cross limit does not prove its general-charge recoupling lemma.',
        elapsed_seconds=time.time()-start,environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__),
        input_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names})
    target.write_text(json.dumps(data,indent=2)+'\n');print('transport passed',data['finite'],flush=True)

if __name__=='__main__':main()
