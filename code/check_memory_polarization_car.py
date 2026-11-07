"""Independent complete CAR check of the weighted polarization identity."""
from pathlib import Path
import hashlib,json,platform,time
import numpy as np
from fractions import Fraction as F
from memory_quadratic_contraction import Contraction
from memory_rank_blocks import casimir_matrix,projectors
from memory_central_polarization import neutral_pair
from cubic_linear_dynamics import weights

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-central-gram'

def main():
    target=OUT/'polarization-car.json'
    if target.exists():raise FileExistsError(target)
    start=time.time();rows=[]
    for N,p,m in [(5,1,3),(6,1,3),(6,1,4)]:
        o=Contraction(N,p,m);count=len(o.rs);n=N-m+1
        pp=projectors(casimir_matrix(N,m),N,min(m,N-m))
        pm=projectors(casimir_matrix(N,m-1),N,min(m-1,N-m+1))
        assert min(len(pp),len(pm))>2*p
        hp=[(o.hs.reshape(count,-1)@num.T/den).reshape(o.hs.shape) for num,den in pp[:p+1]]
        kp=[(o.ks.reshape(count,-1)@num.T/den).reshape(o.ks.shape) for num,den in pm[:p+1]]
        rem=sum(np.kron(c.conj().T,c.T) for c in o.cs)
        gw=o.rs[:,None,None]*o.rs[None,:,None]/2
        def ip(a,b,dim):return float(np.sum(gw*(a.conj()*b).real)/dim)
        def product(arr,a,b):
            val=np.einsum('aij,bjk->abik',arr[a],arr[b])
            if a!=b:val+=np.einsum('aij,bjk->abik',arr[b],arr[a])
            return ((val+val.swapaxes(0,1))/2).reshape(count,count,-1)
        totals=np.zeros(6);max_error=0.;blocks=0
        bs=np.zeros((count,count,o.d*o.d),complex);os=bs.copy();ws=bs.copy();ss=bs.copy()
        for a in range(p+1):
            for b in range(a,p+1):
                upper=product(hp,a,b);lower=product(kp,a,b)
                for r in range(a+b+1):
                    un,ud=pp[r];ln,ld=pm[r]
                    B=upper@un.T/ud/o.C;Bminus=lower@ln.T/ld/o.Cm
                    ell=np.sqrt(float(F(m*(m-r)*(n-r),n)))
                    V=Bminus@rem.T/ell;eps=V-B;t=o.Cm/o.C*ell;hh=m-2*p
                    # Actual operator coefficient from removal, not polarization RHS.
                    O=(lower@ln.T/ld@rem.T+(2*p-m)*upper@un.T/ud)/o.C
                    S=a+b-r;W=( (a+b)*(N+1-a-b)-r*(N+1-r) )/N
                    lhs=W*ip(B,O,o.d)
                    terms=[.5*W*t*(ip(Bminus,Bminus,o.dm)-ip(B,B,o.d)),
                           W*(t-hh)*ip(B,B,o.d),-.5*W*t*ip(eps,eps,o.d)]
                    error=abs(lhs-sum(terms));max_error=max(max_error,error)
                    totals+=np.array([lhs,*terms,W*ip(B,B,o.d),W*ip(Bminus,Bminus,o.dm)])
                    bs+=B;os+=O;ws+=W*B;ss+=S*B;blocks+=1
        fullh=product([o.hs],0,0);fullk=product([o.ks],0,0)
        expectedO=(fullk@rem.T+(2*p-m)*fullh)/o.C
        assert np.max(abs(bs-fullh/o.C))<1e-12
        assert np.max(abs(os-expectedO))<1e-12
        assert abs(ip(ws,os,o.d)-totals[0])<1e-12
        for sector,index in [(m,4),(m-1,5)]:
            omega=weights(N,sector,p)
            expected=sum(omega[a]*omega[b]*neutral_pair(N,sector,a,b) for a in range(p+1) for b in range(p+1))
            assert abs(totals[index]-float(expected))<1e-12
        assert max_error<1e-12
        row=dict(N=N,p=p,m=m,rank_blocks=blocks,ordered_gaussian_pairs=count*count,
                 weighted_polarization_lhs=totals[0],rhs_terms=totals[1:4].tolist(),
                 neutral_moments=totals[4:6].tolist(),max_block_identity_error=max_error,
                 actual_S_O_cross=ip(ss,os,o.d),actual_S_norm_squared=ip(ss,ss,o.d),
                 actual_O_norm_squared=ip(os,os,o.d))
        rows.append(row);print('polarization CAR',N,p,m,'passed',flush=True)
    names=['code/check_memory_polarization_car.py','code/memory_central_polarization.py',
           'code/memory_quadratic_contraction.py','code/memory_rank_blocks.py',
           'code/cubic_linear_dynamics.py','code/complex_syk_dictionary.py','code/singlet_wick_moments.py']
    target.write_text(json.dumps(dict(status='passed',rows=rows,
        scope='All Gaussian coordinate pairs, all input-pair and common output ranks. Weighted polarization and neutral moment are checked against actual CAR coefficients. Square roots and projected coefficient evaluations use float64, tolerance 1e-12. No candidate charged Casimir identity is used.',
        environment=dict(python=platform.python_version(),numpy=np.__version__),elapsed_seconds=time.time()-start,
        input_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names}),indent=2)+'\n')

if __name__=='__main__':main()
