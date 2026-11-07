"""All-sector variance and regression denominators, with exact small CAR checks."""
from pathlib import Path
from fractions import Fraction as F
from math import comb,factorial,prod
import hashlib,json,platform
import numpy as np
from singlet_wick_moments import gue_basis
from cubic_linear_dynamics import weights

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-charge-sum'
def choose(n,k):return comb(n,k) if 0<=k<=n else 0

def small(N,p):
    cs,vs,rs=gue_basis(N,p);charges=np.array([i.bit_count() for i in range(1<<N)])
    ds=np.array([[v@c-c@v for c in cs] for v in vs]);rows=[]
    for m in range(N+1):
        ids=np.flatnonzero(charges==m);d=len(ids);C=choose(m,p)*choose(N-m+p,p)
        bh=F(0);bo=F(0);oo=F(0);amean=F(0)
        for a,va in enumerate(vs):
            aa=sum(z.conj().T@z for z in ds[a])[np.ix_(ids,ids)]
            amean+=F(int(rs[a])*int(round(np.trace(aa).real)),2*d)
            for b,vb in enumerate(vs):
                B=(va@vb+vb@va)[np.ix_(ids,ids)]
                O=sum(u.conj().T@v+v.conj().T@u for u,v in zip(ds[a],ds[b]))[np.ix_(ids,ids)]
                vals=[np.vdot(B,B),np.vdot(B,O),np.vdot(O,O)]
                for z in vals:assert z.imag==0 and z.real==round(z.real) and abs(z)<2**40
                factor=F(int(rs[a])*int(rs[b]),8*d)
                bh+=factor*int(vals[0].real);bo+=factor*int(vals[1].real);oo+=factor*int(vals[2].real)
        assert amean==F(p*(N+1-p)*C,N-m+1)
        if not C:
            assert bh==bo==oo==amean==0;c=F(0);chi=F(0)
        else:
            w=weights(N,m,p);assert sum(w)==1
            chi=sum(w[a]*F(choose(m-a,p)*choose(N-m+p-a,p),C) for a in range(p+1))
            assert 0<=chi<=1 and bh==C*C*(1+chi)
            c=bo/bh;assert oo-c*c*bh>=0
        rows.append(dict(N=N,p=p,m=m,C=C,B_variance=str(bh),chi=str(chi),
                         mean_A=str(amean),O_norm_squared=str(oo),BO=str(bo),c_m=str(c),
                         regression_residual_squared=str(oo-c*c*bh)))
    return rows

def main():
    OUT.mkdir(exist_ok=True);target=OUT/'global-sector-checks.json'
    if target.exists():raise FileExistsError(target)
    checks=0
    for N in range(2,41):
        for p in range(1,min(N,8)+1):
            envelope=prod(F((N+p-2*j)**2,N*N) for j in range(p))
            for m in range(N+1):
                C=choose(m,p)*choose(N-m+p,p)
                normalized=F(4**p*factorial(p)**2*C,N**(2*p))
                assert normalized<=envelope
                if m>=p:
                    cm=choose(m-1,p)*choose(N-m+1+p,p)
                    assert F(cm,C)==F((m-p)*(N-m+p+1),m*(N-m+1))
                    assert F(cm,C)<=p+1
                    assert C==choose(N+p-m,p)*choose(m,p)
                checks+=1
    rows=[]
    for N,p in [(3,1),(4,2),(5,2)]:
        rows+=small(N,p);print('all charge CAR',N,p,'passed',flush=True)
    names=['code/memory_charge_global_checks.py','code/singlet_wick_moments.py',
           'code/cubic_linear_dynamics.py','code/complex_syk_dictionary.py']
    data=dict(status='passed',exact_variance_ratio_cases=checks,exact_CAR_sectors=rows,
              covariance='w=1; real Gaussian-coordinate variances r_a/2; exact integer contractions restored as Fraction.',
              global_bound='C_m <= (Js^2/lambda_N) product_j(1+(p-2j)/N)^2 <= (Js^2/lambda_N) exp(2p/N)',
              scope='All sectors of three small CAR systems; finite variance envelope and adjacent-sector ratios. Analytic all-N bounds and global centered fields are defined in the note.',
              environment=dict(python=platform.python_version(),numpy=np.__version__),
              input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print('passed',checks,'variance cases and',len(rows),'CAR sectors')

if __name__=='__main__':main()
