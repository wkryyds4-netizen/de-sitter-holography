"""Independent polarization checks and removal of sample-dependent scalar parts."""
from pathlib import Path
from fractions import Fraction as F
from math import lcm
import hashlib,json,platform,time
import numpy as np
import sympy as s
from memory_quadratic_contraction import Contraction,inner_exact,adj,anti,ROOT,OUT


def trfrac(x):
    a=8*np.trace(x)
    assert a.imag==0 and a.real==round(a.real)
    return F(int(a.real),8*x.shape[0])


def main():
    target=OUT/'scalar-removed-checks.json'
    if target.exists():raise FileExistsError(target)
    parent=json.loads((OUT/'finite-contractions.json').read_text())
    rows=[]
    for previous in parent['rows']:
        N,p,m=[previous[k] for k in ['N','p','m']];o=Contraction(N,p,m);count=len(o.rs)
        tbar=F(previous['tbar']);den=tbar.denominator;num=tbar.numerator
        def vals(h,k):
            f=den*o.tlinear(h,k)-num*h
            return [h@h,o.A(h,k),anti(h,f)/2,o.contract(h,k)]
        diag=[vals(h,k) for h,k in zip(o.hs,o.ks)]
        G=[[F(0) for _ in range(4)] for _ in range(4)]
        scal=[[F(0) for _ in range(4)] for _ in range(4)]
        verified=0
        # Entire quadratic map checked independently for the two smaller systems.
        # Larger system checks a deterministic 100-point subset of mixed directions.
        for a in range(count):
            for b in range(a,count):
                if a==b:
                    bs=diag[a];h=o.hs[a];k=o.ks[a]
                else:
                    h=o.hs[a]+o.hs[b];k=o.ks[a]+o.ks[b]
                    mix=vals(h,k);bs=[(v-x-y)/2 for v,x,y in zip(mix,diag[a],diag[b])]
                if N<5 or (a*count+b)%53==0 or a==b==0:
                    assert np.array_equal(o.contract(h,k),o.direct_contract(h,k)),(N,p,m,a,b)
                    verified+=1
                factor=F(int(o.rs[a]*o.rs[b]),2)*(1 if a==b else 2)
                means=[trfrac(v) for v in bs]
                for i in range(4):
                    for j in range(4):
                        sc=means[i]*means[j]
                        scal[i][j]+=factor*sc
                        G[i][j]+=factor*(inner_exact(bs[i],bs[j])-sc)
        for i in range(4):
            for j in range(4):
                scale=den**(int(i==2)+int(j==2));G[i][j]/=scale;scal[i][j]/=scale
                assert G[i][j]+scal[i][j]==F(previous['exact_gram'][i][j])
        mat=s.Matrix([[s.Rational(x.numerator,x.denominator) for x in row] for row in G])
        basis=[0,1,2];small=mat.extract(basis,basis);g=mat.extract(basis,[3])
        # Remove redundant old directions if a small system degenerates.
        pivots=list(small.rref()[1]);ids=[basis[i] for i in pivots]
        coeff=mat.extract(ids,ids).inv()*mat.extract(ids,[3])
        residual=s.factor(mat[3,3]-(mat.extract(ids,[3]).T*coeff)[0])
        assert residual>=0
        x=np.array(previous['deterministic_couplings']);h=np.einsum('a,aij->ij',x,o.hs);k=np.einsum('a,aij->ij',x,o.ks)
        values=vals(h,k);values[2]=values[2]/den
        exact=lambda a:s.Matrix(a).applyfunc(lambda z:s.nsimplify(z,rational=True))
        hq=exact(h);vv=[exact(a) for a in values]
        z=vv[3]-sum((c*vv[i] for i,c in zip(ids,coeff)),s.zeros(o.d))
        z=z-s.trace(z)/o.d*s.eye(o.d)
        comm=(hq*z-z*hq).applyfunc(s.expand)
        norm=s.expand(s.trace(comm.conjugate().T*comm)/o.d)
        assert norm>=0
        row=dict(N=N,p=p,m=m,independent_polarization_directions=verified,
                 traceless_gram=[[str(z) for z in row] for row in G],
                 random_scalar_gram=[[str(z) for z in row] for row in scal],
                 independent_old_basis=ids,projection_coefficients=list(map(str,coeff)),
                 residual_norm_squared=str(residual),residual_commutator_norm_squared=str(norm))
        rows.append(row);print(json.dumps(row),flush=True)
    names=['code/check_quadratic_contraction.py','code/memory_quadratic_contraction.py','results/susskind-memory-quadratic-contraction/finite-contractions.json']
    result=dict(status='passed',rows=rows,
        scope='Exact finite-sector covariance, sample-wise trace removed. No large-N or spectral claim.',
        environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
