"""Full finite Gaussian coefficient Gram versus rank-product parity formula."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,platform,time
import numpy as np
import sympy as s
from memory_quadratic_contraction import Contraction,anti
from memory_rank_blocks import casimir_matrix,projectors,exact_integer
from cubic_linear_dynamics import weights
from memory_parity_moments import product_table

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-parity-gram'

def dot(a,b):
    assert max(abs(a).max(),abs(b).max())**2*a.size*2<2**50
    value=np.vdot(a,b).real
    assert value==round(value)
    return int(value)

def main():
    target=OUT/'car-gram.json'
    if target.exists():raise FileExistsError(target)
    start=time.time();rows=[]
    for N,p,m in [(4,1,2),(4,2,2),(6,1,3)]:
        o=Contraction(N,p,m);Ps=projectors(casimir_matrix(N,m),N,m)
        denominator=int(s.ilcm(*[den for _,den in Ps]))
        ranknum=sum(r*num*(denominator//den) for r,(num,den) in enumerate(Ps))
        defnum=p*denominator*np.eye(o.d*o.d,dtype=np.int64)-ranknum
        def values(h):
            square=h@h;dn=(defnum@h.ravel()).reshape(o.d,o.d)
            mixed=anti(h,dn)
            output=2*p*denominator*square-mixed-(ranknum@square.ravel()).reshape(o.d,o.d)
            return exact_integer(np.array([denominator*square,mixed,output]))
        diagonal=[values(h) for h in o.hs];G=[[F(0) for _ in range(3)] for _ in range(3)]
        count=len(o.rs)
        for i in range(count):
            for j in range(i,count):
                if i==j:q=2*diagonal[i]
                else:q=values(o.hs[i]+o.hs[j])-diagonal[i]-diagonal[j]
                factor=F(int(o.rs[i]*o.rs[j])*(1 if i==j else 2),2*(2*denominator)**2*o.d)
                for aa in range(3):
                    for bb in range(aa,3):G[aa][bb]+=factor*dot(q[aa],q[bb])
        for aa in range(3):
            for bb in range(aa):G[aa][bb]=G[bb][aa]
        omega=weights(N,m,p);mean=sum(F(p-a)*omega[a] for a in range(p+1))
        T=s.Matrix([[1,0,0],[-2*s.Rational(mean.numerator,mean.denominator),1,0],[0,0,1]])
        exactG=T*s.Matrix(G)*T.T
        _,_,pi=product_table(N);predict=s.zeros(3)
        for a in range(p+1):
            for b in range(p+1):
                for r,prob in enumerate(pi[a,b]):
                    deficit=a+b-r
                    parity=-1 if deficit%2 else 1
                    f=s.Matrix([1,s.Rational(2*p-a-b)-2*mean,deficit])
                    predict+=o.C**2*omega[a]*omega[b]*(1+parity)*prob*f*f.T
        assert exactG==predict,(N,p,m,exactG,predict)
        assert exactG.det()>=0
        rows.append(dict(N=N,p=p,m=m,covariance_w=1,coordinate_pairs=count*(count+1)//2,
            input_mean_deficit=str(mean),basis=['W2 H^2','W2 {H,D-mean H}','S_N'],
            exact_gram=[[str(x) for x in row] for row in exactG.tolist()]))
        print(N,p,m,'passed',rows[-1]['coordinate_pairs'],'coefficient pairs',flush=True)
    names=['code/check_memory_parity_gram.py','code/memory_parity_moments.py',
        'code/memory_quadratic_contraction.py','code/memory_rank_blocks.py','code/cubic_linear_dynamics.py']
    result=dict(status='passed',rows=rows,elapsed_seconds=time.time()-start,
        method='All symmetric Gaussian coordinate pairs, exact integer-scaled CAR coefficients and Fraction Gram; compared with full rank probability formula.',
        environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
