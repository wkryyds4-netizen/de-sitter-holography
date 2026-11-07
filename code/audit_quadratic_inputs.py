"""Check the older F direction by direct cubic CAR polarization, not its map."""
import hashlib,json,platform
import numpy as np
from memory_quadratic_contraction import Contraction,ROOT,OUT,adj


def main():
    target=OUT/'linear-input-audit.json'
    if target.exists():raise FileExistsError(target)
    rows=[]
    for N,p,m in [(3,1,2),(4,2,3),(5,2,3)]:
        o=Contraction(N,p,m)
        def T(h,k):
            value=np.zeros_like(h)
            for c in o.cs:
                f=c
                for _ in range(3):f=k@f-f@h
                value+=adj(c)@f
            return (value+adj(value))/2
        count=len(o.rs);x=np.array([(7*i+3)%11-5 for i in range(count)])
        cases=list(zip(o.hs,o.ks))+[(np.einsum('a,aij->ij',x,o.hs),np.einsum('a,aij->ij',x,o.ks))]
        for h,k in cases:
            base=T(h,k)
            predicted=sum(r*(T(h+u,k+v)+T(h-u,k-v)-2*base)/4 for r,u,v in zip(o.rs,o.hs,o.ks))
            assert np.array_equal(predicted,o.tlinear(h,k))
        rows.append(dict(N=N,p=p,m=m,directions=len(cases)))
        print('Direct cubic audit',N,p,m,len(cases),flush=True)
    names=['code/audit_quadratic_inputs.py','code/memory_quadratic_contraction.py']
    target.write_text(json.dumps(dict(status='passed',rows=rows,
        scope='Direct CAR cubic T3 and Gaussian polarization verify the F input map on every real coordinate and one integer sample in each system.',
        environment=dict(python=platform.python_version(),numpy=np.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')


if __name__=='__main__':main()
