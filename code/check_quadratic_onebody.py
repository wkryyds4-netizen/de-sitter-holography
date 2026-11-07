"""Analytic one-body GUE control, independent of the charged contraction formula."""
import hashlib,json,platform
import numpy as np
from memory_quadratic_contraction import Contraction,ROOT,OUT


def main():
    target=OUT/'onebody-controls.json'
    if target.exists():raise FileExistsError(target)
    rows=[]
    for N in [2,3,4,5]:
        # Build the same one-particle J used to form dGamma(J) in every charge.
        J=np.array([[((3*i+7*j+2)%9)-4+1j*(((5*i+2*j)%7)-3) for j in range(N)] for i in range(N)])
        J=J+J.conj().T
        for m in range(1,N+1):
            o=Contraction(N,1,m)
            # Recover arbitrary one-body h and k from the explicit CAR blocks.
            h=sum(J[i,j]*adj_ci@cj for i,adj_ci in enumerate([c.conj().T for c in o.cs]) for j,cj in enumerate(o.cs))
            # k is dGamma(J) on charge m-1, recovered by an independent constructor.
            from complex_syk_dictionary import fermions
            cs=fermions(N);full=sum(J[i,j]*cs[i].conj().T@cs[j] for i in range(N) for j in range(N))
            ix=[a for a in range(1<<N) if a.bit_count()==m-1]
            k=full[np.ix_(ix,ix)]
            A=o.A(h,k);q=o.contract(h,k)
            expected=3*N*A+2*np.trace(J)*h+m*np.trace(J@J)*np.eye(o.d)
            assert np.array_equal(q,expected)
            assert np.array_equal(h@q,q@h)
            J4=J@J@J@J
            raw=sum(J4[i,j]*ci.conj().T@cj for i,ci in enumerate(o.cs) for j,cj in enumerate(o.cs))
            assert np.array_equal(o.core(h,k),raw)
            mean=sum(r*o.contract(hh,kk)/4 for r,hh,kk in zip(o.rs,o.hs,o.ks))
            assert np.array_equal(mean,m*(2*N*N+1)*np.eye(o.d))
            rows.append(dict(N=N,m=m,core_mean=m*(2*N*N+1),identities=4))
    names=['code/check_quadratic_onebody.py','code/memory_quadratic_contraction.py','code/complex_syk_dictionary.py']
    target.write_text(json.dumps(dict(status='passed',rows=rows,identities=4*len(rows),
        analytic_result='For p=1: core=dGamma(J^4); Q=3N A+2 tr(J) H+m tr(J^2) I; E core=m(2N^2+1) I. Every realization commutes with H.',
        environment=dict(python=platform.python_version(),numpy=np.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')
    print('One-body controls:',len(rows),'systems,',4*len(rows),'exact identities',flush=True)


if __name__=='__main__':main()
