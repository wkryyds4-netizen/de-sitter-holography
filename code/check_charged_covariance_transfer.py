"""Exhaustively test the universal S Psi identity on rectangular matrix units."""
from pathlib import Path
import hashlib,json,platform,time
import numpy as np
import sympy as s
from memory_quadratic_contraction import Contraction
from memory_rank_blocks import casimir_matrix,projectors,exact_integer
from memory_two_charge_filter import filters

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-two-charge'

def main():
    target=OUT/'transfer-checks.json'
    if target.exists():raise FileExistsError(target)
    start=time.time();rows=[]
    m,n,r=s.symbols('m n r');D=m+n-2*r
    B=(m-r)/(n-r)
    norm=s.factor(m+n*B-2*(m-r))
    identities=[norm-r*D/(n-r),r/norm-(n-r)/D,
        ((m-r)-n*B)/norm+(m-r)/D]
    assert all(s.factor(x)==0 for x in identities)
    for N,p,charge in [(4,1,2),(5,2,3),(6,2,3),(6,2,4),(6,3,4)]:
        o=Contraction(N,p,charge);project=projectors(casimir_matrix(N,charge),N,min(charge,N-charge))
        twiceJ=o.phis['hk'].toarray()
        rem=sum(np.kron(c.conj().T,c.T) for c in o.cs)
        identities_count=0;entries=0;maxentry=0
        for rank,(num,den) in enumerate(project):
            channel=filters(N,p,charge,0,0,rank)['channels'][2]
            alpha,e=channel['alpha'],channel['e']
            common=int(s.ilcm(alpha.denominator,e.denominator))
            aa=int(alpha*common);ee=int(e*common)
            for c in o.cs:
                S=np.kron(np.eye(o.d),c.T)
                T=np.kron(c,np.eye(o.dm))
                left=common*(num@S@twiceJ)
                right=2*(aa*(num@S)-ee*(num@rem@T))
                exact_integer(left);exact_integer(right)
                assert np.array_equal(left,right),(N,p,charge,rank)
                maxentry=max(maxentry,float(abs(left).max()),float(abs(right).max()))
                identities_count+=1;entries+=left.size
        rows.append(dict(N=N,p=p,m=charge,full_rectangular_map_identities=identities_count,
            complex_entry_equalities=entries,max_integer_entry=maxentry))
        print(rows[-1],flush=True)
    names=['code/check_charged_covariance_transfer.py','code/memory_two_charge_filter.py',
        'code/memory_quadratic_contraction.py','code/memory_rank_blocks.py']
    data=dict(status='passed',rows=rows,symbolic_projection_identities=3,
        elapsed_seconds=time.time()-start,
        scope='Complete linear maps for every external fermion index and every output rank; no sampling of rectangular matrices.',
        environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(data,indent=2)+'\n')

if __name__=='__main__':main()
