"""Multiplicity-preserving rank blocks of the existing quadratic memory tensor.

No irreducible basis is chosen. Exact Casimir projectors act on both real
Gaussian input indices and the neutral output matrix. Blocks retain every
intertwiner multiplicity; integer Gram sums use bounded int64 chunks.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import reduce
from math import gcd,comb
import hashlib,json,platform,time
import numpy as np
import sympy as s
from complex_syk_dictionary import fermions
from singlet_wick_moments import gue_basis
from memory_quadratic_contraction import Contraction,anti

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-rank-blocks'


def exact_integer(x):
    assert np.max(abs(x))<2**48
    assert np.array_equal(x.real,np.rint(x.real))
    assert np.array_equal(x.imag,np.rint(x.imag))
    return x


def casimir_matrix(N,m):
    cs=fermions(N);ix=[a for a in range(1<<N) if a.bit_count()==m];d=len(ix);eye=np.eye(d)
    generators={}
    for i in range(N):
        for j in range(N):generators[i,j]=(cs[i].conj().T@cs[j])[np.ix_(ix,ix)]
    out=np.zeros((d*d,d*d),complex)
    ad=lambda e:np.kron(e,eye)-np.kron(eye,e.T)
    for i in range(N):
        for j in range(N):out+=ad(generators[i,j])@ad(generators[j,i])
    assert np.max(abs(out.imag))==0
    return exact_integer(out).real.astype(np.int64)


def projectors(C,N,maxrank):
    eye=np.eye(len(C),dtype=np.int64);roots=[2*j*(N+1-j) for j in range(maxrank+1)]
    out=[]
    for j,root in enumerate(roots):
        num=eye.copy();den=1
        for k,other in enumerate(roots):
            if k!=j:num=num@(C-other*eye);den*=root-other
        if den<0:num=-num;den=-den
        common=reduce(gcd,[den]+[abs(int(x)) for x in num.flat]);num//=common;den//=common
        assert np.array_equal(num@num,den*num)
        assert np.array_equal(C@num,root*num)
        dim=comb(N,j)**2-(comb(N,j-1)**2 if j else 0)
        assert np.trace(num)==den*dim
        out.append((num,den))
    for i,(a,d) in enumerate(out):
        for j,(b,e) in enumerate(out):
            if i!=j:assert not np.any(a@b)
    common=s.ilcm(*[d for _,d in out]);total=sum(a*int(common//d) for a,d in out)
    assert np.array_equal(total,int(common)*eye)
    return out


def tensor(N,p,m):
    o=Contraction(N,p,m);count=len(o.rs)
    prior=json.loads((ROOT/'results/susskind-memory-quadratic-contraction/finite-contractions.json').read_text())
    record=next(r for r in prior['rows'] if (r['N'],r['p'],r['m'])==(N,p,m));tbar=F(record['tbar'])
    assert tbar.denominator==1
    def values(h,k):
        f=o.tlinear(h,k)-tbar.numerator*h
        return np.array([h@h,o.A(h,k),anti(h,f)/2,o.contract(h,k)])
    diag=[values(h,k) for h,k in zip(o.hs,o.ks)]
    B=np.zeros((4,count,count,o.d*o.d),complex)
    for a in range(count):
        B[:,a,a]=8*diag[a].reshape(4,-1)
        for b in range(a+1,count):
            value=4*(values(o.hs[a]+o.hs[b],o.ks[a]+o.ks[b])-diag[a]-diag[b]).reshape(4,-1)
            B[:,a,b]=value;B[:,b,a]=value
    exact_integer(B)
    return o,B,record


def dot_integer(a,b,weight):
    ar=a.real.astype(np.int64).ravel();ai=a.imag.astype(np.int64).ravel()
    br=b.real.astype(np.int64).ravel();bi=b.imag.astype(np.int64).ravel()
    ww=weight.ravel();bound=max(1,int(max(abs(a).max(),abs(b).max())**2)*4)
    chunk=max(1,min(100000,(2**60)//bound));total=0
    assert bound<2**61
    for start in range(0,len(ar),chunk):
        sl=slice(start,start+chunk)
        total+=int(np.sum(ww[sl]*(ar[sl]*br[sl]+ai[sl]*bi[sl]),dtype=np.int64))
    return total


def main():
    OUT.mkdir(exist_ok=True);path=OUT/'rank-blocks.json'
    if path.exists():raise FileExistsError(path)
    N,p,m=5,2,3;start=time.time();o,B,parent=tensor(N,p,m);count=len(o.rs)
    # Convert the p-sector matrix Casimir to the original Hermitian GUE coordinates.
    _,basis,_=gue_basis(N,p);ix=[a for a in range(1<<N) if a.bit_count()==p]
    basis=basis[:,ix[:,None],ix].reshape(count,-1)
    cp=casimir_matrix(N,p)
    coord=(basis.conj()@cp@basis.T)*o.rs[:,None]/2
    assert not np.any(coord.imag);coord=exact_integer(coord).real.astype(np.int64)
    assert np.array_equal(coord,coord.T)
    assert np.array_equal(coord*o.rs[None,:],o.rs[:,None]*coord)
    pin=projectors(coord,N,p);pout=projectors(casimir_matrix(N,m),N,min(m,N-m))
    weight=np.broadcast_to(o.rs[:,None,None]*o.rs[None,:,None],(count,count,o.d*o.d))
    old=json.loads((ROOT/'results/susskind-memory-quadratic-contraction/scalar-removed-checks.json').read_text())['rows'][2]
    coeff=[-s.Rational(v) for v in old['projection_coefficients']]+[s.Integer(1)]
    vec=s.Matrix(coeff);rows=[];full=s.zeros(4);traceless=s.zeros(4)
    max_intermediate=0
    for a,(pa,da) in enumerate(pin):
        # Input projections use Pi^T on coefficient indices; here Pi is real symmetric.
        left=exact_integer(np.einsum('ca,icdz->iadz',pa,B,optimize=True))
        for b in range(a,len(pin)):
            pb,db=pin[b]
            both=exact_integer(np.einsum('iadz,db->iabz',left,pb,optimize=True))
            for r,(pr,dr) in enumerate(pout):
                projected=exact_integer(np.einsum('iabz,yz->iaby',both,pr,optimize=True))
                max_intermediate=max(max_intermediate,float(abs(projected).max()))
                scale=8*da*db*dr;factor=1 if a==b else 2
                G=s.zeros(4)
                for i in range(4):
                    for j in range(i,4):
                        val=s.Rational(factor*dot_integer(projected[i],projected[j],weight),2*scale*scale*o.d)
                        G[i,j]=G[j,i]=val
                full+=G
                if r:traceless+=G
                blocknorm=s.factor((vec.T*G*vec)[0])
                assert blocknorm>=0
                ids=list(G[:3,:3].rref()[1]);local_res=G[3,3]
                fits=[]
                if ids:
                    solve=G.extract(ids,ids).inv()*G.extract(ids,[3])
                    local_res=s.factor(G[3,3]-(G.extract(ids,[3]).T*solve)[0]);fits=list(map(str,solve))
                assert local_res>=0
                row=dict(input_ranks=[a,b],output_rank=r,gram_rank=G.rank(),
                         gram=[[str(x) for x in row] for row in G.tolist()],
                         global_residual_norm_squared=str(blocknorm),
                         old_pivot_basis=ids,local_projection_coefficients=fits,
                         local_residual_norm_squared=str(local_res))
                rows.append(row)
                print(json.dumps({k:row[k] for k in ['input_ranks','output_rank','gram_rank','global_residual_norm_squared','local_residual_norm_squared']}),flush=True)
    expected=s.Matrix([[s.Rational(x) for x in row] for row in parent['exact_gram']])
    expected_tf=s.Matrix([[s.Rational(x) for x in row] for row in old['traceless_gram']])
    assert full==expected and traceless==expected_tf
    residual=s.factor((vec.T*traceless*vec)[0]);assert residual==s.Rational(old['residual_norm_squared'])
    names=['code/memory_rank_blocks.py','code/memory_quadratic_contraction.py','code/singlet_wick_moments.py',
           'code/complex_syk_dictionary.py','results/susskind-memory-quadratic-contraction/finite-contractions.json',
           'results/susskind-memory-quadratic-contraction/scalar-removed-checks.json']
    result=dict(status='passed',N=N,p=p,m=m,covariance_w=1,rows=rows,
                basis=parent['basis'],sum_matches_full_and_traceless_previous_gram=True,
                traceless_global_residual_norm_squared=str(residual),
                input_projectors=[dict(rank=j,denominator=d,numerator=a.tolist()) for j,(a,d) in enumerate(pin)],
                output_projectors=[dict(rank=j,denominator=d,numerator=a.tolist()) for j,(a,d) in enumerate(pout)],
                max_integer_intermediate=max_intermediate,elapsed_seconds=time.time()-start,
                environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__),
                scope='Exact rank blocks preserving all multiplicities, on the previously studied finite system. No DS extrapolation.',
                input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    path.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
