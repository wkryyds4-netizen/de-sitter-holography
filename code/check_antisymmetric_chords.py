"""Independent numerical chord recursion and exact Laguerre checks.

Does not rerun or overwrite the exploratory exact-polynomial probe.
"""
from pathlib import Path
from functools import lru_cache
from math import comb
import json, hashlib, platform
import mpmath as mp
import sympy as sp
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-antisymmetric'

def case(hstr,dstr,prec):
    mp.mp.dps=prec
    h,d=mp.mpf(hstr),mp.mpf(dstr);q=mp.exp(-h)
    def qi(n):return -mp.expm1(-h*n)/(-mp.expm1(-h))
    def lower(v,side):
        w={}
        for (l,r),z in v.items():
            if l:w[l-1,r]=w.get((l-1,r),0)+z*qi(l)*(1 if side=='L' else q**(d+r))
            if r:w[l,r-1]=w.get((l,r-1),0)+z*qi(r)*(1 if side=='R' else q**(d+l))
        return w
    def add(v,w,c=1):
        z=v.copy()
        for k,a in w.items():z[k]=z.get(k,0)+c*a
        return z
    def raise_diff(v):
        z={}
        for (l,r),a in v.items():
            z[l+1,r]=z.get((l+1,r),0)+a
            z[l,r+1]=z.get((l,r+1),0)-a
        return z
    @lru_cache(None)
    def gram(l,r,L,R):
        if l+r!=L+R:return mp.mpf(0)
        v={(L,R):mp.mpf(1)}
        for _ in range(l):v=lower(v,'L')
        for _ in range(r):v=lower(v,'R')
        return v.get((0,0),mp.mpf(0))
    def dot(v,w):return mp.fsum(mp.conj(a)*b*gram(*i,*j) for i,a in v.items() for j,b in w.items())
    rows=[]
    for n in range(1,7):
        v={(l,n-l):mp.mpf((-1)**(n-l)*comb(n,l)) for l in range(n+1)}
        prev={(l,n-1-l):mp.mpf((-1)**(n-1-l)*comb(n-1,l)) for l in range(n)}
        norm=dot(v,v);N=mp.factorial(n)*mp.rf(2*d,n)
        dv=add(lower(v,'L'),lower(v,'R'),-1)
        res=add(dv,prev,-h*n*(2*d+n-1))
        rows.append(dict(n=n,norm_ratio=str(norm/(h**n*N)),aL_normalized_norm=str(mp.sqrt(dot(lower(v,'L'),lower(v,'L'))/norm)),lowering_error=str(mp.sqrt(abs(dot(res,res))/(h**(n+1)*N)))))
    # Direct Q_h moments: only the coefficient of grade zero contributes.
    def Q(v):
        z={k:(d+sum(k))*a for k,a in v.items()}
        z=add(z,add(lower(v,'L'),lower(v,'R'),-1),1/(2j*mp.sqrt(h)))
        return add(z,raise_diff(v),-1/(2j*mp.sqrt(h)))
    v={(0,0):mp.mpf(1)};moments=[]
    for k in range(1,9):
        v=Q(v);a=v.get((0,0),0);target=mp.rf(2*d,k)/2**k
        assert abs(mp.im(a))<mp.mpf('1e-65')
        moments.append(dict(k=k,value=str(mp.re(a)),limit=str(target),relative_error=str(abs(a/target-1))))
    # Whole-grade operator norms: independent of chosen descendants.
    norms=[]
    for n in range(5):
        G=mp.matrix(n+1);Gn=mp.matrix(n+2)
        for i in range(n+1):
            for j in range(n+1):G[i,j]=gram(i,n-i,j,n-j)
        for i in range(n+2):
            for j in range(n+2):Gn[i,j]=gram(i,n+1-i,j,n+1-j)
        W=mp.cholesky(G).T;Wn=mp.cholesky(Gn).T
        C=mp.matrix(n+2,n+1)
        for j in range(n+1):C[j,j]=-1;C[j+1,j]=1
        M=Wn*C*W**-1
        vals=mp.eigsy(M.T*M,eigvals_only=True)
        measured=vals[len(vals)-1]
        bound=2*sum(1-q**(d+j) for j in range(n+1))
        assert measured<=bound+mp.mpf('1e-60')
        norms.append(dict(n=n,norm_squared=str(measured),bound=str(bound)))
    return dict(h=hstr,d=dstr,dps=prec,descendants=rows,Q_moments=moments,whole_grade_bounds=norms)

def symbolic():
    p=sp.symbols('p',positive=True);d=sp.Integer(2)
    def phi(n):return sp.I**n*sp.sqrt(sp.factorial(n)/sp.rf(2*d,n))*2**d/sp.sqrt(sp.gamma(2*d))*p**(d-sp.Rational(1,2))*sp.exp(-p)*sp.assoc_laguerre(n,2*d-1,2*p)
    tests=[]
    for n in range(9):
        f=phi(n);up=sp.sqrt((n+1)*(2*d+n));down=sp.sqrt(n*(2*d+n-1))
        bm=(up*phi(n+1)+(down*phi(n-1) if n else 0))/2
        qm=(d+n)*f+sp.I*up*phi(n+1)/2-(sp.I*down*phi(n-1)/2 if n else 0)
        assert sp.simplify(sp.I*(p*sp.diff(f,p)+f/2)-bm)==0
        assert sp.simplify(p*f-qm)==0
        tests.append(dict(n=n,boost=True,Q=True))
    return tests

def main():
    target=OUT/'validation.json'
    if target.exists():raise FileExistsError(target)
    exact=symbolic();cases=[];refinement=[]
    for d in ['0.7','2']:
        for h in ['0.1','0.02','0.005']:
            a=case(h,d,80);b=case(h,d,110)
            errors=[]
            for ra,rb in zip(a['descendants'],b['descendants']):
                for k in ['norm_ratio','aL_normalized_norm','lowering_error']:
                    errors.append(abs(mp.mpf(ra[k])/mp.mpf(rb[k])-1))
            for ra,rb in zip(a['Q_moments'],b['Q_moments']):errors.append(abs(mp.mpf(ra['value'])/mp.mpf(rb['value'])-1))
            for ra,rb in zip(a['whole_grade_bounds'],b['whole_grade_bounds']):errors.append(abs(mp.mpf(ra['norm_squared'])/mp.mpf(rb['norm_squared'])-1))
            worst=max(errors);assert worst<mp.mpf('1e-50')
            cases.append(b);refinement.append(dict(h=h,d=d,max_relative_difference=str(worst)))
            print(d,h,'passed',flush=True)
    probe=json.loads((OUT/'probe.json').read_text())
    for row in probe['rows']:
        n=row['n'];assert row['norm_leading']==[n,str(sp.factorial(n)*sp.rf(4,n))]
    paths=[Path(__file__),ROOT/'code/antisymmetric_chord_probe.py',OUT/'probe.json']
    data=dict(status='passed',versions=dict(python=platform.python_version(),mpmath=mp.__version__,sympy=sp.__version__),symbolic_laguerre=exact,cases=cases,refinement=refinement,scope='Fixed grades only. Empirical precision refinement, not interval certification; all-n claims are proved in the manuscript.',input_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print('saved',target)
if __name__=='__main__':main()
