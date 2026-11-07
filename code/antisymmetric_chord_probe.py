"""Exact exploratory Gram polynomials for antisymmetric chord descendants."""
from functools import lru_cache
from math import comb
import sympy as s
q=s.symbols('q');d=2
@lru_cache(None)
def qi(n):return sum(q**i for i in range(n))
def lower(v,side):
    out={}
    for (l,r),c in v.items():
        if l:
            coef=qi(l)*(1 if side=='L' else q**(d+r))
            out[l-1,r]=out.get((l-1,r),0)+c*coef
        if r:
            coef=qi(r)*(1 if side=='R' else q**(d+l))
            out[l,r-1]=out.get((l,r-1),0)+c*coef
    return {k:s.expand(v) for k,v in out.items()}
@lru_cache(None)
def gram(l,r,L,R):
    if l+r!=L+R:return s.Integer(0)
    v={(L,R):s.Integer(1)}
    for j in range(l):v=lower(v,'L')
    for j in range(r):v=lower(v,'R')
    return s.expand(v.get((0,0),0))
def vector(n):return {(l,n-l):s.Integer((-1)**(n-l)*comb(n,l)) for l in range(n+1)}
def dot(v,w):return s.expand(sum(a*b*gram(*k,*j) for k,a in v.items() for j,b in w.items()))
def difference(v,w):return {k:s.expand(v.get(k,0)-w.get(k,0)) for k in v.keys()|w.keys()}
def lead(z):
    p=s.Poly(z,q)
    if p.is_zero:return None,'0'
    for k in range(p.degree()+1):
        val=sum(c*s.prod(range(e[0]-k+1,e[0]+1)) for e,c in p.terms() if e[0]>=k)
        if val:return k,str((-1)**k*val/s.factorial(k))
if __name__=='__main__':
    import json
    from pathlib import Path
    out=Path(__file__).resolve().parent.parent/'results/susskind-antisymmetric'
    out.mkdir(exist_ok=True,parents=True)
    target=out/'probe.json'
    if target.exists():raise FileExistsError(target)
    rows=[]
    for n in range(1,7):
        v=vector(n);norm=dot(v,v);a=lower(v,'L');b=lower(v,'R')
        dv=difference(a,b);prev=vector(n-1)
        residual=difference(dv,{k:(1-q)*n*(2*d+n-1)*c for k,c in prev.items()})
        row=dict(n=n,norm_leading=lead(norm),aL_norm_leading=lead(dot(a,a)),lowering_residual_leading=lead(dot(residual,residual)),norm_polynomial=str(norm))
        rows.append(row);print({k:v for k,v in row.items() if k!='norm_polynomial'},flush=True)
    target.write_text(json.dumps(dict(d=d,rows=rows),indent=2)+'\n')
