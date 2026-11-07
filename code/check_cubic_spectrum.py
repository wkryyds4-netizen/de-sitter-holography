"""Higher two-arc moments; compare, rather than assume, a q^3 probe kernel."""
from math import comb
from functools import lru_cache
import json,hashlib
import sympy as s
from cubic_singlet import ROOT,OUT,q,z,marks,moment
from singlet_coefficient import pairings

def probe(gap,end):
    total=0
    for pa in pairings(tuple(range(gap+end))):
        cross=sum(a<c<b<d or c<a<d<b for i,(a,b) in enumerate(pa) for c,d in pa[i+1:])
        bridge=sum(a<gap<=b for a,b in pa)
        total+=q**(cross+3*bridge)
    return s.expand(total)

def main():
    target=OUT/'spectrum-check.json'
    if target.exists():raise FileExistsError(target)
    data=json.loads((OUT/'study.json').read_text())
    parse=lambda v:s.sympify(v.replace('lambda','lam'),locals={'q':q,'lam':z})
    v=s.Matrix([1]+[-parse(c) for c in data['projection']]);V=parse(data['residual_variance'])
    names=['T','S','Z'];rows=[];moments={};crosses=[]
    for n in [0,2,4]:
        full=0
        for k in range(n+1):
            M=s.Matrix(3,3,lambda i,j:moment(names[i],names[j],n-k,k)[0])
            val=s.factor((v.T*M*v)[0]);expected=s.factor(V*probe(n-k,k))
            diff=s.factor(val-expected)
            rows.append(dict(order=n,gap=n-k,end=k,value=str(val),q3_difference=str(diff),matches=diff==0))
            full+=(-1)**k*comb(n,k)*val
            crosses.append(dict(order=n,gap=n-k,end=k,R_S=str(s.factor((v.T*M)[1])),R_Z=str(s.factor((v.T*M)[2]))))
        full=s.factor(full/V/z**(n//2))
        moments[n]=dict(normalized_over_Js_power=str(full),small_lambda=str(s.limit(full.subs(q,s.exp(-4*z)),z,0)))
        print(n,moments[n],'matches',[r['matches'] for r in rows if r['order']==n],flush=True)
    out=dict(status='computed',scope='exact fixed-degree DS two-arc sums; no all-order equivalence claimed',rows=rows,cross_moments=crosses,moments=moments,
      input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['code/check_cubic_spectrum.py','code/cubic_singlet.py','results/susskind-cubic-singlet/study.json']})
    target.write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':main()
