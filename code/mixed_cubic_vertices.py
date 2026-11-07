"""Exact multi-vertex DS polynomials for quadratic and cubic singlets.

Denominators: (1+q) per X (unnormalized Y), D per R; sigma per H slot.
Integer sparse polynomials have keys (lambda degree, q degree).
"""
from collections import Counter
from itertools import combinations
from pathlib import Path
import json,hashlib,platform,time
import sympy as s
from singlet_coefficient import pairings
from singlet_multimark import pmul,ppower,crossing,local
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-mixed-cubic'
q,lam=s.symbols('q lambda')

def poly(expr):
    return {tuple(map(int,k)):int(v) for k,v in s.Poly(s.expand(expr),lam,q).terms() if v}
def clean(a):return {k:v for k,v in a.items() if v}
def plus(*args):
    out=Counter()
    for a in args:out.update(a)
    return clean(out)
def scale(a,c):return {k:c*v for k,v in a.items() if c*v}
def expression(a):return s.factor(sum(v*lam**i*q**j for (i,j),v in a.items()))
D=(1+q*q)*(1+q+q*q);B=2*(1-q)*(1+q**3)
PR=[poly(-q*(1+q+q*q)),poly(-q*(1+q+q*q)),poly(D)]
PX=[poly(1+q)]
CX=poly(2*(1-q));CR=poly(B)
GX=poly(4*lam*(1+q)**2+16*lam*lam*q*(1+q))
GR=poly(lam*4*(1-q*q)*(1+q*q)*(1+q+q*q)**2
    +lam*lam*16*q*q*(1+q)*(1+q*q)*(1+q+q*q)*(2*q*q-q+2))

def blocks(word):
    vertices=[];owners=[];length=0
    for c in word:
        n={'H':1,'X':2,'R':3}[c];slots=list(range(length,length+n))
        owners.extend([len(vertices)]*n if c!='H' else [None])
        if c=='X':vertices.append(dict(kind=c,slots=slots,marks=[tuple(slots)],weights=PX,const=CX))
        elif c=='R':vertices.append(dict(kind=c,slots=slots,marks=[(slots[0],slots[1]),(slots[1],slots[2]),(slots[0],slots[2])],weights=PR,const=CR))
        length+=n
    return length,vertices,owners

def diagram(pa,vertices,verify=False):
    pi={a:b for a,b in pa}|{b:a for a,b in pa}
    cr=sum(crossing(a,b) for a,b in combinations(pa,2))
    alpha=[];single=[];pairs={};checks=0
    for v in vertices:
        al=[2-4*crossing(tuple(sorted((a,pi[a]))),tuple(sorted((b,pi[b])))) for a,b in v['marks']]
        alpha.append(al)
        single.append(pmul({(1,0):1},plus(*(scale(w,a) for w,a in zip(v['weights'],al)),scale(v['const'],-1))))
    for i,j in combinations(range(len(vertices)),2):
        vv=[]
        for k,((a,b),wa) in enumerate(zip(vertices[i]['marks'],vertices[i]['weights'])):
            for l,((c,d),wb) in enumerate(zip(vertices[j]['marks'],vertices[j]['weights'])):
                delta=4*int(pi[a]==d and pi[b]==c)
                if verify:
                    cr0,al,dd=local(2*len(pa),pa,[(a,b),(c,d)])
                    assert (cr0,al,dd[(0,1)])==(cr,[alpha[i][k],alpha[j][l]],delta)
                    checks+=1
                if delta:vv.append(scale(pmul(wa,wb),delta))
        pairs[i,j]=pmul({(1,0):1},plus(*vv))
    def rec(indices):
        if not indices:return {(0,0):1}
        i,*tail=indices;out=Counter(pmul(single[i],rec(tail)))
        for j in tail:
            if pairs.get((i,j)):out.update(pmul(pairs[i,j],rec([a for a in tail if a!=j])))
        return clean(out)
    return {(i,j+cr):v for (i,j),v in rec(list(range(len(vertices)))).items()},checks

def raw(word,verify=False):
    length,vertices,owners=blocks(word);total=Counter();count=checks=0
    if length%2:return {},0,0
    for pa in pairings(tuple(range(length))):
        if any(owners[a] is not None and owners[a]==owners[b] for a,b in pa):continue
        value,c=diagram(pa,vertices,verify);total.update(value);count+=1;checks+=c
    return clean(total),count,checks

def predicted(word):
    slots={c:tuple(i for i,a in enumerate(word) if a==c) for c in 'HXR'}
    if any(len(a)%2 for a in slots.values()):return {}
    out=Counter()
    for hh in pairings(slots['H']):
        for xx in pairings(slots['X']):
            for rr in pairings(slots['R']):
                chords=[(a,1) for a in hh]+[(a,2) for a in xx]+[(a,3) for a in rr]
                cr=sum(crossing(a,b)*i*j for (a,i),(b,j) in combinations(chords,2))
                out[0,cr]+=1
    return pmul(pmul(ppower(GX,len(slots['X'])//2),ppower(GR,len(slots['R'])//2)),dict(out))

def main():
    OUT.mkdir(parents=True,exist_ok=True);target=OUT/'moments.json'
    if target.exists():raise FileExistsError(target)
    rows=[];start=time.monotonic()
    for word in ['XX','RR','XRR','XXRR','XRXR','RRRR','XHRHR','XRRHH','XHRRH','HXRHR','XXRRHH','XHRXHR']:
        actual,n,checks=raw(word,verify=True);expected=predicted(word)
        diff=plus(actual,scale(expected,-1))
        row=dict(word=word,pairings=n,minor_checks=checks,raw_polynomial=str(expression(actual)),
            independent_probe_difference=str(expression(diff)),matches_independent_probe=(not diff),
            coefficients={f'{i},{j}':v for (i,j),v in actual.items()})
        rows.append(row);print(word,n,'difference',str(expression(diff)),round(time.monotonic()-start,2),flush=True)
    names=['code/mixed_cubic_vertices.py','code/singlet_multimark.py','code/singlet_coefficient.py']
    target.write_text(json.dumps(dict(status='computed',scope='Exact DS marked-minor polynomials, not finite-N matrix checks; q and lambda algebraically independent',
        records=rows,environment=dict(python=platform.python_version(),sympy=s.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')
if __name__=='__main__':main()
