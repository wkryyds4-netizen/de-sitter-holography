"""Normal and reverse two-slot derivative marks, with exact finite CAR checks."""
from pathlib import Path
from itertools import combinations
from fractions import Fraction as F
import hashlib,json,platform
import numpy as np
import sympy as s
from check_multimark_car import determinant
from singlet_multimark import nmul
from singlet_coefficient import pairings,dual_product
from complex_syk_dictionary import fermions,creators

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-oriented-kernel'


def graph_oriented(length,pa,marks,reverse):
    pi={a:b for a,b in pa}|{b:a for a,b in pa}
    def entry(a,b):
        values=[1 if a<=pi[b] else -1,0,0,0]
        for j,(x,y) in enumerate(marks):
            values[j+1]=(2 if reverse[j] else -2)*int((pi[b],a)==((y,x) if reverse[j] else (x,y)))
        return tuple(values)
    singles=[entry(i,i) for i in range(length)]
    assert all(v[1:]==(0,0,0) for v in singles)
    out=[-length*(length-1)//2-length//2,0,0,0]
    for i,j in combinations(range(length),2):
        aa=dual_product(entry(i,i),entry(j,j));bb=dual_product(entry(i,j),entry(j,i))
        den=singles[i][0]*singles[j][0]
        out=[v+(a-b)*den for v,a,b in zip(out,aa,bb)]
    cross=lambda a,b,c,d: min(a,b)<min(c,d)<max(a,b)<max(c,d) or min(c,d)<min(a,b)<max(c,d)<max(a,b)
    cr=sum(cross(a,b,c,d) for i,(a,b) in enumerate(pa) for c,d in pa[i+1:])
    assert out[0]==-4*cr
    return dict(crossings=cr,alpha=out[1],beta=out[2],delta=out[3])


def local_oriented(N,p,sets,pa,marks,reverse):
    pi={a:b for a,b in pa}|{b:a for a,b in pa};out={0:1}
    for mode in range(N):
        slots=[i for i in range(len(sets)) if mode in sets[i]];vals=[]
        for mask in range(1<<len(marks)):
            def ent(a,b):
                v=1 if a<=pi[b] else -1
                for j,(x,y) in enumerate(marks):
                    if mask>>j&1:v+=(2 if reverse[j] else -2)*int((pi[b],a)==((y,x) if reverse[j] else (x,y)))
                return v
            vals.append(determinant([[ent(a,b) for b in slots] for a in slots]))
        for j in range(len(marks)):
            for mask in range(1<<len(marks)):
                if mask>>j&1:vals[mask]-=vals[mask^(1<<j)]
        out=nmul(out,{k:v for k,v in enumerate(vals) if v})
    return [F((-1)**(len(sets)//2*p)*out.get(mask,0),2**(len(sets)*p)) for mask in range(1<<len(marks))]


def main():
    OUT.mkdir(exist_ok=True);target=OUT/'oriented-marks.json'
    if target.exists():raise FileExistsError(target)
    rules=0
    for length in [4,6,8]:
        for gap in range(length-3):
            marks=((0,1),(2+gap,3+gap))
            for pa in pairings(tuple(range(length))):
                if any(mark in pa for mark in marks):continue
                pi={a:b for a,b in pa}|{b:a for a,b in pa}
                for rev in [(False,False),(False,True),(True,False),(True,True)]:
                    g=graph_oriented(length,pa,marks,rev)
                    base=graph_oriented(length,pa,marks,(False,False))
                    assert (g['alpha'],g['beta'])==(base['alpha'],base['beta'])
                    (a,b),(c,d)=marks
                    delta=(4*int(pi[a]==d and pi[b]==c) if rev[0]==rev[1]
                           else -4*int(pi[a]==c and pi[b]==d))
                    assert g['delta']==delta
                    rules+=1
    rng=np.random.default_rng(110686);records=[];nonzero=negative=0
    for N,p,length in [(4,1,4),(4,2,4),(5,2,6),(6,2,8)]:
        cs=fermions(N);ss=list(combinations(range(N),p));bs=dict(zip(ss,creators(cs,p)))
        marks=((0,1),(length-2,length-1))
        allowed=[pa for pa in pairings(tuple(range(length))) if not any(mm in pa for mm in marks)]
        for case in range(40):
            pa=allowed[case%len(allowed)];pi={a:b for a,b in pa}|{b:a for a,b in pa}
            chosen=([ss[case%len(ss)]]*length if case<12 else [ss[rng.integers(len(ss))] for _ in range(length)])
            ops=[bs[chosen[i]]@bs[chosen[pi[i]]].conj().T for i in range(length)]
            for rev in [(False,False),(False,True),(True,False),(True,True)]:
                marked=[]
                for (a,b),rr in zip(marks,rev):
                    X,Y=ops[a],ops[b]
                    if rr:val=-sum((X@c-c@X)@(Y@c.conj().T-c.conj().T@Y) for c in cs)
                    else:val=-sum((X@c.conj().T-c.conj().T@X)@(Y@c-c@Y) for c in cs)
                    marked.append(val)
                actual=[]
                for mask in range(4):
                    word=np.eye(1<<N,dtype=complex);i=0
                    while i<length:
                        if i in [mm[0] for mm in marks]:
                            j=[mm[0] for mm in marks].index(i)
                            word=word@(marked[j] if mask>>j&1 else ops[i]@ops[i+1]);i+=2
                        else:word=word@ops[i];i+=1
                    tr=np.trace(word);assert tr.imag==0 and tr.real==round(tr.real)
                    actual.append(F(int(tr.real),1<<N))
                expected=local_oriented(N,p,chosen,pa,marks,rev)
                assert actual==expected,(N,p,pa,chosen,rev,actual,expected)
                nonzero+=int(actual[-1]!=0);negative+=int(actual[-1]<0)
                records.append(dict(N=N,p=p,length=length,pairing=pa,sets=chosen,reverse=rev,values=list(map(str,actual))))
    assert nonzero>100 and negative>0
    names=['code/memory_oriented_marks.py','code/check_multimark_car.py','code/singlet_multimark.py',
           'code/singlet_coefficient.py','code/complex_syk_dictionary.py']
    data=dict(status='passed',local_minor_rules=rules,finite_words=len(records),finite_comparisons=4*len(records),
              full_mark_nonzero=nonzero,full_mark_negative=negative,records=records,
              rule='Same alpha for both orientations. Equal orientations: delta=4 for reversed two bridges. Opposite orientations: delta=-4 for parallel two bridges.',
              environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__),
              input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print('passed',rules,'minor rules;',len(records),'CAR words;',nonzero,'nonzero;',negative,'negative')


if __name__=='__main__':main()
